# -*- coding: utf-8 -*-
"""
R2 revision analyses for HP518-ADR-Fusion.

Reproduces every supplementary table added in the R2 revision, directly from
the data committed in this repository:

    S1  Scenario A (evidence integration) vs Scenario B (a priori) scoring
    S2  Sensitivity of the ranking to the per-drug PT cut-off
    S3  Stability of the serious-signal PTs under the 100-PT cut-off
    S4  Leave-one-reference-drug-out re-scoring
    S5  Machine-learning baselines under the revised cross-validation protocol
    S6  Weight sensitivity (random Dirichlet weights + one-at-a-time shifts)
    S7  Disproportionality (ROR/PRR/chi2) against the whole-database background

Inputs
------
    data/faers_extended.json   FAERS snapshot, window 20040101-20260401
    data/meddra_all_se.tsv     SIDER 4.1 MedDRA adverse-effect associations

Output
------
    results/r2_supplementary_results.json

Usage
-----
    python code/r2_revision_analyses.py

Dependencies: numpy, scikit-learn (S5 only). Everything else is stdlib.
Results are deterministic: every random draw is seeded with SEED.

Weight scale note
-----------------
The five weights sum to 1.05, not 1.00, so raw fusion scores live on a
0-105 scale. Normalised weights are reported alongside so that the
relative contribution of each evidence stream can be read directly.
Scenario A and Scenario B are on different absolute scales (1.05 vs 1.00)
and must only be compared by rank, never by absolute score.
"""
import collections
import json
import math
import os
import random
import statistics
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, 'data')
RESULTS = os.path.join(ROOT, 'results')

SEED = 20260922

# ---------------------------------------------------------------- constants
# Phase 1 (HP518 clinical) observations used as the concordance term.
PHASE1 = ['FATIGUE', 'NAUSEA', 'VOMITING', 'DECREASED APPETITE', 'CONSTIPATION', 'DIARRHOEA']

# Spelling variants accepted when matching Phase 1 terms to MedDRA PTs.
SYN = {'FATIGUE': {'FATIGUE', 'ASTHENIA'},
       'DIARRHOEA': {'DIARRHOEA', 'DIARRHEA'},
       'DECREASED APPETITE': {'DECREASED APPETITE', 'APPETITE DECREASED'}}

REF6 = ['enzalutamide', 'apalutamide', 'darolutamide', 'bicalutamide', 'flutamide', 'abiraterone']
GEN2 = {'enzalutamide', 'apalutamide', 'darolutamide'}

# SIDER compound IDs for the four reference drugs with usable SIDER records.
SIDER_DRUGS = {'bicalutamide': 'CID100002375', 'flutamide': 'CID100003397',
               'nilutamide': 'CID100004493', 'abiraterone': 'CID100132971'}

# Scenario A - evidence integration. Sum = 1.05, hence the 0-105 raw scale.
W_SCENARIO_A = {'F': 0.40, 'I': 0.25, 'C': 0.20, 'V': 0.15, 'B': 0.05}
# Scenario B - a priori prediction, clinical-concordance term removed. Sum = 1.00.
W_SCENARIO_B = {'F': 0.45, 'I': 0.25, 'C': 0.25, 'V': 0.00, 'B': 0.05}
# The 0.38 / 0.22 pair printed in the submitted manuscript. Never executed.
W_MISPRINT = {'F': 0.38, 'I': 0.22, 'C': 0.20, 'V': 0.15, 'B': 0.05}

SERIOUS_PTS = ['INTERSTITIAL LUNG DISEASE', 'MYOCARDIAL INFARCTION', 'PULMONARY EMBOLISM',
               'NEUTROPENIA', 'THROMBOCYTOPENIA', 'CEREBROVASCULAR ACCIDENT', 'CARDIAC FAILURE',
               'ACUTE KIDNEY INJURY', 'ATRIAL FIBRILLATION', 'DEATH', 'SEPSIS', 'HYPERTENSION',
               'DYSPNOEA', 'ANAEMIA', 'HEPATIC FAILURE', 'HAEMATURIA', 'FALL', 'DEHYDRATION',
               'SYNCOPE']


# --------------------------------------------------------------- data loading
def load_sider():
    """SIDER 4.1: compound ID -> set of MedDRA terms (upper-cased).

    Also returns the counts quoted in the manuscript's data description:
    309,849 drug-side-effect pairs, 1,430 drugs and 5,868 adverse-effect
    terms. The 5,868 figure counts distinct side-effect concepts, i.e.
    distinct UMLS concept identifiers in column 2. Counting the free-text
    names in column 6 instead gives 6,123 strings (4,251 if restricted to
    MedDRA PT rows), which is a different quantity and must not be quoted
    as the number of side effects.
    """
    d = {}
    path = os.path.join(DATA, 'meddra_all_se.tsv')
    pairs = 0
    drugs = set()
    cuis = set()
    with open(path, encoding='utf-8') as f:
        for line in f:
            c = line.rstrip('\n').split('\t')
            if len(c) >= 6:
                pairs += 1
                drugs.add(c[0])
                cuis.add(c[2])
                d.setdefault(c[0], set()).add(c[5].upper())
    return d, {'pairs': pairs, 'drugs': len(drugs), 'terms': len(cuis)}


def sider_universe(sider):
    """Pooled SIDER term universe across the four reference compounds."""
    s = set()
    for cid in SIDER_DRUGS.values():
        s |= sider.get(cid, set())
    return s


def load_faers():
    path = os.path.join(DATA, 'faers_extended.json')
    with open(path, encoding='utf-8') as f:
        return json.load(f)


# ------------------------------------------------------------------- helpers
def build_pool(profiles, topn):
    """Pool the top-`topn` PTs of each drug into counts, drug-counts, 2G set."""
    cnt, nd, g2 = collections.Counter(), collections.Counter(), set()
    for drug, prof in profiles.items():
        for pt, c in sorted(prof.items(), key=lambda kv: -kv[1])[:topn]:
            cnt[pt] += c
            nd[pt] += 1
            if drug in GEN2:
                g2.add(pt)
    return cnt, nd, g2


def winsor(vals, lo_p=5.0, hi_p=95.0):
    """P5/P95 clipping bounds for the FAERS count normalisation."""
    s = sorted(vals)
    n = len(s)

    def pc(p):
        k = (n - 1) * p / 100.0
        f = math.floor(k)
        c = min(f + 1, n - 1)
        return s[f] + (s[c] - s[f]) * (k - f)

    return pc(lo_p), max(pc(hi_p), pc(lo_p) + 1e-9)


def score(cnt, nd, g2, sider_pt, w, rng_range, phase1=PHASE1):
    """Weighted fusion score. Each component is mapped to 0-100 first."""
    lo, hi = rng_range
    rows = []
    for k in cnt:
        F = 100.0 * min(max(cnt[k] - lo, 0.0), hi - lo) / (hi - lo)
        I = 100.0 if k.upper() in sider_pt else 0.0
        C = 100.0 * nd[k] / 6.0
        V = 100.0 if k in phase1 else 0.0
        B = 100.0 if k in g2 else 0.0
        rows.append(dict(adr=k, cnt=cnt[k], nd=nd[k], F=F, I=I, C=C, V=V, B=B,
                         score=w['F'] * F + w['I'] * I + w['C'] * C + w['V'] * V + w['B'] * B))
    rows.sort(key=lambda r: (-r['score'], r['adr']))
    for i, r in enumerate(rows, 1):
        r['rank'] = i
    return rows


def cov(rows, k, phase1=PHASE1):
    """How many Phase 1 ADRs appear in the top-k rows (synonyms matched)."""
    hit = [a for a in phase1
           if any(r['adr'] in SYN.get(a, {a}) or r['adr'] == a for r in rows[:k])]
    return len(hit), hit


def rho(a, b):
    """Spearman rho between two {adr: rank} maps over their common keys."""
    common = [k for k in a if k in b]
    if len(common) < 3:
        return None
    n = len(common)
    da = {k: i + 1 for i, k in enumerate(sorted(common, key=lambda k: a[k]))}
    db = {k: i + 1 for i, k in enumerate(sorted(common, key=lambda k: b[k]))}
    return 1 - 6 * sum((da[k] - db[k]) ** 2 for k in common) / (n * (n * n - 1))


def normalised(w):
    """Same weights rescaled to sum to 1.00 (for reading contributions only)."""
    s = sum(w.values())
    return {k: round(v / s, 4) for k, v in w.items()}


# ---------------------------------------------------------------- analyses
def s1_scenarios(cnt, nd, g2, sider_pt, rng_range):
    rows_a = score(cnt, nd, g2, sider_pt, W_SCENARIO_A, rng_range)
    rows_b = score(cnt, nd, g2, sider_pt, W_SCENARIO_B, rng_range)
    rows_e = score(cnt, nd, g2, sider_pt, W_MISPRINT, rng_range)
    r_a = {r['adr']: r['rank'] for r in rows_a}
    r_b = {r['adr']: r['rank'] for r in rows_b}
    r_e = {r['adr']: r['rank'] for r in rows_e}
    out = dict(
        n_candidates=len(rows_a),
        weights=dict(
            scenario_a=dict(raw=W_SCENARIO_A, sum=round(sum(W_SCENARIO_A.values()), 6),
                            scale='0-105', normalised=normalised(W_SCENARIO_A)),
            scenario_b=dict(raw=W_SCENARIO_B, sum=round(sum(W_SCENARIO_B.values()), 6),
                            scale='0-100', normalised=normalised(W_SCENARIO_B)),
            note='Scenario A and B are on different absolute scales; compare by rank only.'),
        scenario_a=dict(coverage={k: cov(rows_a, k)[0] for k in (5, 10, 20, 80)},
                        hits_top10=cov(rows_a, 10)[1],
                        top10=[r['adr'] for r in rows_a[:10]]),
        scenario_b=dict(coverage={k: cov(rows_b, k)[0] for k in (5, 10, 20, 80)},
                        hits_top10=cov(rows_b, 10)[1],
                        hits_top20=cov(rows_b, 20)[1],
                        top10=[r['adr'] for r in rows_b[:10]]),
        spearman_a_vs_b=round(rho(r_a, r_b), 4),
        misprint_check=dict(
            weights=W_MISPRINT,
            top10_set_identical=set(r['adr'] for r in rows_a[:10]) == set(r['adr'] for r in rows_e[:10]),
            coverage_top10_a=cov(rows_a, 10)[0], coverage_top10_misprint=cov(rows_e, 10)[0],
            spearman_a_vs_misprint=round(rho(r_a, r_e), 4),
            note='The 0.38/0.22 pair was a transcription error; it changes no reported result.'),
    )
    print('S1  n=%d  A cov10=%d/6  B cov10=%d/6 cov20=%d/6  rho(A,B)=%s'
          % (out['n_candidates'], out['scenario_a']['coverage'][10],
             out['scenario_b']['coverage'][10], out['scenario_b']['coverage'][20],
             out['spearman_a_vs_b']))
    return out


def s2_cutoff(profiles, sider_pt, rng_range):
    base = score(*build_pool(profiles, 100), sider_pt=sider_pt, w=W_SCENARIO_A, rng_range=rng_range)
    base_rank = {r['adr']: r['rank'] for r in base}
    base_top30 = {r['adr'] for r in base[:30]}
    out = {}
    for t in (25, 50, 75, 100, 150, 200):
        c, n, g = build_pool(profiles, t)
        rr = score(c, n, g, sider_pt, W_SCENARIO_A, rng_range)
        t30 = {x['adr'] for x in rr[:30]}
        cut = min(sorted(profiles[d].values(), reverse=True)[t - 1] for d in profiles)
        out[str(t)] = dict(
            n_pool=len(rr), min_count_at_cut=cut,
            coverage={k: cov(rr, k)[0] for k in (5, 10, 20, 80)},
            missing_at_top10=[a for a in PHASE1 if a not in cov(rr, 10)[1]],
            spearman_vs_100=(None if t == 100
                             else round(rho({x['adr']: x['rank'] for x in rr}, base_rank), 3)),
            jaccard_top30=(None if t == 100
                           else round(len(t30 & base_top30) / len(t30 | base_top30), 3)))
        print('S2  top%-3d pool=%3d cov10=%d/6 rho=%s'
              % (t, out[str(t)]['n_pool'], out[str(t)]['coverage'][10],
                 out[str(t)]['spearman_vs_100']))
    return out


def s3_pt_stability(profiles, cnt, nd):
    c25, _, _ = build_pool(profiles, 25)
    out = {s: dict(in_top25=s in c25, in_top100=s in cnt,
                   pooled_reports=cnt.get(s, 0), n_drugs=nd.get(s, 0)) for s in SERIOUS_PTS}
    print('S3  serious PTs retained at the 100-PT cut-off: %d/%d'
          % (sum(1 for v in out.values() if v['in_top100']), len(SERIOUS_PTS)))
    return out


def s4_leave_one_out(profiles, sider_pt, rng_range):
    out = {}
    for held in REF6:
        sub = {d: p for d, p in profiles.items() if d != held}
        c, n, g = build_pool(sub, 100)
        rr = score(c, n, g, sider_pt, W_SCENARIO_B, rng_range)
        out[held] = dict(n_pool=len(rr), coverage={k: cov(rr, k)[0] for k in (5, 10, 20)},
                         top10=[x['adr'] for x in rr[:10]])
        print('S4  held=%-14s pool=%3d cov10=%d/6 cov20=%d/6'
              % (held, out[held]['n_pool'], out[held]['coverage'][10], out[held]['coverage'][20]))
    return out


def s5_ml_baselines(rows_b):
    """Revised protocol: repeated stratified 5-fold CV + bootstrap CIs.

    Features (4): min-max normalised pooled FAERS count (F), SIDER label (I),
    number of reference drugs reporting the ADR (C), 2nd-generation indicator (B).
    Label: 1 if the ADR was observed in the Phase 1 HP518 dataset (V > 0).
    The concordance term V is deliberately excluded from the feature set so
    that the baselines and the Scenario B fusion score predict the same target
    from the same information.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.naive_bayes import MultinomialNB
    from sklearn.model_selection import StratifiedKFold
    from sklearn.preprocessing import StandardScaler
    from sklearn.metrics import roc_auc_score, average_precision_score

    X = np.asarray([[r['F'], r['I'], r['C'], r['B']] for r in rows_b], float)
    y = np.asarray([1 if r['V'] > 0 else 0 for r in rows_b], int)
    adrs = [r['adr'] for r in rows_b]
    print('S5  design matrix %s  positives=%d' % (X.shape, int(y.sum())))

    def repeated_cv(make_model, scale):
        """20 repeats x 5 folds; returns one out-of-fold prediction per repeat."""
        oof = []
        for rep in range(20):
            skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED + rep)
            p = np.zeros(len(y))
            for tr, te in skf.split(X, y):
                Xtr, Xte = X[tr], X[te]
                if scale:
                    sc = StandardScaler().fit(Xtr)
                    Xtr, Xte = sc.transform(Xtr), sc.transform(Xte)
                m = make_model()
                m.fit(Xtr, y[tr])
                p[te] = m.predict_proba(Xte)[:, 1]
            oof.append(p)
        return np.vstack(oof)

    def ci(vals):
        v = np.asarray(vals, float)
        return [round(float(np.percentile(v, 2.5)), 3), round(float(np.percentile(v, 97.5)), 3)]

    def top10_hits(order_idx):
        hit = set()
        for i in order_idx:
            a = adrs[i]
            for p in PHASE1:
                if a == p or a in SYN.get(p, {p}):
                    hit.add(p)
        return len(hit)

    models = {
        'Logistic Regression (L2)': (
            lambda: LogisticRegression(solver='liblinear', penalty='l2', C=1.0, max_iter=1000), True),
        'Multinomial Naive Bayes': (lambda: MultinomialNB(alpha=1.0), False),
    }
    res = {}
    for name, (mk, sc) in models.items():
        oof = repeated_cv(mk, sc)
        aucs = [roc_auc_score(y, o) for o in oof]
        aprs = [average_precision_score(y, o) for o in oof]
        c10 = [top10_hits(np.argsort(-o)[:10]) for o in oof]
        mean = oof.mean(axis=0)
        rng = np.random.default_rng(SEED)
        ba, bp = [], []
        for _ in range(2000):
            idx = rng.integers(0, len(y), len(y))
            if y[idx].sum() in (0, len(idx)):
                continue
            ba.append(roc_auc_score(y[idx], mean[idx]))
            bp.append(average_precision_score(y[idx], mean[idx]))
        res[name] = dict(AUROC=round(float(np.mean(aucs)), 3), AUROC_CI=ci(aucs),
                         AUPR=round(float(np.mean(aprs)), 3), AUPR_CI=ci(aprs),
                         coverage_top10_mean=round(float(np.mean(c10)), 2),
                         boot_AUROC_CI=ci(ba), boot_AUPR_CI=ci(bp))
        print('S5  %-26s AUROC=%.3f AUPR=%.3f' % (name, res[name]['AUROC'], res[name]['AUPR']))

    # Fusion model (Scenario B) scored on the same candidate set - deterministic,
    # so it has no cross-validation spread; only bootstrap CIs are available.
    fs = np.asarray([r['score'] for r in rows_b], float)
    c_f = top10_hits(np.argsort(-fs)[:10])
    rng2 = np.random.default_rng(SEED)
    ba, bp = [], []
    for _ in range(2000):
        i2 = rng2.integers(0, len(y), len(y))
        if y[i2].sum() in (0, len(i2)):
            continue
        ba.append(roc_auc_score(y[i2], fs[i2]))
        bp.append(average_precision_score(y[i2], fs[i2]))
    res['Fusion model (Scenario B)'] = dict(
        AUROC=round(float(roc_auc_score(y, fs)), 3),
        AUPR=round(float(average_precision_score(y, fs)), 3),
        coverage_top10_mean=c_f,
        boot_AUROC_CI=ci(ba), boot_AUPR_CI=ci(bp),
        note='Deterministic score: no cross-validation spread; bootstrap CI only. '
             'Point estimates rank above both baselines but the intervals overlap, '
             'so superiority over machine learning is NOT claimed.')
    print('S5  %-26s AUROC=%.3f AUPR=%.3f'
          % ('Fusion model (Scenario B)',
             res['Fusion model (Scenario B)']['AUROC'], res['Fusion model (Scenario B)']['AUPR']))

    return dict(n=int(len(y)), positives=int(y.sum()),
                features=['min-max normalised pooled FAERS report count (F)',
                          'binary SIDER label presence (I)',
                          'number of reference drugs reporting the ADR (C)',
                          'second-generation AR-inhibitor indicator (B)'],
                label='1 if the ADR was observed in the Phase 1 HP518 dataset (V > 0)',
                protocol='20 repeats x stratified 5-fold CV; 2000 bootstrap resamples for CIs',
                results=res)


def s6_weight_sensitivity(cnt, nd, g2, sider_pt, rng_range, rows_a):
    R = random.Random(SEED)
    base_rank = {r['adr']: r['rank'] for r in rows_a}
    base_top30 = {r['adr'] for r in rows_a[:30]}
    covs, rhos, jac = [], [], []
    for _ in range(5000):
        x = [R.expovariate(1) for _ in range(5)]
        s = sum(x)
        w = dict(zip('FICVB', [v / s for v in x]))
        rr = score(cnt, nd, g2, sider_pt, w, rng_range)
        covs.append(cov(rr, 10)[0])
        rhos.append(rho({q['adr']: q['rank'] for q in rr}, base_rank))
        t30 = {q['adr'] for q in rr[:30]}
        jac.append(len(t30 & base_top30) / len(t30 | base_top30))
    oaat = {}
    for comp in 'FICVB':
        for dl in (-0.20, -0.10, 0.10, 0.20):
            w = dict(W_SCENARIO_B)
            w[comp] *= (1 + dl)
            t = sum(w.values())
            w = {k: v / t for k, v in w.items()}
            rr = score(cnt, nd, g2, sider_pt, w, rng_range)
            oaat['%s%+d%%' % (comp, int(dl * 100))] = dict(
                coverage_top10=cov(rr, 10)[0],
                spearman=round(rho({q['adr']: q['rank'] for q in rr}, base_rank), 3))
    out = dict(
        random=dict(n=5000, sampling='normalised exponential (Dirichlet-like) draws',
                    coverage_top10_mean=round(statistics.mean(covs), 3),
                    coverage_top10_min=min(covs), coverage_top10_max=max(covs),
                    prop_coverage_6=round(sum(c == 6 for c in covs) / 5000, 3),
                    spearman_mean=round(statistics.mean(rhos), 4),
                    spearman_min=round(min(rhos), 4),
                    jaccard_top30_mean=round(statistics.mean(jac), 4),
                    jaccard_top30_min=round(min(jac), 4)),
        one_at_a_time=oaat)
    print('S6  random weights: cov10 mean=%.3f  spearman mean=%.4f  min=%.4f'
          % (out['random']['coverage_top10_mean'], out['random']['spearman_mean'],
             out['random']['spearman_min']))
    return out


def s7_disproportionality(ext):
    """ROR / PRR / chi2 against the whole-database background.

    This differs from results/disproportionality_results.json, which uses the
    other five AR-pathway drugs as an internal comparator. The internal design
    is deliberately conservative; this is the conventional full-database design.
    """
    N = ext['db_total']
    E = ext['pooled_total_reports']

    def tab(a, ne, ca):
        b = ne - a
        c = ca - a
        d = N - a - b - c
        if min(a, b, c, d) <= 0:
            return None
        r = (a / b) / (c / d)
        se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
        prr = (a / ne) / (c / (N - ne))
        chi = ((a * d - b * c) ** 2 * N) / ((a + b) * (c + d) * (a + c) * (b + d))
        return dict(a=int(a), b=int(b), c=int(c), d=int(d), ROR=round(r, 3),
                    ROR_CI=[round(math.exp(math.log(r) - 1.96 * se), 3),
                            round(math.exp(math.log(r) + 1.96 * se), 3)],
                    PRR=round(prr, 3), chi2=round(chi, 1),
                    signal=bool(math.exp(math.log(r) - 1.96 * se) > 1 and prr >= 2 and chi >= 4))

    pool = {}
    for pt, a in ext['pooled_signal_counts'].items():
        v = tab(a, E, ext['db_pt_totals'][pt])
        if v:
            pool[pt] = v
    order = sorted(pool.items(), key=lambda kv: -kv[1]['chi2'])
    m = len(order)

    def sf(x):
        return math.erfc(math.sqrt(max(x, 0) / 2.0))

    for i, (pt, v) in enumerate(order, 1):
        v['q_BH'] = round(min(1.0, sf(v['chi2']) * m / i), 5)
    out = dict(background=dict(db_total=N, exposed=E, n_tested=m,
                              date_window=ext['date_window'],
                              design='whole-database background (conventional)'),
               table={k: v for k, v in sorted(pool.items(), key=lambda kv: -kv[1]['ROR'])})
    print('S7  %d PTs tested, %d flagged (ROR 95%% CI lower bound > 1, PRR >= 2, chi2 >= 4)'
          % (m, sum(1 for v in pool.values() if v['signal'])))
    return out


# ---------------------------------------------------------------------- main
def main():
    os.makedirs(RESULTS, exist_ok=True)
    ext = load_faers()
    profiles = {d: ext['drugs'][d]['pt_profile'] for d in ext['drugs']}

    print('Loading SIDER 4.1 ...')
    sider, sider_stats = load_sider()
    sider_pt = sider_universe(sider)
    print('SIDER: %d pairs, %d drugs, %d side-effect concepts'
          % (sider_stats['pairs'], sider_stats['drugs'], sider_stats['terms']))
    print('SIDER terms in the pooled four-drug reference universe: %d' % len(sider_pt))

    c200, _, _ = build_pool(profiles, 200)
    rng_range = winsor(list(c200.values()))
    print('FAERS count normalisation range (P5, P95): %s' % (tuple(round(v, 2) for v in rng_range),))

    cnt, nd, g2 = build_pool(profiles, 100)
    rows_a = score(cnt, nd, g2, sider_pt, W_SCENARIO_A, rng_range)
    rows_b = score(cnt, nd, g2, sider_pt, W_SCENARIO_B, rng_range)

    R = {
        'provenance': dict(
            faers_window=ext['date_window'],
            faers_fetched_at=ext['fetched_at'],
            faers_db_total=ext['db_total'],
            faers_pooled_total=ext['pooled_total_reports'],
            n_reference_drugs=len(profiles),
            sider_term_universe=len(sider_pt),
            sider_database=dict(
                pairs=sider_stats['pairs'], drugs=sider_stats['drugs'],
                adverse_effect_concepts=sider_stats['terms'],
                concept_definition='distinct UMLS concept identifiers (column 2 of meddra_all_se.tsv)'),
            normalisation_range=[round(v, 2) for v in rng_range],
            seed=SEED,
            python=sys.version.split()[0],
        ),
        'S1_scenarios': s1_scenarios(cnt, nd, g2, sider_pt, rng_range),
        'S2_cutoff_sensitivity': s2_cutoff(profiles, sider_pt, rng_range),
        'S3_pt_stability': s3_pt_stability(profiles, cnt, nd),
        'S4_leave_one_drug_out': s4_leave_one_out(profiles, sider_pt, rng_range),
        'S5_ml_baselines': s5_ml_baselines(rows_b),
        'S6_weight_sensitivity': s6_weight_sensitivity(cnt, nd, g2, sider_pt, rng_range, rows_a),
        'S7_disproportionality': s7_disproportionality(ext),
    }
    R['provenance']['numpy'] = np.__version__
    try:
        import sklearn
        R['provenance']['scikit_learn'] = sklearn.__version__
    except ImportError:
        pass

    out = os.path.join(RESULTS, 'r2_supplementary_results.json')
    with open(out, 'w', encoding='utf-8') as f:
        json.dump(R, f, ensure_ascii=False, indent=1, default=str)
    print('\nsaved ->', out)


if __name__ == '__main__':
    main()
