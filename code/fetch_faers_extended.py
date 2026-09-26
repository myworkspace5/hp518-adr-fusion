# -*- coding: utf-8 -*-
"""
Fetch the extended FAERS snapshot used by the R2 revision analyses.

Output: data/faers_extended.json

Contents
--------
  db_total             total FAERS case reports inside the date window
  drugs                per-drug report total + top-200 MedDRA PT profile
  pooled_total_reports reports for the six-drug AR-pathway pool
  pooled_signal_counts reports for each serious-signal PT within the pool
  db_pt_totals         reports for each serious-signal PT in the whole database
  or_check             sanity check that the OR syntax is additive

The date window is pinned to 20040101-20260401 so that re-running the script
reproduces the snapshot reported in the manuscript (20,344,365 reports).
Changing DATE produces a different snapshot and will change every downstream
number; do not change it without re-running the full analysis chain.

Usage
-----
    python code/fetch_faers_extended.py            # writes data/faers_extended.json
    python code/fetch_faers_extended.py --out x.json

Requires network access to https://api.fda.gov (no authentication).
Behind a TLS-inspecting proxy, export REQUESTS_CA_BUNDLE / SSL_CERT_FILE
pointing at the proxy CA instead of disabling certificate verification.
"""
import argparse
import json
import os
import sys
import time
import urllib.error
import urllib.request

BASE = 'https://api.fda.gov/drug/event.json'
DATE = 'receivedate:[20040101+TO+20260401]'

DRUGS = {
    'enzalutamide': 'ENZALUTAMIDE',
    'apalutamide': 'APALUTAMIDE',
    'darolutamide': 'DAROLUTAMIDE',
    'bicalutamide': 'BICALUTAMIDE',
    'flutamide': 'FLUTAMIDE',
    'abiraterone': 'ABIRATERONE',
}

# Serious-signal preferred terms (MedDRA PT level) tracked in the R2 revision.
SIGNALS = [
    'INTERSTITIAL LUNG DISEASE', 'MYOCARDIAL INFARCTION', 'PULMONARY EMBOLISM',
    'NEUTROPENIA', 'THROMBOCYTOPENIA', 'CEREBROVASCULAR ACCIDENT',
    'CARDIAC FAILURE', 'ACUTE KIDNEY INJURY', 'ATRIAL FIBRILLATION',
    'DEATH', 'FALL', 'FEBRILE NEUTROPENIA', 'HAEMATURIA', 'SEPSIS', 'DYSPNOEA',
    'ANAEMIA', 'HYPERTENSION', 'HEPATIC FAILURE',
]

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def get(url, tries=4):
    """GET with bounded retries. Raises RuntimeError if all attempts fail."""
    last = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'hp518-adr-fusion/1.0'})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read().decode('utf-8'))
        except Exception as e:                      # noqa: BLE001 - report and retry
            last = e
            sys.stderr.write('retry %d: %s\n' % (k, e))
            time.sleep(2 + 3 * k)
    raise RuntimeError('failed: %s (%s)' % (url, last))


def total(search):
    """Total number of matching reports (meta.results.total)."""
    j = get('%s?search=%s&limit=1' % (BASE, search))
    return j['meta']['results']['total']


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', default=os.path.join(ROOT, 'data', 'faers_extended.json'))
    args = ap.parse_args()

    res = {'date_window': DATE,
           'fetched_at': time.strftime('%Y-%m-%dT%H:%M:%S'),
           'drugs': {}}

    # --- OR-syntax sanity check -------------------------------------------
    t1 = total('patient.drug.medicinalproduct:%22ENZALUTAMIDE%22+AND+' + DATE)
    t2 = total('patient.drug.medicinalproduct:%22ABIRATERONE%22+AND+' + DATE)
    t3 = total('patient.drug.medicinalproduct:(%22ENZALUTAMIDE%22+%22ABIRATERONE%22)+AND+' + DATE)
    res['or_check'] = {'enza': t1, 'abi': t2, 'or_paren': t3, 'sum': t1 + t2}
    print('OR check:', res['or_check'])

    # --- per-drug top-200 PT profile + total reports ------------------------
    for key, name in DRUGS.items():
        s = 'patient.drug.medicinalproduct:%%22%s%%22+AND+%s' % (name, DATE)
        j = get('%s?search=%s&count=patient.reaction.reactionmeddrapt.exact&limit=200' % (BASE, s))
        prof = {d['term']: d['count'] for d in j['results']}
        tot = total(s)
        res['drugs'][key] = {'generic': name, 'total_reports': tot, 'pt_profile': prof}
        print('%-14s total=%7d pts=%d' % (key, tot, len(prof)))
        time.sleep(0.4)

    # --- whole-database total ---------------------------------------------
    res['db_total'] = total(DATE)
    print('db total =', res['db_total'])
    time.sleep(0.4)

    # --- pooled six-drug exposure + per-signal counts, DB-wide PT totals ----
    pool = '+'.join('%%22%s%%22' % n for n in DRUGS.values())
    sp = 'patient.drug.medicinalproduct:(%s)+AND+%s' % (pool, DATE)
    res['pooled_total_reports'] = total(sp)
    print('pooled 6-drug total =', res['pooled_total_reports'])
    time.sleep(0.4)

    res['pooled_signal_counts'] = {}
    res['db_pt_totals'] = {}
    for pt in SIGNALS:
        pte = pt.replace(' ', '+')
        n1 = total('patient.drug.medicinalproduct:(%s)+AND+patient.reaction.reactionmeddrapt:%%22%s%%22+AND+%s'
                   % (pool, pte, DATE))
        res['pooled_signal_counts'][pt] = n1
        time.sleep(0.3)
        n2 = total('patient.reaction.reactionmeddrapt:%%22%s%%22+AND+%s' % (pte, DATE))
        res['db_pt_totals'][pt] = n2
        print('  %-32s pooled=%6d  DB=%8d' % (pt, n1, n2))
        time.sleep(0.3)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, 'w', encoding='utf-8') as f:
        json.dump(res, f, ensure_ascii=False, indent=1)
    print('saved ->', args.out)


if __name__ == '__main__':
    main()
