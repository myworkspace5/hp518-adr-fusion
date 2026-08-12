# -*- coding: utf-8 -*-
"""
ROR & PRR Disproportionality Analysis (真实计算)
=================================================
数据局限说明：
  - FDA OpenAPI 仅返回 6 个 AR 药物的 FAERS 数据
  - "其他药物" = 另外 5 个 AR 药物（而非全 FAERS 背景）
  - 这是内部比较（internal comparator），结果更保守
  - S1.4 将如实披露此局限性

阈值标准（FDA 惯例）:
  - ROR: 95% CI lower bound > 1
  - PRR: PRR >= 2 且 chi^2 >= 4
  - 同时满足两项 → 显著 disproportionality signal
"""
import sys, os, json, math, zipfile, shutil
sys.stdout.reconfigure(encoding='utf-8')

# ============================================================
# 1. 加载数据
# ============================================================
FAERS_PATH = r"F:\2025-2026-2科研\HP518 ADR\Figures&Supplementary Files\hp518-adr-fusion\data\ar_prostate_cancer_faers_data.json"

with open(FAERS_PATH, 'r', encoding='utf-8') as f:
    raw = json.load(f)

drugs_data = raw['data']  # {drug_name: {total_reports, adr_profile: {adr: count}}}
drugs = list(drugs_data.keys())

# 9 个目标信号（Table 4）
TARGET_SIGNALS = [
    "Interstitial Lung Disease",
    "Myocardial Infarction",
    "Pulmonary Embolism",
    "Neutropenia",
    "Thrombocytopenia",
    "Cerebrovascular Accident",
    "Cardiac Failure",
    "Acute Kidney Injury",
    "Atrial Fibrillation",
]

# ADR 别名映射（处理命名差异）
ALIASES = {
    "Interstitial Lung Disease": ["Interstitial Lung Disease", "Interstitial lung disease", "LUNG INFILTRATE"],
    "Myocardial Infarction": ["Myocardial Infarction", "Myocardial infarction", "HEART ATTACK"],
    "Pulmonary Embolism": ["Pulmonary Embolism", "Pulmonary embolism", "Pulmonary embolism venous"],
    "Neutropenia": ["Neutropenia", "Neutropenia", "Neutrophil count decreased"],
    "Thrombocytopenia": ["Thrombocytopenia", "Thrombocytopenia", "Platelet count decreased"],
    "Cerebrovascular Accident": ["Cerebrovascular Accident", "Cerebrovascular accident", "CVA", "STROKE"],
    "Cardiac Failure": ["Cardiac Failure", "Cardiac failure", "Heart failure", "Cardiac failure congestive"],
    "Acute Kidney Injury": ["Acute Kidney Injury", "Acute kidney injury", "RENAL FAILURE", "Kidney injury"],
    "Atrial Fibrillation": ["Atrial Fibrillation", "Atrial fibrillation", "AF", "Auricular fibrillation"],
}

def find_adr_count(adr_profile, signal_name):
    """在 adr_profile 中查找信号计数（尝试多个别名）"""
    aliases = ALIASES.get(signal_name, [signal_name])
    for alias in aliases:
        # 精确匹配
        if alias in adr_profile:
            return adr_profile[alias]
        # 模糊匹配（小写）
        alias_lower = alias.lower()
        for k, v in adr_profile.items():
            if k.lower() == alias_lower:
                return v
    return 0

def get_all_adr_counts(drugs_data, signal_name):
    """获取所有药物中目标信号的计数"""
    counts = {}
    for drug, info in drugs_data.items():
        profile = info.get('adr_profile', {})
        counts[drug] = find_adr_count(profile, signal_name)
    return counts

# ============================================================
# 2. ROR 计算
# ============================================================
def compute_ror(target_drug, signal_name, drugs_data):
    """
    计算单个药物×信号的 ROR（内部比较池）
    2×2 表:
              Signal    Other ADR
    Target      a          b
    Other       c          d
    a = target_drug & signal
    b = target_drug & other_ADR = total_reports_target - a
    c = other_drugs & signal (sum over 5 other drugs)
    d = other_drugs & other_ADR = sum(total_reports_other) - c
    """
    target_info = drugs_data[target_drug]
    target_total = target_info['total_reports']
    target_profile = target_info.get('adr_profile', {})
    a = find_adr_count(target_profile, signal_name)

    # 其他药物
    other_total = 0
    c = 0
    for drug, info in drugs_data.items():
        if drug != target_drug:
            other_total += info['total_reports']
            other_profile = info.get('adr_profile', {})
            c += find_adr_count(other_profile, signal_name)

    b = target_total - a
    d = other_total - c

    if any(x <= 0 for x in [a, b, c, d]):
        return None  # 无法计算（计数为0）

    # ROR = (a/b) / (c/d) = ad/bc
    ror = (a * d) / (b * c)

    # SE(log ROR) = sqrt(1/a + 1/b + 1/c + 1/d)
    se_log_ror = math.sqrt(1/a + 1/b + 1/c + 1/d)

    # 95% CI for log ROR: log(ROR) ± 1.96 * SE
    log_ror = math.log(ror)
    ci_lower = math.exp(log_ror - 1.96 * se_log_ror)
    ci_upper = math.exp(log_ror + 1.96 * se_log_ror)

    # ROR 显著：95% CI lower bound > 1
    ror_sig = ci_lower > 1.0

    return {
        'drug': target_drug,
        'signal': signal_name,
        'a': a, 'b': b, 'c': c, 'd': d,
        'target_total': target_total,
        'other_total': other_total,
        'ror': ror,
        'ci_lower': ci_lower,
        'ci_upper': ci_upper,
        'ror_sig': ror_sig,
    }

# ============================================================
# 3. PRR 计算
# ============================================================
def compute_prr(target_drug, signal_name, drugs_data):
    """
    计算 PRR
    PRR = (a/(a+b)) / (c/(c+d))
    chi^2 = sum((O-E)^2/E) for 2x2 table
    显著：PRR >= 2 AND chi^2 >= 4
    """
    target_info = drugs_data[target_drug]
    target_total = target_info['total_reports']
    target_profile = target_info.get('adr_profile', {})
    a = find_adr_count(target_profile, signal_name)

    other_total = 0
    c = 0
    for drug, info in drugs_data.items():
        if drug != target_drug:
            other_total += info['total_reports']
            other_profile = info.get('adr_profile', {})
            c += find_adr_count(other_profile, signal_name)

    b = target_total - a
    d = other_total - c

    if any(x <= 0 for x in [a, b, c, d]):
        return None

    # PRR
    prr = (a / (a + b)) / (c / (c + d)) if (c + d) > 0 else None

    # Chi^2 (Yates correction)
    # chi^2 = ((|ad - bc| - 0.5)^2 * (a+b+c+d)) / ((a+b)(c+d)(a+c)(b+d))
    numerator = (abs(a * d - b * c) - 0.5 * (a + b + c + d)) ** 2
    denominator = (a + b) * (c + d) * (a + c) * (b + d)
    chi2 = numerator / denominator if denominator > 0 else 0

    prr_sig = (prr is not None and prr >= 2.0 and chi2 >= 4.0)

    return {
        'drug': target_drug,
        'signal': signal_name,
        'a': a, 'b': b, 'c': c, 'd': d,
        'prr': prr,
        'chi2': chi2,
        'prr_sig': prr_sig,
    }

# ============================================================
# 4. 主循环：计算所有 9 信号 × 6 药物
# ============================================================
print("=" * 70)
print("ROR & PRR 真实计算结果")
print("=" * 70)
print(f"\n数据来源: {FAERS_PATH}")
print(f"药物数: {len(drugs)} ({', '.join(drugs)})")
print(f"目标信号: {len(TARGET_SIGNALS)}")
print(f"注意: '其他药物' = 另外 5 个 AR 药物（内部比较池，非全 FAERS 背景）\n")

results = {}  # {signal: {drug: {ror_result, prr_result}}}

for signal in TARGET_SIGNALS:
    results[signal] = {}
    print(f"\n{'─'*65}")
    print(f"▶ {signal}")
    print(f"{'─'*65}")

    # 汇总：跨6药总计
    total_a = 0
    total_b = 0
    total_c = 0
    total_d = 0

    for drug in drugs:
        ror_r = compute_ror(drug, signal, drugs_data)
        prr_r = compute_prr(drug, signal, drugs_data)
        results[signal][drug] = {'ror': ror_r, 'prr': prr_r}

        if ror_r:
            total_a += ror_r['a']
            total_b += ror_r['b']
            total_c += ror_r['c']
            total_d += ror_r['d']

    # 池化总计（6药合计 vs 池外背景=0）
    # 对于汇总分析，计算跨药总体 ROR
    if total_a > 0 and total_b > 0 and total_c > 0 and total_d > 0:
        pool_ror = (total_a * total_d) / (total_b * total_c)
        se_log = math.sqrt(1/total_a + 1/total_b + 1/total_c + 1/total_d)
        pool_ci_lower = math.exp(math.log(pool_ror) - 1.96 * se_log)
        pool_ci_upper = math.exp(math.log(pool_ror) + 1.96 * se_log)
        pool_ror_sig = pool_ci_lower > 1.0
    else:
        pool_ror = pool_ci_lower = pool_ci_upper = None
        pool_ror_sig = False

    # 打印每药结果
    print(f"  {'药物':<15} {'a':>6} {'b':>7} {'ROR':>8} {'95%CI Lower':>11} {'ROR显著':>9} {'PRR':>7} {'χ²':>7} {'PRR显著':>8}")
    print(f"  {'─'*75}")

    ror_sig_count = 0
    prr_sig_count = 0

    for drug in drugs:
        ror_r = results[signal][drug]['ror']
        prr_r = results[signal][drug]['prr']

        if ror_r:
            ror_str = f"{ror_r['ror']:.3f}"
            ci_str = f"{ror_r['ci_lower']:.3f}"
            ror_s = "✓" if ror_r['ror_sig'] else "✗"
            if ror_r['ror_sig']:
                ror_sig_count += 1
        else:
            ror_str = ci_str = "N/A"
            ror_s = "N/A"

        if prr_r:
            prr_str = f"{prr_r['prr']:.3f}" if prr_r['prr'] else "N/A"
            chi2_str = f"{prr_r['chi2']:.2f}"
            prr_s = "✓" if prr_r['prr_sig'] else "✗"
            if prr_r['prr_sig']:
                prr_sig_count += 1
        else:
            prr_str = chi2_str = "N/A"
            prr_s = "N/A"

        a_val = ror_r['a'] if ror_r else 0
        b_val = ror_r['b'] if ror_r else 0
        print(f"  {drug:<15} {a_val:>6} {b_val:>7} {ror_str:>8} {ci_str:>11} {ror_s:>9} {prr_str:>7} {chi2_str:>7} {prr_s:>8}")

    print(f"\n  ★ 池化(6药总计): a={total_a}, b={total_b}, c={total_c}, d={total_d}")
    if pool_ror:
        print(f"    Pooled ROR={pool_ror:.3f} (95%CI: {pool_ci_lower:.3f}–{pool_ci_upper:.3f}) {'✓ 显著' if pool_ror_sig else '✗ 不显著'}")

    print(f"  → ROR显著药物数: {ror_sig_count}/6, PRR显著药物数: {prr_sig_count}/6")

# ============================================================
# 5. 汇总表：9个信号的总体判断
# ============================================================
print("\n" + "=" * 70)
print("汇总：9 个高优先级信号的 Disproportionality 分析")
print("=" * 70)
print(f"\n{'信号':<30} {'池化ROR':>10} {'95%CI下界':>10} {'ROR显著':>8} {'池化PRR':>10} {'χ²':>7} {'PRR显著':>8} {'双达标':>8}")
print(f"{'─'*90}")

summary = []
for signal in TARGET_SIGNALS:
    total_a = sum((results[signal][drug]['ror']['a'] if results[signal][drug]['ror'] else 0) for drug in drugs)
    total_b = sum((results[signal][drug]['ror']['b'] if results[signal][drug]['ror'] else 0) for drug in drugs)
    total_c = sum((results[signal][drug]['ror']['c'] if results[signal][drug]['ror'] else 0) for drug in drugs)
    total_d = sum((results[signal][drug]['ror']['d'] if results[signal][drug]['ror'] else 0) for drug in drugs)

    if all(x > 0 for x in [total_a, total_b, total_c, total_d]):
        pool_ror = (total_a * total_d) / (total_b * total_c)
        se_log = math.sqrt(1/total_a + 1/total_b + 1/total_c + 1/total_d)
        pool_ci = math.exp(math.log(pool_ror) - 1.96 * se_log)
        ror_sig = pool_ci > 1.0

        pool_prr = (total_a / (total_a + total_b)) / (total_c / (total_c + total_d))
        num_chi = (abs(total_a * total_d - total_b * total_c) - 0.5 * (total_a + total_b + total_c + total_d)) ** 2
        den_chi = (total_a + total_b) * (total_c + total_d) * (total_a + total_c) * (total_b + total_d)
        pool_chi2 = num_chi / den_chi if den_chi > 0 else 0
        prr_sig = pool_prr >= 2.0 and pool_chi2 >= 4.0
    else:
        pool_ror = pool_ci = pool_prr = pool_chi2 = None
        ror_sig = prr_sig = False

    both = "✓✓" if (ror_sig and prr_sig) else ("✓" if (ror_sig or prr_sig) else "✗")

    summary.append({
        'signal': signal,
        'pool_ror': pool_ror,
        'ci_lower': pool_ci,
        'ror_sig': ror_sig,
        'pool_prr': pool_prr,
        'chi2': pool_chi2,
        'prr_sig': prr_sig,
        'both': both,
    })

    ror_str = f"{pool_ror:.3f}" if pool_ror else "N/A"
    ci_str = f"{pool_ci:.3f}" if pool_ci else "N/A"
    prr_str = f"{pool_prr:.3f}" if pool_prr else "N/A"
    chi2_str = f"{pool_chi2:.2f}" if pool_chi2 is not None else "N/A"

    print(f"  {signal:<28} {ror_str:>10} {ci_str:>10} {'✓' if ror_sig else '✗':>8} {prr_str:>10} {chi2_str:>7} {'✓' if prr_sig else '✗':>8} {both:>8}")

# 达标统计
sig_both = sum(1 for s in summary if s['ror_sig'] and s['prr_sig'])
sig_ror_only = sum(1 for s in summary if s['ror_sig'] and not s['prr_sig'])
sig_prr_only = sum(1 for s in summary if s['prr_sig'] and not s['ror_sig'])
sig_neither = sum(1 for s in summary if not s['ror_sig'] and not s['prr_sig'])

print(f"\n★ 双达标(ROR+PRR): {sig_both}/9")
print(f"  仅 ROR 显著: {sig_ror_only}/9")
print(f"  仅 PRR 显著: {sig_prr_only}/9")
print(f"  均不显著: {sig_neither}/9")
print(f"\n  结论: {sig_both} 个信号同时满足 ROR 和 PRR 显著性标准")

# ============================================================
# 6. 保存结果 JSON
# ============================================================
OUT_PATH = r"F:\2025-2026-2科研\HP518 ADR\Figures&Supplementary Files\hp518-adr-fusion\results\disproportionality_results.json"
out = {
    'metadata': {
        'data_source': str(FAERS_PATH),
        'limitation': '内部比较池：其他药物=另外5个AR药物，非全FAERS背景，结果更保守',
        'threshold_ror': '95%CI lower bound > 1',
        'threshold_prr': 'PRR >= 2 AND chi^2 >= 4',
        'drugs': drugs,
        'target_signals': TARGET_SIGNALS,
    },
    'summary': [
        {
            'signal': s['signal'],
            'pooled_ror': round(s['pool_ror'], 4) if s['pool_ror'] else None,
            'ror_95ci_lower': round(s['ci_lower'], 4) if s['ci_lower'] else None,
            'ror_significant': s['ror_sig'],
            'pooled_prr': round(s['pool_prr'], 4) if s['pool_prr'] else None,
            'chi2': round(s['chi2'], 4) if s['chi2'] else None,
            'prr_significant': s['prr_sig'],
            'both_significant': s['ror_sig'] and s['prr_sig'],
        }
        for s in summary
    ],
    'detailed': [
        {
            'signal': signal,
            'drug': drug,
            'a': results[signal][drug]['ror']['a'] if results[signal][drug]['ror'] else 0,
            'b': results[signal][drug]['ror']['b'] if results[signal][drug]['ror'] else 0,
            'c': results[signal][drug]['ror']['c'] if results[signal][drug]['ror'] else 0,
            'd': results[signal][drug]['ror']['d'] if results[signal][drug]['ror'] else 0,
            'ror': round(results[signal][drug]['ror']['ror'], 4) if results[signal][drug]['ror'] else None,
            'ror_ci_lower': round(results[signal][drug]['ror']['ci_lower'], 4) if results[signal][drug]['ror'] else None,
            'ror_significant': results[signal][drug]['ror']['ror_sig'] if results[signal][drug]['ror'] else None,
            'prr': round(results[signal][drug]['prr']['prr'], 4) if results[signal][drug]['prr'] else None,
            'chi2': round(results[signal][drug]['prr']['chi2'], 4) if results[signal][drug]['prr'] else None,
            'prr_significant': results[signal][drug]['prr']['prr_sig'] if results[signal][drug]['prr'] else None,
        }
        for signal in TARGET_SIGNALS
        for drug in drugs
    ]
}

with open(OUT_PATH, 'w', encoding='utf-8') as f:
    json.dump(out, f, ensure_ascii=False, indent=2)

print(f"\n结果已保存: {OUT_PATH}")
