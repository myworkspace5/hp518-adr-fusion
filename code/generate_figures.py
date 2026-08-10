# -*- coding: utf-8 -*-
"""
AR-PROTAC ADR 预测 - Part 4: 可视化图表生成（修复版）
修复：字体问题、重叠问题、英文坐标轴
"""
import json
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

sys.stdout.reconfigure(encoding='utf-8')

workspace = r"C:\Users\Bazinga\.qclaw\workspace\arptpred"
output_dir = os.path.join(workspace, "results")
fig_dir = os.path.join(output_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

# 使用 DejaVu Sans（支持更多 Unicode 字符，避免方格）
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
plt.rcParams['axes.unicode_minus'] = False

# 加载数据
fusion_path = os.path.join(output_dir, "hp518_faers_fusion_report.json")
with open(fusion_path, 'r', encoding='utf-8') as f:
    fusion_data = json.load(f)

faers_raw_path = os.path.join(workspace, "data", "faers", "ar_prostate_cancer_faers_data.json")
with open(faers_raw_path, 'r', encoding='utf-8') as f:
    faers_raw = json.load(f)

protac_path = os.path.join(output_dir, "hp518_protac_risk_model.json")
protac_data = None
if os.path.exists(protac_path):
    with open(protac_path, 'r', encoding='utf-8') as f:
        protac_data = json.load(f)

print("=" * 60)
print("Part 4: Visualization (Fixed)")
print("=" * 60)

# ============================================================
# Figure 1: HP518 ADR Top 30 Bar Plot
# ============================================================
print("\n[Figure 1] Top 30 ADR Fusion Score Bar Plot...")

top30 = fusion_data["top80_predictions"][:30]
adr_names = [p["adr"][:35] for p in top30]
scores = [p["fusion_score"] for p in top30]
colors = ['#d62728' if p.get("hp518_validated") else '#1f77b4' for p in top30]

fig, ax = plt.subplots(figsize=(16, 10))
bars = ax.barh(range(len(adr_names)), scores, color=colors, edgecolor='white', linewidth=0.5)
ax.set_yticks(range(len(adr_names)))
ax.set_yticklabels(adr_names, fontsize=9)
ax.invert_yaxis()
ax.set_xlabel('Fusion Score (0-100)', fontsize=11)
ax.set_title('HP518 ADR Prediction - FAERS+SIDER Fusion Ranking (Top 30)\nRed = Clinically Validated in HP518 Phase 1', 
             fontsize=13, fontweight='bold')

for i, (bar, score) in enumerate(zip(bars, scores)):
    ax.text(score + 1, i, f'{score:.1f}', va='center', fontsize=7.5)

from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='#d62728', label='HP518 Clinical Validated'),
                   Patch(facecolor='#1f77b4', label='Predicted')]
ax.legend(handles=legend_elements, loc='lower right')
ax.set_xlim(0, max(scores) * 1.15)
plt.tight_layout()
fig.savefig(os.path.join(fig_dir, "fig1_top30_adr_barplot.png"), dpi=200, bbox_inches='tight')
plt.close()
print("  OK fig1_top30_adr_barplot.png")

# ============================================================
# Figure 2: Drug-ADR Heatmap (Fixed: English labels, no subscript)
# ============================================================
print("\n[Figure 2] Drug-ADR Heatmap...")

all_adrs = set()
drug_adr_matrix = {}
for drug_name, data in faers_raw["data"].items():
    drug_adr_matrix[drug_name] = data["adr_profile"]
    all_adrs.update(data["adr_profile"].keys())

adr_counts_total = Counter()
for profile in drug_adr_matrix.values():
    for adr, cnt in profile.items():
        adr_counts_total[adr] += cnt

top25_adrs = [adr for adr, _ in adr_counts_total.most_common(25)]
drug_names = list(drug_adr_matrix.keys())
# Use English drug names
en_names_map = {
    "enzalutamide": "Enzalutamide", "apalutamide": "Apalutamide", "darolutamide": "Darolutamide",
    "bicalutamide": "Bicalutamide", "flutamide": "Flutamide", "abiraterone": "Abiraterone"
}

matrix = np.zeros((len(top25_adrs), len(drug_names)))
for j, drug in enumerate(drug_names):
    for i, adr in enumerate(top25_adrs):
        matrix[i, j] = drug_adr_matrix[drug].get(adr, 0)

matrix_log = np.log10(matrix + 1)

fig, ax = plt.subplots(figsize=(14, 10))
im = ax.imshow(matrix_log, cmap='YlOrRd', aspect='auto')

ax.set_xticks(range(len(drug_names)))
en_labels = [en_names_map.get(d, d) for d in drug_names]
ax.set_xticklabels(en_labels, fontsize=10, rotation=30, ha='right')
ax.set_yticks(range(len(top25_adrs)))
ax.set_yticklabels([a[:40] for a in top25_adrs], fontsize=8.5)

ax.set_xlabel('Prostate Cancer Drugs', fontsize=11)
ax.set_ylabel('Adverse Reaction (MedDRA PT)', fontsize=11)
ax.set_title('FAERS ADR Frequency Heatmap\nTop 25 ADRs x AR/CYP17 Inhibitors (log10 scale)', 
             fontsize=12, fontweight='bold')

cbar = plt.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label('log10(Report Count + 1)', fontsize=9)

plt.tight_layout()
fig.savefig(os.path.join(fig_dir, "fig2_drug_adr_heatmap.png"), dpi=200, bbox_inches='tight')
plt.close()
print("  OK fig2_drug_adr_heatmap.png")

# ============================================================
# Figure 3: Multi-source Evidence Stacked Bar Chart
# ============================================================
print("\n[Figure 3] Multi-source Evidence Stacked Bar Chart...")

top15 = fusion_data["top80_predictions"][:15]
adr_names_15 = [p["adr"][:28].replace(" ", "\n") for p in top15]

faers_vals = []
sider_vals = []
sharing_vals = []
hp518_vals = []

for p in top15:
    faers_vals.append(p.get("faers_total_count", 0))
    sider_vals.append(25 if p.get("sider_present") else 0)
    sharing_vals.append(min(p.get("faers_drug_count", 0) / 7 * 20, 20))
    hp518_vals.append(15 if p.get("hp518_validated") else 0)

x = np.arange(len(adr_names_15))
width = 0.65

fig, ax = plt.subplots(figsize=(16, 9))
bars1 = ax.bar(x, faers_vals, width, label='FAERS Freq (norm)', color='#e74c3c', alpha=0.85)
bars2 = ax.bar(x, sider_vals, width, bottom=faers_vals, label='SIDER Label', color='#3498db', alpha=0.85)
bottom2 = [f + s for f, s in zip(faers_vals, sider_vals)]
bars3 = ax.bar(x, sharing_vals, width, bottom=bottom2, label='AR Sharing', color='#2ecc71', alpha=0.85)
bottom3 = [b + sh for b, sh in zip(bottom2, sharing_vals)]
bars4 = ax.bar(x, hp518_vals, width, bottom=bottom3, label='HP518 Clinical', color='#f39c12', alpha=0.95)

ax.set_ylabel('Evidence Score Component', fontsize=11)
ax.set_title('Multi-Source Evidence Decomposition for HP518 Top 15 Predicted ADRs\n(FAERS + SIDER + AR Drug Sharing + Clinical Validation)', 
             fontsize=11, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(adr_names_15, fontsize=7.5, rotation=45, ha='right')
ax.legend(loc='upper right', fontsize=9)

totals = [f + s + sh + h for f, s, sh, h in zip(faers_vals, sider_vals, sharing_vals, hp518_vals)]
for i, t in enumerate(totals):
    if t > 50:
        ax.text(i, t + 2, f'{t:.0f}', ha='center', va='bottom', fontsize=7.5, fontweight='bold')

plt.tight_layout()
fig.savefig(os.path.join(fig_dir, "fig3_multisource_evidence.png"), dpi=200, bbox_inches='tight')
plt.close()
print("  OK fig3_multisource_evidence.png")

# ============================================================
# Figure 4: Risk Priority Matrix (Fixed: no overlap, clean labels)
# ============================================================
print("\n[Figure 4] Risk Priority Matrix (Scatter Plot)...")

top80 = fusion_data["top80_predictions"]

def get_severity(adr):
    adr_l = adr.lower()
    critical_kw = ["death", "fatal", "cardiac arrest"]
    high_kw = ["heart failure", "cardiac failure", "myocardial", "pneumonitis",
               "interstitial lung", "hepatitis", "liver failure", "thrombo",
               "embolism", "neutropenia", "thrombocytopenia", "sepsis",
               "anaphylactic", "cerebrovascular"]
    moderate_kw = ["hypertension", "nausea", "vomiting", "diarrhoea", "fatigue",
                   "anaemia", "fall", "fracture", "depression", "rash",
                   "pruritus", "headache", "dizziness", "insomnia", "kidney",
                   "atrial"]
    
    for kw in critical_kw:
        if kw in adr_l: return 4
    for kw in high_kw:
        if kw in adr_l: return 3
    for kw in moderate_kw:
        if kw in adr_l: return 2
    return 1

severities = [get_severity(p["adr"]) for p in top80]
scores = [p["fusion_score"] for p in top80]
adr_labels = [p["adr"][:20] for p in top80]  # Shorter labels
validated = [p.get("hp518_validated", False) for p in top80]

# Deduplicate validated entries: keep only highest-scoring one per unique name
# (Top80 has 10 validated=True: 5 high-score + 5 low-score SIDER-only duplicates)
validated_indices = [i for i, v in enumerate(validated) if v]
validated_indices.sort(key=lambda i: scores[i], reverse=True)
seen_names = set()
deduped_indices = []
for i in validated_indices:
    norm = adr_labels[i].lower().strip()
    if norm not in seen_names:
        seen_names.add(norm)
        deduped_indices.append(i)
deduped_indices.sort(key=lambda i: scores[i])
validated_set = set(deduped_indices)

severity_colors = {1: '#95a5a6', 2: '#f1c40f', 3: '#e67e22', 4: '#c0392b'}
severity_labels = {1: 'Mild', 2: 'Moderate', 3: 'High', 4: 'Critical'}
marker_sizes = [220 if i in validated_set else 20 for i in range(len(scores))]

fig, ax = plt.subplots(figsize=(16, 10))

# Apply jitter to separate overlapping points with same score+severity
import random
random.seed(42)
jitter_x = [random.uniform(-0.8, 0.8) for _ in range(len(scores))]
jitter_y = [random.uniform(-0.04, 0.04) for _ in range(len(severities))]

# Plot all points by severity (with jitter)
for sev in [4, 3, 2, 1]:
    mask = [s == sev for s in severities]
    x_sc = [scores[i] + jitter_x[i] for i in range(len(scores)) if mask[i]]
    y_sc = [severities[i] + jitter_y[i] for i in range(len(severities)) if mask[i]]
    ms = [marker_sizes[i] for i in range(len(marker_sizes)) if mask[i]]
    
    sc = ax.scatter(x_sc, y_sc, c=severity_colors[sev], s=ms, alpha=0.75,
                    label=f'{severity_labels[sev]} (n={sum(mask)})', edgecolors='white', linewidth=0.5)

# Annotate only the 5 deduped validated points

# Place each label in a visually empty zone with curved leader line
# All offsets chosen to stay WITHIN plot bounds (x:0-115, y:0.5-4.5)
# Annotate only the 5 deduped validated points
# Use name-based label placement (reliable regardless of sort order)
label_zones = {
    'fatigue':    {'xytext': (-105, 48), 'ha': 'right',  'va': 'bottom', 'conn': 'arc3,rad=0.16'},
    'nausea':     {'xytext': (-8, 52),   'ha': 'center', 'va': 'bottom', 'conn': 'arc3,rad=0.08'},
    'decreased':  {'xytext': (30, 52),   'ha': 'left',   'va': 'bottom', 'conn': 'arc3,rad=0.10'},
    'constipation':{'xytext': (-82, -20), 'ha': 'right',  'va': 'top',    'conn': 'arc3,rad=0.10'},
    'vomiting':   {'xytext': (0, -32),   'ha': 'center', 'va': 'top',    'conn': 'arc3,rad=0.06'},
}
default_zone = {'xytext': (15, 15), 'ha': 'center', 'va': 'center', 'conn': 'arc3,rad=0'}

for idx, i in enumerate(deduped_indices):
    xi = scores[i] + jitter_x[i]
    yi = severities[i] + jitter_y[i]
    li = adr_labels[i]
    cfg = default_zone
    for kw, zone in label_zones.items():
        if kw in li.lower():
            cfg = zone
            break
    sev_color = severity_colors[severities[i]]
    
    ax.annotate(li, (xi, yi), textcoords="offset points", xytext=cfg['xytext'],
                fontsize=9, fontweight='bold', alpha=0.95,
                ha=cfg['ha'], va=cfg['va'],
                bbox=dict(boxstyle='round,pad=0.35', facecolor='white', edgecolor=sev_color, alpha=0.93, linewidth=1.1),
                arrowprops=dict(arrowstyle='-', color=sev_color, lw=1.2, alpha=0.45,
                                connectionstyle=cfg['conn']),
                zorder=10,
                clip_on=False)  # Allow labels to extend slightly beyond axes

ax.set_xlabel('Fusion Prediction Score', fontsize=12)
ax.set_ylabel('Severity Level', fontsize=12)
ax.set_yticks([1, 2, 3, 4])
ax.set_yticklabels(['Mild', 'Moderate', 'High', 'Critical'], fontsize=11)
ax.set_title('HP518 ADR Risk Priority Matrix\n(Large markers = Clinically Validated in Phase 1; Only validated points are labeled)', 
             fontsize=13, fontweight='bold')
ax.legend(loc='upper left', fontsize=10, framealpha=0.95)
ax.grid(True, alpha=0.3, axis='x')
ax.set_xlim(-5, 115)

plt.tight_layout()
fig.savefig(os.path.join(fig_dir, "fig4_risk_priority_matrix.png"), dpi=200, bbox_inches='tight')
plt.close()
print("  OK fig4_risk_priority_matrix.png")

# ============================================================
# Figure 5: Method Comparison Coverage Chart (Fixed: English labels)
# ============================================================
print("\n[Figure 5: Method Comparison Coverage Chart]...")

methods = ['SIDER Only\n(AR Similarity)', 'FAERS+SIDER\n(Fusion)']
coverage = [67, 83]
hp518_aes = ['NAUSEA', 'VOMITING', 'FATIGUE', 'CONSTIPATION', 'DIARRHOEA', 'DECREASED\nAPPETITE']

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Left: Coverage comparison
ax1 = axes[0]
methods_en = ['SIDER Only', 'FAERS+SIDER']
bars = ax1.bar(methods_en, coverage, color=['#3498db', '#e74c3c'], width=0.5, edgecolor='white', linewidth=2)
for bar, c in zip(bars, coverage):
    ax1.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 1, f'{c}%',
             ha='center', va='bottom', fontsize=14, fontweight='bold')
ax1.set_ylabel('Coverage Rate (%)', fontsize=11)
ax1.set_title('HP518 Clinical AE Coverage Comparison', fontsize=12, fontweight='bold')
ax1.set_ylim(0, 100)
ax1.axhline(y=100, color='gray', linestyle='--', alpha=0.5)

# Right: Per-AE hit/miss
ax2 = axes[1]
x_pos = np.arange(len(hp518_aes))
width = 0.35
sider_hits = [1, 1, 1, 0, 0, 1]  # NAUSEA, VOMITING, FATIGUE, CONSTIPATION, DIARRHOEA, DECREASED APPETITE
fusion_hits = [1, 1, 1, 0, 1, 1]  # Fusion captures DIARRHOEA (5/6=83%), still misses CONSTIPATION

bars1 = ax2.bar(x_pos - width/2, sider_hits, width, label='SIDER Only', color='#3498db', alpha=0.8)
bars2 = ax2.bar(x_pos + width/2, fusion_hits, width, label='FAERS+SIDER Fusion', color='#e74c3c', alpha=0.8)

ax2.set_xticks(x_pos)
ae_en = ['NAUSEA', 'VOMITING', 'FATIGUE', 'CONSTIPATION', 'DIARRHOEA', 'DECREASED\nAPPETITE']
ax2.set_xticklabels(ae_en, fontsize=9, rotation=15, ha='right')
ax2.set_ylabel('Predicted (1=Hit / 0=Miss)', fontsize=11)
ax2.set_title('Per-AE Hit/Miss by Method', fontsize=12, fontweight='bold')
ax2.set_ylim(0, 1.3)
ax2.legend(fontsize=9)
ax2.set_yticks([0, 1])
ax2.set_yticklabels(['Miss', 'Hit'], fontsize=10)

# Annotate with text instead of symbols
for bar, h in zip(bars1, sider_hits):
    ax2.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.05,
             'Hit' if h else 'Miss', ha='center', va='bottom', fontsize=9, fontweight='bold',
             color='green' if h else 'red')
for bar, h in zip(bars2, fusion_hits):
    ax2.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.05,
             'Hit' if h else 'Miss', ha='center', va='bottom', fontsize=9, fontweight='bold',
             color='green' if h else 'red')

plt.tight_layout()
fig.savefig(os.path.join(fig_dir, "fig5_method_comparison.png"), dpi=200, bbox_inches='tight')
plt.close()
print("  OK fig5_method_comparison.png")

# ============================================================
# Figure 6: ADR Category Radar Chart (if PROTAC data exists)
# ============================================================
if protac_data:
    print("\n[Figure 6: ADR Category Radar Chart]...")
    
    categories = ['GI', 'Fatigue/Asthenia', 'Cardiac Risk', 
                  'Pulmonary Risk', 'Hematologic', 'Hepatic/Renal',
                  'Neurologic', 'Skin/Rash']
    
    cat_keywords = {
        'GI': ['nausea', 'vomiting', 'diarrhoea', 'constipation', 'abdominal', 'appetite'],
        'Fatigue/Asthenia': ['fatigue', 'asthenia', 'malaise'],
        'Cardiac Risk': ['cardiac', 'atrial', 'myocardial', 'heart failure', 'hypertension'],
        'Pulmonary Risk': ['pulmonary', 'pneumonitis', 'interstitial lung', 'dyspnoea', 'embolism'],
        'Hematologic': ['anaemia', 'neutropenia', 'thrombocytopenia', 'hemorrhage'],
        'Hepatic/Renal': ['hepatic', 'liver', 'renal', 'kidney', 'aminotransferase', 'alkaline phosphatase', 'bilirubin'],
        'Neurologic': ['headache', 'dizziness', 'neuropathy', 'confusion', 'seizure', 'insomnia'],
        'Skin/Rash': ['rash', 'pruritus', 'dermatitis'],
    }
    
    scores_by_cat = {}
    for cat, kws in cat_keywords.items():
        cat_score = 0
        count = 0
        for p in fusion_data["top80_predictions"][:80]:
            adr_l = p["adr"].lower()
            if any(k in adr_l for k in kws):
                cat_score += p["fusion_score"]
                count += 1
        scores_by_cat[cat] = cat_score / max(count, 1)
    
    values = [scores_by_cat.get(cat, 20) for cat in categories]
    values += values[:1]  # Close the radar
    
    angles = np.linspace(0, 2 * np.pi, len(categories), endpoint=False).tolist()
    angles += angles[:1]
    
    fig, ax = plt.subplots(figsize=(9, 9), subplot_kw=dict(polar=True))
    ax.fill(angles, values, color='#e74c3c', alpha=0.3)
    ax.plot(angles, values, color='#c0392b', linewidth=2.5, marker='o', markersize=6)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, fontsize=11)
    # Expand rmax to give category labels enough room outside the axis
    ax.set_rmax(max(values) * 1.45)
    
    ax.set_title('HP518 ADR Category Profile\n(Radar View of Fusion Scores by Organ System)', 
                 fontsize=12, fontweight='bold', pad=20)
    
    for angle, value, cat in zip(angles[:-1], values[:-1], categories):
        # Place value label INSIDE (toward center) to avoid overlap with outer category labels
        inner_ratio = 0.85
        label_r = value * inner_ratio
        ax.annotate(f'{int(value)}', xy=(angle, label_r),
                   fontsize=10, ha='center', va='center', fontweight='bold', color='#c0392b',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='none', alpha=0.92))
    
    plt.tight_layout()
    fig.savefig(os.path.join(fig_dir, "fig6_adr_category_radar.png"), dpi=200, bbox_inches='tight', pad_inches=0.5)
    plt.close()
    print("  OK fig6_adr_category_radar.png")

# ============================================================
# Summary
# ============================================================
print()
print("=" * 60)
print(f"All figures saved to: {fig_dir}")
print("=" * 60)
for f in sorted(os.listdir(fig_dir)):
    fpath = os.path.join(fig_dir, f)
    size = os.path.getsize(fpath) / 1024
    print(f"  {f} ({size:.1f} KB)")
