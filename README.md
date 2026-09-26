# HP518-ADR-Fusion

Multi-source evidence fusion pipeline for adverse drug reaction (ADR) prediction of the AR-PROTAC candidate HP518.

This repository accompanies the R2 revision of the manuscript under review at
*Frontiers in Pharmacology*. Every number reported in the revised manuscript and
its supplementary tables S1–S7 can be regenerated from the data committed here.

## Repository Structure

```
hp518-adr-fusion/
├── code/
│   ├── paths.py                     # shared, repo-relative path resolution
│   ├── faers_fusion_analysis.py     # Part 2: FAERS retrieval, deduplication & MedDRA mapping
│   ├── protac_risk_modeling.py      # Part 3: PROTAC risk model & severity classification
│   ├── generate_figures.py          # Part 4: Figure generation (bar, heatmap, radar, scatter)
│   ├── build_baseline_model.py      # Baseline ML models, original protocol (single split)
│   ├── adr_prediction_no_structure.py # Full pipeline entry point
│   ├── expanded_drug_analysis.py    # Extended analysis across the AR-PROTAC drug class
│   ├── analyze_sider_data.py        # SIDER database processing
│   ├── download_sider.py            # SIDER data download utility
│   ├── disproportionality_analysis.py # ROR/PRR, internal comparator (6 AR drugs)
│   ├── fetch_faers_extended.py      # R2: fetch the extended FAERS snapshot
│   └── r2_revision_analyses.py      # R2: reproduces supplementary tables S1-S7
├── data/
│   ├── ar_prostate_cancer_faers_data.json   # Processed FAERS data (AR prostate cancer drug reports)
│   ├── faers_extended.json                  # R2 FAERS snapshot (window 20040101-20260401)
│   ├── drug_names.tsv                       # Drug name mapping (SIDER)
│   └── meddra_all_se.tsv                    # Full MedDRA adverse events (SIDER, raw download)
└── results/
    ├── hp518_faers_fusion_report.json       # Final fusion scores for all candidate ADRs
    ├── hp518_expanded_drug_fusion_report.json
    ├── hp518_protac_risk_model.json         # PROTAC mechanism-specific risk model
    ├── baseline_model_results.json          # Baseline ML performance, original protocol
    ├── disproportionality_results.json      # ROR results, internal comparator design
    └── r2_supplementary_results.json        # R2: supplementary tables S1-S7
```

## Requirements

- Python 3.8+ (developed and verified on Python 3.9.13)
- `matplotlib`, `numpy`, `scikit-learn` (NumPy 1.26.4 / scikit-learn 1.6.1 were used for the R2 analyses)
- API: OpenFDA API (free, no authentication required)

```bash
pip install matplotlib numpy scikit-learn
```

## Configuration

All paths default to the repository root, so a fresh clone runs without setup:

```
data/     inputs
results/  generated reports
```

To point the pipeline at a different data/results tree, set:

```bash
export HP518_WORKSPACE=/path/to/other/workspace    # Linux/macOS
set HP518_WORKSPACE=C:\path\to\other\workspace     # Windows
```

## Pipeline Overview

### Part 1 — Data Retrieval
FAERS reports are retrieved via the **FDA OpenFDA API** (`faers_fusion_analysis.py`, Part 1).
Reports for enzalutamide, apalutamide, darolutamide, bicalutamide, flutamide and
abiraterone are collected and deduplicated.

The R2 snapshot in `data/faers_extended.json` covers
`receivedate:[20040101 TO 20260401]` and contains **20,344,365** case reports in
total. The date window is pinned in `fetch_faers_extended.py`; changing it
produces a different snapshot and changes every downstream number.

### Part 2 — Evidence Fusion & Scoring
- **FAERS signal mining**: disproportionality (ROR) at the MedDRA PT level
- **SIDER integration**: approved indications and adverse event frequencies from SIDER 4.1
- **Cross-database validation**: matching against clinical trial records
- **Severity classification**: PT-based severity tiers (critical / serious / non-serious)
- **Fusion scoring**: weighted combination of five components

| Component | Weight | Normalised | Cap |
|---|---|---|---|
| FAERS reporting frequency | 0.40 | 0.381 | 40 |
| SIDER label | 0.25 | 0.238 | 25 |
| Cross-drug sharing | 0.20 | 0.190 | 20 |
| Clinical concordance (V_match) | 0.15 | 0.143 | 15 |
| Second-generation AR bonus | 0.05 | 0.048 | 5 |
| **Sum** | **1.05** | 1.000 | **105** |

The weights sum to **1.05**, so raw fusion scores live on a **0–105 scale**; the
normalised column is provided only to make each stream's relative contribution
readable. The pair 0.38/0.22 printed in Table 1 of the submitted manuscript was
a transcription error and was never executed; re-scoring the identical candidate
pool with either pair leaves the top ten unchanged (Spearman ρ = 0.998), so no
reported result is affected.

**Two scoring scenarios** are reported (`r2_revision_analyses.py`, Table S1):

- **Scenario A — evidence integration.** The five weights above, with the
  clinical-concordance term included. Answers how strongly all available
  evidence, including what is already known about HP518 itself, converges on an ADR.
- **Scenario B — a priori prediction.** `V_match` is dropped (w = 0) and the
  freed weight is redistributed: FAERS 0.45, SIDER 0.25, cross-drug 0.25,
  second-generation 0.05 (sum = 1.00). Uses no HP518 outcome information and
  therefore gives the unbiased estimate of prospective performance.

Scenario A and Scenario B are on **different absolute scales** (1.05 vs 1.00).
They must be compared by rank only, never by absolute score. Rank correlation
between the two is ρ = 0.997.

### Part 3 — PROTAC-Specific Risk Modeling
Mechanism-specific ADR analysis for CRBN-based AR-PROTACs:
- E3 ligase connection analysis (CRBN/VHL)
- Warhead and linker toxicity signals
- AR degraders vs. AR inhibitors risk comparison

### Part 4 — Visualization
Figure generation: Top-30 ADR bar chart, Drug–ADR heatmap, Method comparison (AUPR), Risk radar, Priority scatter.

### Baseline Models (Supplementary Section S1.5)

Two protocols are reported. They answer different questions and both are retained.

**Original protocol** — `build_baseline_model.py`, single train/test split:
Logistic Regression AUROC 0.78 / AUPR 0.61; Naive Bayes AUROC 0.72 / AUPR 0.55.
The split was not stratified and no resampling was performed, so these figures
carry no variance estimate.

**Revised protocol** — `r2_revision_analyses.py` (Table S5), 20 repeats of
stratified 5-fold cross-validation with 2,000 bootstrap resamples for CIs:

| Model | AUROC | AUPR |
|---|---|---|
| Logistic Regression (L2) | 0.964 | 0.311 |
| Multinomial Naive Bayes | 0.940 | 0.236 |
| Fusion model (Scenario B) | 0.970 | 0.345 |

Features: min–max normalised pooled FAERS count, binary SIDER label, number of
reference drugs reporting the ADR, second-generation AR-inhibitor indicator.
Label: 1 if the ADR was observed in the Phase 1 HP518 dataset.

**No superiority claim is made.** The fusion model ranks above both baselines by
point estimate, but the bootstrap confidence intervals overlap almost entirely,
so these data cannot establish that the fusion framework outperforms machine
learning. The original protocol's larger AUPR values reflect a different (and
less rigorous) evaluation design, not a better model.

### Sensitivity and robustness (Tables S2, S3, S4, S6)
- **S2** ranking stability across per-drug cut-offs of 25–200 PTs
- **S3** retention of the serious-signal PTs under the 100-PT cut-off (19/19 retained)
- **S4** leave-one-reference-drug-out re-scoring
- **S6** 5,000 random (Dirichlet-like) weight draws plus one-at-a-time ±10%/±20% shifts

### Disproportionality (Table S7)
Two designs are reported and must not be confused:

- `disproportionality_analysis.py` → **internal comparator**: the other five
  AR-pathway drugs serve as background. Deliberately conservative, because
  comparator drugs share similar ADR profiles.
- `r2_revision_analyses.py` (S7) → **whole-database background**: the
  conventional design, using all 20,344,365 reports as the comparator.

## Usage

Run the complete pipeline:

```bash
python code/adr_prediction_no_structure.py   # Step 1: fusion analysis
python code/protac_risk_modeling.py          # Step 2: PROTAC risk model
python code/generate_figures.py              # Step 3: figures
python code/build_baseline_model.py          # Step 4: baseline ML (original protocol)
python code/disproportionality_analysis.py   # Step 5: ROR (internal comparator)
```

Reproduce the R2 supplementary tables S1–S7:

```bash
python code/r2_revision_analyses.py          # -> results/r2_supplementary_results.json
```

Re-fetch the FAERS snapshot (requires network; rewrites `data/faers_extended.json`):

```bash
python code/fetch_faers_extended.py
```

All R2 results are deterministic: every random draw is seeded (`SEED = 20260922`).

## Data Description

| File | Description | Size |
|------|-------------|------|
| `data/ar_prostate_cancer_faers_data.json` | Processed FAERS reports for AR drugs | 21 KB |
| `data/faers_extended.json` | R2 FAERS snapshot: 6 drugs × top-200 PTs, 18 serious-signal PTs | 36 KB |
| `data/drug_names.tsv` | Drug name mapping (generic ↔ brand) | 34 KB |
| `data/meddra_all_se.tsv` | Full MedDRA PT–ADR associations from SIDER | 19 MB |
| `results/hp518_faers_fusion_report.json` | Final ADR fusion scores (top candidates ranked) | 24 KB |
| `results/hp518_protac_risk_model.json` | PROTAC mechanism-specific risk model | 20 KB |
| `results/baseline_model_results.json` | Baseline ML performance (original protocol) | 3 KB |
| `results/disproportionality_results.json` | ROR results, internal comparator | 20 KB |
| `results/r2_supplementary_results.json` | R2 supplementary tables S1–S7 | — |

### SIDER counts quoted in the manuscript

SIDER 4.1 (`data/meddra_all_se.tsv`) contains **309,849** drug–side-effect pairs
covering **1,430** small-molecule drugs and **5,868** adverse-effect terms.

The 5,868 figure counts **distinct side-effect concepts**, i.e. distinct UMLS
concept identifiers in column 2. Counting the free-text names in column 6 gives
6,123 strings (4,251 if restricted to MedDRA PT rows). These are different
quantities; only the concept count should be quoted as the number of side effects.

## Notes

- All data in this repository are **fully anonymized** with no patient-identifiable information.
- FAERS data are retrieved via the **public FDA OpenFDA API** (aggregate counts only, no individual case data).
- SIDER data are publicly available at [SIDER 4.1](http://sideeffects.embl.de/media/download/mirror/).
- This repository is anonymized for double-blind peer review.

## Citation

Manuscript under review at Frontiers in Pharmacology.
