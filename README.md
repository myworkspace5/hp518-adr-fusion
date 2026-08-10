# HP518-ADR-Fusion

Multi-source evidence fusion pipeline for adverse drug reaction (ADR) prediction of the AR-PROTAC candidate HP518.

## Repository Structure

```
hp518-adr-fusion/
├── code/
│   ├── faers_fusion_analysis.py       # Part 2: FAERS real-world data retrieval, deduplication & MedDRA mapping
│   ├── protac_risk_modeling.py        # Part 3: PROTAC risk model & severity classification
│   ├── generate_figures.py            # Part 4: Figure generation (bar, heatmap, radar, scatter, etc.)
│   ├── build_baseline_model.py        # Baseline ML models (Logistic Regression & Naive Bayes)
│   ├── adr_prediction_no_structure.py # Full pipeline entry point
│   ├── expanded_drug_analysis.py      # Extended analysis across AR-PROTAC drug class
│   ├── analyze_sider_data.py          # SIDER database processing
│   └── download_sider.py              # SIDER data download utility
├── data/
│   ├── ar_prostate_cancer_faers_data.json   # Processed FAERS data (AR prostate cancer drug reports)
│   ├── drug_names.tsv                 # Drug name mapping (SIDER)
│   └── meddra_all_se.tsv              # Full MedDRA adverse events (SIDER, raw download)
└── results/
    ├── hp518_faers_fusion_report.json       # Final fusion scores for all candidate ADRs
    ├── hp518_expanded_drug_fusion_report.json
    ├── hp518_protac_risk_model.json         # PROTAC mechanism-specific risk model
    └── baseline_model_results.json          # Baseline ML model performance metrics
```

## Requirements

- Python 3.8+
- Required packages: `matplotlib`, `numpy`, `scikit-learn` (for baseline models)
- API: OpenFDA API (free, no authentication required)

Install dependencies:
```bash
pip install matplotlib numpy scikit-learn
```

## Pipeline Overview

### Part 1 — Data Retrieval
FAERS reports are retrieved via the **FDA OpenFDA API** (`faers_fusion_analysis.py`, Part 1).
Reports for enzalutamide, apalutamide, darolutamide, and other AR drugs are collected and deduplicated.

### Part 2 — Evidence Fusion & Scoring
- **FAERS signal mining**: disproportionality analysis (ROR, PRR) at the MedDRA PT level
- **SIDER integration**: approved indications and adverse event frequencies from SIDER 4.1
- **Cross-database validation**: matching against clinical trial records
- **Severity classification**: PT-based severity tiers (critical / serious / non-serious)
- **Fusion scoring**: weighted combination of FAERS (38%), SIDER (22%), cross-validation (20%), VigiMatch (15%)
  - Bonus (+5%) applied when supported by ≥2 independent sources

### Part 3 — PROTAC-Specific Risk Modeling
Mechanism-specific ADR analysis for CRBN-based AR-PROTACs:
- E3 ligase connection analysis (CRBN/VHL)
- Warhead and linker toxicity signals
- AR degraders vs. AR inhibitors risk comparison

### Part 4 — Visualization
Figure generation: Top-30 ADR bar chart, Drug–ADR heatmap, Method comparison (AUPR), Risk radar, Priority scatter.

### Baseline Models (Supplementary Section S1.5)
- Logistic Regression with L2 regularization
- Multinomial Naive Bayes
- Features: FAERS normalized frequency per ADR
- Performance: Logistic AUROC=0.78 / AUPR=0.61; NB AUROC=0.72 / AUPR=0.55

## Usage

Run the complete pipeline:

```bash
# Step 1: Retrieve FAERS data and run fusion analysis
python code/adr_prediction_no_structure.py

# Step 2: Build PROTAC risk model
python code/protac_risk_modeling.py

# Step 3: Generate figures
python code/generate_figures.py

# Step 4: Run baseline ML models
python code/build_baseline_model.py
```

Or run individual parts:

```bash
python code/faers_fusion_analysis.py   # Parts 1+2: FAERS retrieval & fusion
python code/protac_risk_modeling.py     # Part 3: PROTAC risk model
python code/generate_figures.py          # Part 4: Visualization
```

## Data Description

| File | Description | Size |
|------|-------------|------|
| `data/ar_prostate_cancer_faers_data.json` | Processed FAERS reports for AR drugs | 21 KB |
| `data/drug_names.tsv` | Drug name mapping (generic ↔ brand) | 34 KB |
| `data/meddra_all_se.tsv` | Full MedDRA PT–ADR associations from SIDER | 19 MB |
| `results/hp518_faers_fusion_report.json` | Final ADR fusion scores (top candidates ranked) | 24 KB |
| `results/hp518_protac_risk_model.json` | PROTAC mechanism-specific risk model | 20 KB |
| `results/baseline_model_results.json` | Baseline ML performance metrics | 3 KB |

## Notes

- All data in this repository are **fully anonymized** with no patient-identifiable information.
- FAERS data are retrieved via the **public FDA OpenFDA API** (aggregate counts only, no individual case data).
- SIDER data are publicly available at [SIDER 4.1](http://sideeffects.embl.de/media/download/mirror/).
- This repository is anonymized for double-blind peer review.

## Citation

Manuscript under review at Frontiers in Pharmacology (ID: 1932530).
