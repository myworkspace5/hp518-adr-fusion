# -*- coding: utf-8 -*-
"""
Shared path resolution for the HP518-ADR-Fusion pipeline.

By default every path is resolved relative to the repository root, so a fresh
clone runs without any configuration. Set the environment variable
HP518_WORKSPACE to point the whole pipeline at a different data/results tree
(for example the original working tree used to produce the submission).

Layout accepted in this repository:
    data/     meddra_all_se.tsv, drug_names.tsv, *_faers_data.json
    results/  all generated JSON reports
"""
import os

WORKSPACE = os.environ.get("HP518_WORKSPACE") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))

DATA_DIR = os.path.join(WORKSPACE, "data")
RESULTS_DIR = os.path.join(WORKSPACE, "results")

# SIDER and FAERS files live in data/sider and data/faers in the original
# working tree but directly in data/ here. Accept either layout.
_SIDER_SUBDIR = os.path.join(DATA_DIR, "sider")
SIDER_DIR = _SIDER_SUBDIR if os.path.isdir(_SIDER_SUBDIR) else DATA_DIR

_FAERS_SUBDIR = os.path.join(DATA_DIR, "faers")
FAERS_DIR = _FAERS_SUBDIR if os.path.isdir(_FAERS_SUBDIR) else DATA_DIR

SIDER_MEDDRA_TSV = os.path.join(SIDER_DIR, "meddra_all_se.tsv")
SIDER_DRUG_NAMES_TSV = os.path.join(SIDER_DIR, "drug_names.tsv")
FAERS_PROSTATE_CANCER_JSON = os.path.join(FAERS_DIR, "ar_prostate_cancer_faers_data.json")


def ensure_dirs():
    """Create the output directories the pipeline writes into."""
    os.makedirs(RESULTS_DIR, exist_ok=True)
    os.makedirs(FAERS_DIR, exist_ok=True)
