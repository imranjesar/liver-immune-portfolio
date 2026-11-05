#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
# Requires a successful renv::restore(); run_r.R selects the restored project library.
Rscript scripts/run_r.R scripts/sc_strict500.R
Rscript scripts/run_r.R scripts/sc_doublets.R
Rscript scripts/run_r.R scripts/sc_consensus_recluster.R
Rscript scripts/run_r.R scripts/sc_annotate_sensitivity.R
Rscript scripts/run_r.R scripts/sc_doublet_figures.R
Rscript scripts/run_r.R tests/validate_sc_sensitivity.R
Rscript scripts/run_r.R tests/validate_sc_sensitivity_annotations.R
python3 scripts/report_sc_sensitivity.py
