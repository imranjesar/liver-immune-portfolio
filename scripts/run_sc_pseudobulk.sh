#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
Rscript scripts/run_r.R scripts/sc_pseudobulk_prepare.R
Rscript scripts/run_r.R scripts/sc_pseudobulk_de.R
Rscript scripts/run_r.R scripts/sc_pseudobulk_figures.R
Rscript scripts/run_r.R scripts/sc_pseudobulk_audit.R
Rscript scripts/run_r.R tests/validate_pseudobulk.R
python3 scripts/report_sc_pseudobulk.py
