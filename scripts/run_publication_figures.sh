#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
Rscript scripts/run_r.R scripts/sc_hallmark.R
Rscript scripts/run_r.R scripts/prepare_publication_data.R
Rscript scripts/run_r.R scripts/publication_figures.R
Rscript scripts/run_r.R tests/validate_hallmark_and_figure_data.R
python3 scripts/validate_publication_exports.py
