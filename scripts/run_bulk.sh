#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
Rscript scripts/run_r.R scripts/bulk_complete.R
Rscript scripts/run_r.R tests/validate_bulk.R
Rscript scripts/run_r.R scripts/bulk_figures.R
python3 scripts/validate_publication_exports.py
