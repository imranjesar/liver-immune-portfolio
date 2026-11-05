#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
Rscript scripts/run_r.R scripts/sc_review.R
Rscript scripts/run_r.R scripts/sc_qc_followup.R
python3 scripts/tag_call_qc.py
Rscript scripts/run_r.R scripts/sc_annotate.R
Rscript scripts/run_r.R tests/validate_sc_review.R
python3 scripts/report_sc_review.py
