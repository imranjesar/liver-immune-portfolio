.PHONY: download convert sc bulk test
download:
	python3 scripts/download.py
convert:
	python3 scripts/convert_sc.py
sc:
	python3 scripts/validate_inputs.py
	Rscript scripts/run_r.R scripts/sc_analysis.R
bulk:
	bash scripts/run_bulk.sh
test:
	python3 -m unittest discover -s tests -v

figures:
	bash scripts/run_publication_figures.sh
