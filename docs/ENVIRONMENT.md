# Environment and reproducibility status

Baseline execution uses Ubuntu, R 4.5.2 and Python 3.11.3. The locked package inventory is in `metadata/environment-packages.tsv`; each analysis saves `sessionInfo.txt`. `renv.lock` captures the packages used, with the caveat below. The repository does not auto-activate renv, so existing system packages remain usable during baseline development. Run `renv::restore()` explicitly in a fresh local clone.

## Environment correction and remaining release check

The preinstalled environment initially mixed Bioconductor 3.21 with S4Vectors 0.48.0 from 3.22. A separate local library was built with S4Vectors 0.46.0 from the Bioconductor 3.21 source mirror, without changing shared system packages. The bulk workflow was rerun with that library. The current lockfile was regenerated with ordinary validation enabled and records the consistent 3.21 package set, Seurat 5.4.0 and presto at its recorded GitHub commit. The single-cell workflow uses the preinstalled Seurat stack. Its saved-object verification session is recorded separately; see `metadata/VALIDATION.md` for execution details.

A fresh full `renv::restore()` and end-to-end rerun in a clean environment are still pending. Passing snapshot validation is not proof of a successful restore. Bioconductor's main host returned HTTP 403 during source download here; its Posit mirror supplied the S4Vectors source. If that restriction recurs during restore, configure an official Bioconductor mirror and preserve the locked versions. Do not update shared system packages as a shortcut.

## Checks implemented

- Offline tests check sparse orientation, exact counts, rejection of negative counts and manifest/sample structure.
- Input validation checks SHA256 and byte counts and tests each cell barcode against its library's sample-tag calls, rejecting overlap across mouse tags within a library.
- Analysis rejects missing samples and malformed/noninteger bulk counts.
- Result checks verify every expected sample survives QC, finite embedding coordinates, the intended bulk contrast, and valid adjusted P values.
- CI currently runs offline Python tests only; CI does not download GEO data or validate a fresh R environment.
- Exact reruns can compare result tables, while embeddings are assessed for integrity rather than byte-for-byte coordinates.

## Doublet-analysis environment

The model stack is now aligned to Bioconductor 3.21: scDblFinder 1.22.0, scater 1.36.0, scran 1.36.0, S4Vectors 0.46.0 and their matching dependencies, plus XGBoost 1.7.11.1. The ordinary lockfile validation passed without forcing it. Preinstalled 3.22 dependency versions were superseded in a separate local library; shared system libraries were unchanged. Exact model dependencies appear in `results/sc/doublets/run_parameters.json`; source archive URLs and SHA256 values are in `metadata/doublet-package-sources.json`.

A restore alone does not guarantee that a later plain Rscript command uses the project library. Use `Rscript scripts/run_r.R scripts/analysis.R` after restoring; the launcher calls `renv::load()` explicitly. The current analyses were executed with the isolated library selected explicitly. A fresh full restore and complete wrapper rerun remain unverified.

## Latest isolated-library and bulk reproduction check

All 219 locked package versions were verified in a new separate library populated by `renv::hydrate(packages=names(lock$Packages), library=..., sources=...)` from installed packages; `renv::restore(library=..., lockfile=...)` then reported synchronization. Explicit package/library arguments are essential: an empty directory with no analysis scripts and default hydration can discover no packages, and a restore using default paths may inspect the active library rather than the intended new library. That preliminary no-op was discarded and is not counted as a restore check.

Using the explicit isolated library, an independent bulk rerun reproduced all 14 CSV tables byte-for-byte, and existing single-cell output checks passed. This verifies locked versions and numerical reproducibility for that bulk rerun. It does not verify a new source download/build of all packages or an end-to-end rerun of single-cell models. See `docs/RELEASE_STATUS.md` and `metadata/reproduction_check.json` for the precise scope. CI now also checks publication export files; neither CI job runs the full R workflow.
