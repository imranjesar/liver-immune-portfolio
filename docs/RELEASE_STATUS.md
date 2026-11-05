# Release status

This repository is a local publication draft. It has not been pushed, assigned a remote, tagged or published. Multimodal AI work is on hold.

## Completed verification

- Public GEO input checksums, headers and sample-tag assignments verified.
- Python offline tests passed.
- Separate, initially empty R library hydrated from existing installed packages: all **219** versions exactly match `renv.lock`. A subsequent `renv::restore()` against that explicit library reported synchronization.
- Bulk analysis rerun into an independent output directory using that library. All **14 CSV result tables** were byte-identical, including annotations, QC, DE, nine omission fits, Hallmark ranks/results and descriptive cross-study context.
- Existing saved single-cell results passed doublet/cell-conservation, dependency-signature, annotation-fingerprint, pseudobulk raw-count/replicate/FDR and Hallmark/source-data checks under the isolated library. These are verification of existing outputs, not a new end-to-end single-cell computation.
- Six figure sets visually reviewed from rendered PDFs; all have one-page vector PDFs, embedded fonts, SVG and 600 dpi PNG. Export QA checks nominal dimensions within one PDF point (0.353 mm), allowing Cairo integer-point rounding, and records checksums.
- The public repository contains newly authored workflows and public data/results; unpublished personal research inputs are excluded. Application planning is outside the repository.

Machine-readable records: `metadata/reproduction_check.json`, `metadata/isolated-library-check.csv`, `metadata/publication_figure_manifest.json` and `metadata/environment-packages.tsv`.

## Remaining portability checks

Hydration reused installed binaries; it does **not** establish that every pinned source package can be downloaded and built on a fresh Ubuntu machine. No full source-download installation or complete single-cell rerun was performed in this check. Earlier execution interruptions and corrections remain in `metadata/VALIDATION.md`. Annotation fingerprints deliberately require review if recomputation changes cluster assignments.

GitHub Actions is configured for offline Python and figure-export checks, but has not run remotely. It does not currently perform R installation, GEO downloads or a complete analysis. Before describing this as a fully clean-install validated release, run the documented workflow on a fresh machine or equivalent CI runner and preserve logs.

## Publication preparation

The README, source tables, formal figure captions, study citations, gene-set license attribution and environment lock are included. Review the repository name, account, authorship and code/documentation license before creating a public repository. Third-party GEO inputs and MSigDB annotations retain their own terms; no license for them is implied by any future code license. No account identity or authorship has been invented.

A source archive supports review without publishing. Large raw matrices and serialized analysis objects are intentionally excluded and can be regenerated from public inputs.
