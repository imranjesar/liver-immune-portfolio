# Baseline validation — 3 November 2026

- Three offline Python tests passed: sparse orientation/count preservation, negative-count rejection and manifest/sample structure.
- All 13 downloaded public files match captured SHA256 values and sizes.
- Every cell in all nine matrices matches the expected tag-call assignment; no cell overlaps tags within its GEO library. Total input cells: 57,039.
- Bulk workflow completed with consistent Bioconductor 3.21/S4Vectors 0.46.0. The DE result CSV was byte-identical on repeat execution. SHA256: `74672e4caa48036e4e60f27805f1d06ac8c8185ba505e7db2c53b7709d749d84`.
- Single-cell analysis produced 55,634 post-QC cells, 26 clusters, cluster markers, composition tables and a readable Seurat object. During the initial execution an edit to the live script caused a trailing parsing error after output writing. The source was corrected; saved-object verification, figure export and result integrity checks subsequently completed successfully. A clean end-to-end execution of the final source revision remains pending.
- Result integrity checks passed for both workflows: intended bulk sample groups, valid adjusted P values, all nine single-cell sample tags retained, unique cell IDs, finite UMAP coordinates and explicit QC bounds.
- `results/sc/sessionInfo.txt` describes saved-object verification; the lockfile and package inventory record the analysis package stack. It is not presented as an uninterrupted original execution log.
- Figures were visually inspected. Cell-type annotation and mouse-level condition inference are not complete.
- Fresh-environment renv restore and end-to-end rerun are not yet verified. CI currently validates Python smoke tests only.

## QC and annotation review

The new review workflow exports marker detection fractions, sample/cluster QC, threshold sensitivity, tag-call counts, same-cell lineage conflicts and provisional labels for all 26 clusters. Review integrity checks cover every baseline barcode, complete mapping, composition denominators, sensitivity monotonicity and the observed B/T conflict in cluster 11. Five clusters (12,849 cells) are held from condition analysis. This is manual marker-based review, not statistical doublet calling or final treatment inference.

Review components completed individually, and the final annotation integrity checks passed. A repeat invocation of the combined review wrapper was terminated with exit 143 before completion; this is not recorded as an uninterrupted end-to-end success. Fresh-environment and complete-wrapper execution remain release checks.

## Completed sensitivity stage

Separate ≥500-gene clustering completed with 46,954 cells. Capture-aware scDblFinder models completed for three explicit rate scenarios using matched Bioconductor 3.21 dependencies; preliminary mixed-version outputs were superseded. Consensus doublet exclusion removed 2,766 candidates and left 44,188 cells across all nine mouse tags. The retained branch completed reclustering, marker generation and provisional re-annotation.

Tests passed for exact input membership, capture mapping, scenario uniqueness, majority-vote accounting, retained-cell conservation, all dependency signatures, finite embeddings in both saved Seurat objects, complete labels, preservation of review flags, cluster-assignment fingerprint and recorded model dependency versions versus the lockfile. The ordinary lockfile validation also passed. No treatment inference or fresh-environment restoration is claimed.

## Completed exploratory pseudobulk stage

A frozen conservative selection includes 30,171 cells across nine tags and eight prespecified high-confidence populations. Raw counts were summed into 72 population-by-mouse aggregates. Six populations met ≥20 cells and ≥10,000 UMIs in all six ACLF comparison mice; the reviewed Ltf/Ngp neutrophil subset and basophils were skipped as sparse. Six primary DESeq2 fits and 36 leave-one-mouse-out fits completed.

Tests passed for gene-wise UMI conservation, cell-selection/count reconciliation, independent marker-count spot checks against the raw Seurat layer, biological tag units, the CCS1477-versus-vehicle coefficient, rank deficiency of condition plus capture, BH calculation across tested populations, and omission-fit accounting. PCA, sample-distance, ranked-gene heatmaps and feature-class audits were generated and visually reviewed. No treatment efficacy, validated mechanism or fresh restore is claimed.

## Hallmark analysis and curated figure set

Mouse MSigDB MH v2026.2.Mm (50 sets) was analyzed with fgsea 1.34.2 using signed Wald ranks and a documented feature-exclusion sensitivity variant. Tests passed for valid rank membership, gene-set source SHA256, BH adjustment across populations within each variant, atlas cell counts, six-mouse PCA units and effect-table source values. The ordinary dependency snapshot passed with the added pathway and graphics packages.

Three main figures and one technical supplement were exported at nominal 180 mm width as vector PDF/SVG and 600 dpi PNG, with 240 dpi GitHub previews. All PDF fonts are embedded, and PDF inspection confirms no raster image objects, including color bars. Dimensions and PNG resolution passed automated checks; the figure manifest records final checksums. All four final PDFs were rendered and visually checked for clipping, overlap and readable labels. Initial checks exposed Cairo's whole-point page-size rounding and rasterized default color-bar legends; the tolerance now accounts for the recorded rounding and all color bars use vector rectangles.

The GitHub figure-export CI job is configured to check committed exports and manifest without rerunning models. It has not run on a remote repository in this session. A complete clean-environment restore remains unverified. Application scheduling was moved outside the prospective public repository. Explanatory diagnostics remain separate from the curated figure set.

## Curated bulk, publication figures and isolated-library check

- Bulk primary contrast: liver 72 h (five libraries) versus 24 h (four), no sample exclusion, fixed count >=10 in >=3 filter, 7,327 features. Explicit Ensembl-version preservation and all annotation mappings exported. Individual-mouse provenance/pooling remains unverified.
- 1,002 BH q<0.05 features (364 higher, 638 lower at 72 h); all retain direction with converged fits in all nine leave-one-library-out analyses. Primary and feature-filtered Hallmark rankings each yield no q<0.05 set (45 eligible sets; minimum primary q=0.148). Context comparison with single-cell is descriptive, with distinct contrasts and populations.
- The first pathway run failed because the worker backend attempted to open a local port. It was superseded by a successful complete run explicitly using `BiocParallel::SerialParam()`; no results from the failed run are presented as complete.
- Bulk tests passed raw count conservation, sample mapping, recomputed contrast/coefficient/adjusted P, annotation uniqueness, omission accounting and pathway BH checks.
- Two bulk figure sets exported as vector PDF/SVG and 600 dpi PNG. Both final PDFs rendered and visually inspected. All six figure sets pass font/vector/resolution/dimension checks. Dimension tolerance is now one PDF point (0.353 mm), reflecting integer-point Cairo page rounding, rather than the earlier insufficient 0.2 mm tolerance.
- Isolated library: 219 exact lockfile package versions hydrated from installed packages and synchronized by explicit-library restore. Initial default-path no-op excluded. An independent bulk rerun using that library produced 14 byte-identical CSV tables. Existing SC sensitivity, annotation, pseudobulk and Hallmark/source-data checks passed in the same environment. Full fresh source installation and end-to-end single-cell recomputation remain unverified.
- Public inputs and Python smoke tests rechecked successfully. GitHub Actions is configured but not executed remotely. No commit, remote, public repository or release was created. Multimodal work remains on hold.


## GitHub publication: November 5, 2026

- Published the prepared 365-file initial commit to https://github.com/imranjesar/liver-immune-portfolio on main; raw inputs, processed matrices, serialized R objects and authentication files were excluded. Publication audit found no private-key or GitHub-token patterns in selected text files.
- Initial hosted workflow passed both jobs: offline Python smoke tests and six-set figure-export verification. Run: https://github.com/imranjesar/liver-immune-portfolio/actions/runs/41953888049 . This does not certify a fresh R source installation or full analysis rerun.
- Confirmed the public README renders and all four main embedded figure previews load successfully. Existing clean-install/single-cell recomputation limitations remain documented. No code license or versioned release was chosen.
