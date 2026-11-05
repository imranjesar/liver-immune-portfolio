# Hallmark analysis and publication figure methods

Mouse MSigDB MH v2026.2.Mm was downloaded unchanged from the Broad release archive. The collection has 50 orthology-mapped Hallmark sets; exact URL, version and SHA256 are in `metadata/hallmark_source.json`. The executed fgsea version is 1.34.2.

## Analysis

Each population is ranked by signed DESeq2 Wald statistic over finite, converged model results with a nonmissing P value. There is no significant-gene prefilter for enrichment. fgseaMultilevel uses overlaps of 15-500 genes, `eps=0`, 10,000 simple initialization permutations, and seed 20261103. BH is applied across all six populations separately within each rank variant. The sensitivity variant removes mitochondrial, operational ribosomal-protein, immunoglobulin-variable symbols and Gm42418 from the ranks; it is not a DE model refit or a batch correction.

## Results

| Population | Primary q <0.05 | Same direction and q <0.05 in both variants |
|---|---:|---:|
| Antibody_secreting_cells | 6 | 6 |
| B_cells | 9 | 9 |
| CD4_T_cells | 13 | 13 |
| CD8_T_cells | 21 | 21 |
| NK_cells | 6 | 5 |
| pDC | 9 | 9 |

Interferon-alpha/gamma and other immune signatures remain condition-associated under the rank-feature filter. Several growth/cell-cycle signatures have opposite enrichment directions. These are gene-set associations in this design, not validated changes in pathway activity or a mechanism of CCS1477. The allograft-rejection signature name does not imply a graft experiment.

The one primary association that loses the adjusted-P threshold is MYC targets V1 in NK cells: NES 1.682, global q=0.00229 in the primary ranking, versus NES 1.383, global q=0.0737 after the rank-feature filter. Its direction is unchanged but its support is sensitivity-dependent; Figure 3 marks it accordingly. Full primary and sensitivity tables, including all nonsignificant sets, are retained.

## Figure methods

Three main figures and one technical supplement are separated from explanatory diagnostics. Each has a nominal 180 mm canvas width, matched Liberation Sans typography, embedded-font vector PDF, editable SVG and 600 dpi PNG. GitHub previews use 240 dpi. Cairo rounds PDF page sizes to whole points, within 0.2 mm of the nominal dimensions. Legends and test definitions are supplied in [captions](../figures/CAPTIONS.md).

Figure source tables allow fast regeneration without loading the full Seurat object. UMAP and marker evidence originate from the retained-branch object. Expression estimates and confidence intervals come from the existing mouse-level DESeq2 models. The six illustrated genes are selected transparently from interferon leading edges in these same data; they are not a prespecified validation panel. Figure 3 shows every primary-associated gene set rather than a selected top-N subset. Its example curve is the smallest primary global-q pair.

PDF QA checks dimensions, font embedding, absence of raster images and PNG resolution. All PDF pages are rendered and visually inspected for clipping, overlap and legibility. Exact output checksums are recorded in `metadata/publication_figure_manifest.json`. Source-data and statistical contracts are checked in `tests/validate_hallmark_and_figure_data.R`.

## Reproduce

After restoring `renv.lock`, run `bash scripts/run_publication_figures.sh` to run enrichment, prepare figure source data, export figures and validate them. To redraw only, use `Rscript scripts/run_r.R scripts/publication_figures.R`. Poppler utilities are needed for export QA.

## Interpretation and attribution

Condition and capture library remain confounded. Shared mitochondrial/ribosomal features, unresolved annotation and within-lineage state composition merit further review; rank-feature stability does not resolve these limitations. Independent data are required before asserting a generalizable drug effect.

[Mouse MSigDB collections](https://www.gsea-msigdb.org/gsea/msigdb/mouse/collections.jsp) and [license terms](https://gsea-msigdb.org/gsea/msigdb_license_terms.jsp); Liberzon et al. (2015), [doi:10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004); [fgsea](https://bioconductor.org/packages/fgsea/).
