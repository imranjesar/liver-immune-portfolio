Latest analysis: [exploratory mouse-level pseudobulk](SC_PSEUDOBULK.md). Earlier baseline outputs are summarized below.

# Baseline results

## Bulk RNA-seq

Processed public counts: 48,714 gene rows, 12 samples. The primary analysis uses 9 liver samples (24h n=4; 72h n=5). After expression filtering, 7,327 genes remain; 7,325 have nonmissing BH adjusted P values. 1,002 meet FDR < 0.05. This count is descriptive, not a measure of model quality or causal evidence. Positive log2 fold changes indicate higher expression at 72h.

The PCA separates the time groups along PC1 but shows substantial variation among 72h samples. Inspect sample-level metadata and sensitivity to individual samples before interpreting pathways. No sample has been removed on the basis of its PCA position. Fold changes are unshrunk.

See `results/bulk/pca.pdf`, `MA.pdf`, `72h_vs_24h.csv` and `sessionInfo.txt`.

## Single-cell RNA-seq

Initial QC retained 55,634 cells across 9 sample tags. Provisional marker-based annotation and QC sensitivity review are available in `SC_QC_ANNOTATION.md`. Formal doublet modeling, sensitivity reclustering and mouse-level condition analysis remain pending. See `results/sc/`.
