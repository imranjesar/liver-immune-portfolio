# Dataset validation — 3 November 2026

Primary evidence is the archived GEO family SOFT files in `metadata/`; each contains sample identities, design, processing and original supplementary URLs. They were downloaded directly from NCBI, not inferred from accession numbers. `config/files.tsv` lists the 13 individually downloadable files; `metadata/downloads-all.json` records observed byte counts and SHA256 checksums. Counts remain outside version control.

## Single-cell: GSE338157

[Series record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE338157) · [SOFT source](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE338nnn/GSE338157/soft/GSE338157_family.soft.gz)

Mouse liver CD45+ immune cells in fibrosis, ACLF plus vehicle, and ACLF plus CCS1477. Public since 14 July 2026. Three GEO samples (GSM9868022–GSM9868024) correspond to three condition-level multiplexed libraries. The design describes three individually sample-tagged biological replicate mice per group. Nine per-tag RSEC UMI CSV files and three tag-call CSV files are available (~195 MiB compressed altogether). The per-tag files supply the candidate mouse-level IDs; they are not nine independent library batches. Verify tag assignments against the tag-call files before inferential pseudobulk analysis.

Observed CSV orientation is **cells × genes**, with `Cell_Index` first, comment preamble, and 27,190 gene symbols in the first inspected matrix. The converter transposes into sparse gene × cell Matrix Market. Cell indexes are prefixed by tag before merging to avoid collisions across samples. Reference is GRCm38/mm10. GEO states BD Rhapsody pipeline v2.0; the first downloaded matrix header states v1.11 and January 2023 analysis. Preserve this discrepancy in reporting.

Initial output is unintegrated, exploratory QC, clustering and marker expression. Examine sample composition, doublets and ambient RNA before final annotation. Library and condition are confounded; neither Harmony nor any correction can recover an independently estimable library effect. Cell-level marker P values describe clusters and are not evidence of treatment efficacy. Condition comparisons require mouse-level pseudobulk after annotation and metadata review. No publication identifier was listed in the retrieved Series record; cite accession and contributors, check again before release.

## Bulk: GSE247373

[Series record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247373) · [SOFT source](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE247nnn/GSE247373/soft/GSE247373_family.soft.gz) · [Associated publication, PMID 38896080](https://pubmed.ncbi.nlm.nih.gov/38896080/)

Mouse neutrophils: spleen control n=3, liver 24h after CCl4 n=4, liver 72h n=5. Public since 15 June 2024. The space-delimited processed matrix contains 48,714 Ensembl gene rows and 12 named columns, ~382 KiB compressed. Column order differs from the GEO sample order; `config/bulk_samples.tsv` explicitly maps observed headers to GSM accessions. GRCm38 reference; counts generated from STAR alignments according to GEO.

Primary DESeq2 contrast is liver 72h versus liver 24h (positive log2 fold change means higher at 72h). Filter: at least 10 counts in at least 3 samples. Design: `~ condition`, Wald test, BH adjusted P values, alpha=0.05. No unreported covariates or pairing are invented. Inspect the paper for pooling and experimental-unit details before interpreting replicate independence. Spleen comparisons conflate tissue, time and treatment and are excluded from the primary contrast. DESeq2 output log2 fold changes are unshrunk. Cross-study comparisons are contextual, not a jointly batch-corrected experiment.

## Attribution and scope

This is an independent educational reanalysis of submitter-processed public data, not raw FASTQ reprocessing, original data generation or a replication claim. Cite both datasets, their investigators, the associated paper where available, GEO, and software. Public availability does not establish a blanket license for submitter data; retain original terms and attribution. No unpublished personal research materials are used.
