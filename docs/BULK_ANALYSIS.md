# Bulk neutrophil RNA-seq: a within-liver time comparison

## Dataset and experimental unit

Public [GSE247373](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247373) and [Maali et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11186811/) describe FACS-sorted mouse neutrophils after CCl4 injury. The paper reports male C57BL/6N mice, RNA integrity >7, stranded mRNA libraries and NovaSeq paired-end sequencing. This portfolio starts from the submitter's processed GRCm38 count table; it cannot assess read quality, alignment or sorting purity directly.

The explicit header-to-GSM mapping is `config/bulk_samples.tsv`. Four liver 24 h and five liver 72 h libraries enter the primary model. Three spleen 0 h libraries appear only in descriptive library QC. Tissue, injury exposure and time are inseparable in spleen-to-liver comparisons. No samples were discarded based on PCA. The 72 h libraries are heterogeneous.

**Replication uncertainty:** retrieved GEO records and the main-text methods do not establish individual-mouse provenance or pooling per RNA-seq library. The linked supplementary methods could not be retrieved in this session. We therefore report n=4 and n=5 **libraries**, rather than claiming independently verified mouse replication. Unrecorded batches cannot be estimated or adjusted. The time contrast is observational and cannot establish that a gene causes repair.

## Counts, identifiers and model

The source table has 48,714 features and 12 libraries. Versioned Ensembl IDs are retained throughout; version suffixes are removed only for lookup in `org.Mm.eg.db` 3.21.0. All identifier-to-symbol mappings are exported. Of 48,714 features, 30,138 receive one unambiguous symbol; 225 have multiple candidate symbols. Unmapped and ambiguous rows remain in gene-level DE. Distinct count rows are never combined by symbol.

DESeq2 1.48.2 fits `~ condition`, with 24 h as reference; positive log2 fold change means higher at 72 h. Features require counts >=10 in at least three of the nine liver libraries (7,327 tested features). Default median-ratio size factors, Wald tests, independent filtering and Cook's handling are retained; automatic count replacement is explicitly disabled. Fold changes are unshrunk. Missing adjusted values remain missing. Sample size factors and maximum Cook's distances are exported as diagnostics, not used for ad hoc sample removal.

A blinded variance-stabilizing transformation and the 500 most variable retained features define the PCA. Count-library totals are sums of supplied gene counts, not measured sequencing read depth. Mitochondrial fractions use unambiguously mapped `mt-` features and cannot recover unmapped mitochondrial genes.

## Gene results and sensitivity

At BH q<0.05, **1,002 features** are associated with time: 364 higher and 638 lower at 72 h. Nine additional fits omit each library in turn, keeping the primary feature universe fixed. All 1,002 retain the primary direction in every omission fit, with model convergence. This is a direction check, not proof of consistent effect magnitude, independent validation or a correction for unrecorded batch.

`results/bulk/curated/` contains full DE coefficients/SE/statistics/p-values, raw feature mapping, library QC, size factors, Cook's diagnostics, PCA, VST sample distances and all omission coefficients. Existing baseline results are preserved separately.

## Hallmark analysis

The same unchanged MSigDB mouse MH v2026.2.Mm GMT used by the single-cell project supplies 50 sets. Among converged genes with finite Wald statistic/p-value and an unambiguous symbol, duplicate symbols retain the feature with highest primary baseMean (ties: original ID). Selection does not use significance. Signed Wald statistics are ranked deterministically. fgseaMultilevel 1.34.2 uses seed 20261103, serial execution, eps=0, nPermSimple=10,000 and overlap sizes 15–500. BH correction is across the 45 eligible sets, separately for each rank variant.

The primary ranking and a sensitivity ranking excluding `^mt-`, `^Rp[sl][0-9]`, `^Ig[hkl]v` and Gm42418 both yield **zero Hallmark sets with BH q<0.05**. The minimum primary q is 0.148. The main figure shows the ten smallest primary q values for descriptive inspection and explicitly states that none passes the threshold. These results do not establish pathway activation or absence of biological programs; they reflect this contrast, feature coverage and gene-set definition. Gene-level DE and pathway enrichment ask different questions.

## Cross-study context

`cross_study_context.csv` joins pathway names and NES values to the six single-cell population results. This is descriptive only: bulk compares time in sorted neutrophils after acute CCl4 injury; single-cell compares CCS1477 with vehicle in ACLF and is capture-confounded. Neither matching direction nor shared gene-set names establishes replication. We do not merge expression matrices, batch-correct the two studies together or conduct a pooled treatment test. No robust primary bulk Hallmark association is available to claim pathway replication.

## Run and inspect

After restoring the locked R library and downloading public inputs:

```bash
bash scripts/run_bulk.sh
```

The context table additionally requires `results/sc/pathways/hallmark_all_results.csv` from the preceding single-cell workflow. Figure exports are vector PDF/SVG and 600 dpi PNG, with source tables and formal legends in `figures/CAPTIONS.md`. Tests verify raw-count conservation, library mapping, sample sizes, contrast direction, exact result recomputation, annotation policy, omission fits and pathway BH correction.
