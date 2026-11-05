"""Build the review report from generated evidence, with explicit inference limits."""
import csv,json,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[1]
OUT=ROOT/'results/sc/review'
def read(name):
    with (OUT/name).open() as f:return list(csv.DictReader(f))
q=read('sample_qc.csv');s=read('qc_threshold_sensitivity.csv');a=read('annotation_summary.csv');co=read('lineage_coexpression.csv')
hold=sum(int(r['cells']) for r in a if r['qc_status'].startswith('hold_'))
retained=sum(int(r['cells_after']) for r in q)
lines=['# Single-cell QC and provisional annotation — 3 November 2026','',
'## Decision','',
'Keep the original 200–6,000 detected-gene / ≤20% mitochondrial baseline unchanged. The 500-gene threshold is the preferred **candidate for a separate QC rerun**, because it primarily removes low-complexity unresolved populations while retaining almost all cells in most marker-defined lineages. Do not overwrite the baseline, assume cluster IDs survive reclustering, or compare treatment groups before reviewing that rerun and within-tag doublets.','',
f'Baseline: {retained:,} cells across nine mouse tags. Five clusters (1, 4, 11, 12, 25), totaling {hold:,} cells ({100*hold/retained:.1f}%), are held from condition analysis because of low complexity, conflicting lineage, or unresolved cycling lineage. Other labels are provisional lineage assignments, not proof of singlet status.','',
'## Sample QC','',
'| Mouse tag | Before | Retained | Retained % | Median genes | Median UMIs | Median mitochondrial % |','|---|---:|---:|---:|---:|---:|---:|']
for r in q:lines.append(f"| {r['group']} | {int(r['cells_before']):,} | {int(r['cells_after']):,} | {float(r['retained_pct']):.1f} | {float(r['median_genes_after']):.0f} | {float(r['median_umi_after']):.0f} | {float(r['median_mt_after']):.2f} |")
lines+=['','All nine samples retain ≥96.6% of cells under the initial thresholds. High retention does not establish good quality: the initial cutoff admits large populations with very low RNA complexity. Mitochondrial medians alone do not identify those populations.','',
'## Threshold sensitivity','',
'| Minimum detected genes | Retained cells | Sample retention range |','|---|---:|---|']
for t in ['200','500','750']:
    rows=[r for r in s if r['min_features']==t]
    lines.append(f"| {t} | {sum(int(r['retained']) for r in rows):,} | {min(float(r['retained_pct']) for r in rows):.1f}–{max(float(r['retained_pct']) for r in rows):.1f}% |")
lines+=['',
'Cluster 1 has median 359 genes / 825 UMIs; cluster 4 has median 427.5 genes / 938.5 UMIs. Only 15.3% and 34.4%, respectively, pass 500 genes. They are concentrated in the ACLF vehicle and CCS1477 libraries, so a stricter threshold changes group representation. Remaining low-complexity cells should not automatically inherit an apparently normal lineage label. The 500-gene threshold also affects cycling cells and antibody-secreting cells and must be rechecked after reclustering.','',
'## Conflicting-lineage evidence','',
'The stronger B/T screen requires at least two of Cd79a/Cd79b/Ms4a1 and at least two of Cd3d/Cd3e/Trac each to have ≥2 raw UMIs in the **same barcode**. Cluster 11 passes at 86.7%, versus 0.15% in the main B cluster. Its median marker sums are 39 B-lineage and 22 T-lineage UMIs. This supports a mixed-lineage/doublet candidate annotation, not an activated CD4 identity. It does not prove doublets: physical cell conjugates or other mixed-cell events remain possible.','',
'Cluster 25 has 41 cells, strong B markers and a granulocyte transcript program (median 38 B marker and 33 granulocyte marker UMIs). Hold it as unresolved: ambient contamination and mixed barcodes cannot be distinguished here. The broad ≥1-UMI discordant-lineage screen is deliberately sensitive and can flag ambient RNA, especially in high-RNA myeloid populations; it is **not** a calibrated doublet classifier and does not justify automatic removal.','',
'## Provisional annotation','',
'Confidence describes lineage-marker agreement, not doublet-free status. Numeric clusters remain unchanged. Dot-plot colors are scaled within each gene and should not be interpreted as absolute expression comparisons between genes. Detailed support and conflicts are in `config/celltypes.tsv`; quantitative detection fractions and mean log-expression are in `results/sc/review/panel_evidence.csv`. The exported cell metadata can be joined to the baseline Seurat object by `cell_id`.','',
'| Cluster | Cells | Label | Confidence | Review status |','|---|---:|---|---|---|']
for r in sorted(a,key=lambda x:int(x['cluster'])):
    lines.append(f"| {r['cluster']} | {int(r['cells']):,} | {r['cell_type']} | {r['confidence']} | {r['qc_status']} |")
lines+=['',
'![Provisional lineage map](../results/sc/review/annotated_umap.png)','',
'![Cluster complexity](../results/sc/review/cluster_complexity.png)','',
'## Tag calls and doublet limits','',
'GEO tag-call tables include Multiplet and Undetermined calls. These are already absent from the nine supplied per-tag count matrices; no additional tag-derived barcode removal was performed. Counts below refer to the tag-call tables, not necessarily all physical capture events.','',
'| GEO library | Listed barcodes | Tag multiplets | Undetermined |','|---|---:|---:|---:|']
for r in json.loads((OUT/'tag_call_summary.json').read_text()):
    lines.append(f"| {r['file'].split('_')[0]} | {r['total_called_barcodes']:,} | {r['calls'].get('Multiplet',0):,} | {r['calls'].get('Undetermined',0):,} |")
lines+=['',
'Hashing cannot rule out within-tag doublets. Statistical doublet calling has **not** been run: scDblFinder is not installed in the locked environment. If added, use a compatible isolated environment, model the **three capture libraries**, and account for the already tag-filtered BD Rhapsody input. Do not treat nine mouse tags as nine independent capture batches, assume a 10x default doublet rate, or interpret tag multiplet fractions as the full true doublet rate. Expression profiles for discarded tag multiplets are not present in these per-tag matrices.','',
'## Remaining work before inference','',
'1. Run a separately named ≥500-gene sensitivity workflow, recluster and re-annotate. Check which conclusions survive; preserve baseline outputs.',
'2. Add capture-aware doublet modeling with BD-appropriate assumptions; inspect same-cell conflicts and suspicious high-complexity barcodes.',
'3. Assess ambient RNA limitations. Without an unfiltered droplet matrix, do not claim a validated soup correction.',
'4. Resolve cycling lymphocytes, NK/ILC-like and innate-like T subtypes using focused subclustering or an appropriate mouse reference. Treg-enriched is not synonymous with all-Treg.',
'5. Freeze reviewed cell types, then aggregate raw counts by mouse tag and cell type for pseudobulk. Library and condition remain confounded. No treatment-effect test is included in this review.','',
'## Reproduce','',
'Run `bash scripts/run_sc_review.sh` from the repository. The script reads the saved baseline Seurat object, writes review evidence and annotations, and tests metadata integrity. It preserves the baseline object and uses no unpublished inputs.','',
'## References and scope','',
'Manual labels were inferred from the observed public-data markers; they are not automated reference-mapping results. [ImmGen](https://www.immgen.org/) provides the mouse immune-lineage reference context. A [mouse lung immune atlas](https://elifesciences.org/articles/56890/figures) supplies independent tissue-context examples of DC and granulocyte marker patterns; cross-tissue agreement is supporting context, not a validation of this liver dataset. [Basophil differentiation study](https://www.nature.com/articles/s41467-023-38356-1) supports Mcpt8/Cd200r3 with little Kit as basophil evidence. [scDblFinder documentation](https://www.bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html) explains capture-aware handling of multiplexed samples. Consult documentation matching any version actually installed.','']
if (ROOT/'docs/SC_SENSITIVITY.md').exists(): lines[2:2]=['This is the first-pass review snapshot. See the [completed doublet sensitivity and reclustering stage](SC_SENSITIVITY.md) for the latest status.','']
(ROOT/'docs/SC_QC_ANNOTATION.md').write_text('\n'.join(lines))
print('QC and annotation review written; held cells:',hold)
