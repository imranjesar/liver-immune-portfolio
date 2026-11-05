"""Generate post-QC/doublet comparisons from completed, separate analysis branches."""
import csv,collections,json,pathlib,gzip
ROOT=pathlib.Path(__file__).resolve().parents[1]
SC=ROOT/'results/sc'
def read(path):
    opener=gzip.open if path.suffix==".gz" else open
    with opener(path,"rt") as f:return list(csv.DictReader(f))
def main():
    for stage in ['strict500','doublets','strict500_consensus']:
        if not (SC/stage/'COMPLETE').exists():raise RuntimeError(f'Incomplete analysis: {stage}')
    d=read(SC/'doublets/cell_doublet_summary.csv');base=read(SC/'qc_after.csv');post=read(SC/'strict500_consensus/cell_metadata.csv');strict=read(SC/'strict500/cell_metadata.csv');rates=read(SC/'doublets/scenario_summary.csv')
    gone=[r for r in d if r['consensus_doublet']=='TRUE'];uncertain=[r for r in d if r['rate_sensitive_call']=='TRUE']
    byold=collections.defaultdict(list)
    for r in d:byold[r['baseline_cluster']].append(r)
    lines=['# Stricter QC, doublet sensitivity and reclustering','',
    '## Analysis choices','',
    'The baseline is preserved. A separate branch requires 500–6,000 detected genes and ≤20% mitochondrial counts, using the baseline QC gene universe. It is normalized and clustered again from raw counts. Statistical doublet detection then runs independently in the three GEO capture libraries, retaining the nine mouse tags as biological samples.','',
    'scDblFinder 1.22.0 with XGBoost 1.7.11.1 evaluates assumed residual doublet rates of 2.5%, 5% and 10%. These are stress scenarios after sample-tag filtering, **not measured BD Rhapsody loading rates**. Each uses the same seed, 1,352 selected features, 20 dimensions and 3,000 random artificial doublets per capture. Rate uncertainty uses the version-pinned default ±40% of the assumed rate. Scores are prediction scores, not calibrated probabilities. Calling thresholds optimize classification and need not yield the assumed percentage of calls.','',
    'The conservative exclusion branch removes only cells classified as doublets in **all three** scenarios. Cells called in a subset of scenarios remain flagged as rate-sensitive. A negative model call does not prove singlet status, especially for homotypic doublets, rare populations and ambient contamination. No original cluster is deleted wholesale.','',
    '## Cell accounting','',
    '| Stage | Cells |','|---|---:|',f'| Initial baseline QC | {len(base):,} |',f'| ≥500-gene branch before doublet exclusion | {len(d):,} |',f'| Consensus predicted doublets excluded | {len(gone):,} |',f'| Retained after consensus exclusion | {len(post):,} |',f'| Retained rate-sensitive calls | {len(uncertain):,} |','',
    f'The stricter branch has {len(set(r["seurat_clusters"] for r in strict))} clusters; the consensus-filtered branch has {len(set(r["seurat_clusters"] for r in post))}. New cluster IDs are independent of baseline IDs. All nine mouse tags remain.','',
    '## Capture and rate scenarios','',
    '| Capture library | Assumed residual rate | Cells tested | Predicted doublets | Predicted % |','|---|---:|---:|---:|---:|']
    for r in rates:lines.append(f"| {r['capture']} | {100*float(r['assumed_rate']):.1f}% | {int(r['cells']):,} | {int(r['predicted_doublets']):,} | {100*int(r['predicted_doublets'])/int(r['cells']):.2f}% |")
    lines+=['','## Previously concerning clusters','',
    '| Baseline cluster | Passed ≥500 genes | Consensus doublets | Consensus % |','|---|---:|---:|---:|']
    for cl in ['1','4','11','12','20','25']:
        rs=byold[cl];n=sum(r['consensus_doublet']=='TRUE'for r in rs)
        lines.append(f'| {cl} | {len(rs):,} | {n:,} | {100*n/len(rs):.1f}% |')
    lines+=['',
    'The mixed B/T cluster was flagged independently by same-cell lineage-marker co-expression before modeling; compare its model enrichment here without treating the manual screen as ground truth. Statistical doublet detection is not a low-complexity filter. Low-complexity unresolved barcodes and rare B/granulocyte-program barcodes can remain even when not called doublets. These review flags must survive into subsequent annotation decisions.','',
    '## Sample retention','',
    '| Mouse tag | Baseline QC | ≥500 genes | Consensus excluded | Retained |','|---|---:|---:|---:|---:|']
    for sample in sorted({r['sample_id'] for r in d}):
        lines.append(f"| {sample} | {sum(r['sample_id']==sample for r in base):,} | {sum(r['sample_id']==sample for r in d):,} | {sum(r['sample_id']==sample for r in gone):,} | {sum(r['sample_id']==sample for r in post):,} |")
    lines+=['',
    '## New clusters and marker review','',
    'Baseline labels shown in overlap plots are provenance only. They are not automatically transferred final annotations. Inspect new marker evidence and the new mapping document before freezing labels. Review remaining low-complexity cells, cycling lineage, NK/ILC-like cells and marker conflicts; do not interpret changing cluster numbers or UMAP geometry as treatment effects.','',
    '![New clusters and baseline lineage overlap](../results/sc/strict500_consensus/clusters_and_baseline_labels.png)','',
    '![Marker panel after consensus doublet exclusion](../results/sc/strict500_consensus/lineage_marker_panel.png)','',
    '## Reproduce and audit','',
    'Restore `renv.lock` in a project-local environment. Run `bash scripts/run_sc_sensitivity.sh` from the repository root after a successful restore. The wrapper uses `scripts/run_r.R` to select the restored project library; direct `Rscript scripts/analysis.R` commands alone do not guarantee that the restored packages are selected. Apply the reviewed mapping with `Rscript scripts/run_r.R scripts/sc_annotate_sensitivity.R` when the saved cluster-assignment signature matches. The doublet stage checkpoints each capture/rate result and only reuses it when its input/configuration/dependency signature matches. `COMPLETE` files distinguish completed stages. Large Seurat objects remain outside version control.','',
    'The matching dependency versions were installed in an isolated library; shared R packages were unchanged. Preliminary models using a mixed preinstalled Bioconductor stack were superseded. The final package inventory, model session information and source archive checksums document the execution environment. A full fresh-environment restore remains a release check.','',
    '## Remaining inference limits','',
    'The residual doublet rate is unknown, discarded tag-multiplet expression is unavailable, and no unfiltered droplet matrix is supplied for validated ambient-RNA correction. Condition and capture library remain confounded. Within-mouse biological heterogeneity and injury-associated cell damage can affect retention. These analyses support QC and annotation; they do not establish a treatment effect.','',
    'Method: [scDblFinder](https://bioconductor.org/packages/scDblFinder/); capture-aware multiplexing and expected-rate handling are documented in the version-pinned package manual and source.','']
    mapping=SC/'strict500_consensus/annotation_mapping.csv'
    annotation=SC/'strict500_consensus/cell_annotations.csv.gz'
    if mapping.exists() and annotation.exists():
        maps=read(mapping);ann=read(annotation)
        counts=collections.Counter(r['seurat_clusters'] for r in ann)
        flagged=sum(r['needs_further_review']=='TRUE' for r in ann)
        lines += ['## Reviewed provisional mapping','',f'All {len(ann):,} retained cells have provisional annotations. {flagged:,} retain at least one review flag: an earlier unresolved-cell flag, a newly unresolved-cluster flag, or a rate-sensitive model call. These are review flags, not additional automatic cell exclusions. Prior flags survive even when cells move into a clearer new cluster.','',
        'The new cycling cluster has Cd3e in 82%, Trac in 85% and Mki67 in 92% of cells, supporting a T-enriched cycling identity. Low-complexity cells remain unresolved (new cluster 5; median 638 genes). The rare B/granulocyte-program population persists as new cluster 23 (36 cells), so it remains held. CD4/Treg, innate-like T and NK/ILC boundaries are still tentative; no treatment-effect conclusion follows from these state labels.','',
        '| New cluster | Cells | Provisional label | Confidence | Review status |','|---|---:|---|---|---|']
        for r in sorted(maps,key=lambda x:int(x['cluster'])):
            lines.append(f"| {r['cluster']} | {counts[r['cluster']]:,} | {r['cell_type']} | {r['confidence']} | {r['qc_status']} |")
        lines += ['','The mapping is guarded by a SHA256 fingerprint of the cell-to-cluster assignments. If reclustering changes assignments, it must be reviewed again rather than silently reused.','',
        '![Provisional annotation after sensitivity analysis](../results/sc/strict500_consensus/annotated_umap.png)','']
    (ROOT/'docs/SC_SENSITIVITY.md').write_text('\n'.join(lines))
if __name__=='__main__':main()
