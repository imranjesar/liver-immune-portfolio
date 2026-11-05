import csv,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[1]
rows=list(csv.DictReader((ROOT/'results/bulk/72h_vs_24h.csv').read_text().splitlines()))
padj=[float(r['padj']) for r in rows if r['padj'] not in ('NA','')]
s=f"""# Baseline results

## Bulk RNA-seq

Processed public counts: 48,714 gene rows, 12 samples. The primary analysis uses 9 liver samples (24h n=4; 72h n=5). After expression filtering, {len(rows):,} genes remain; {len(padj):,} have nonmissing BH adjusted P values. {sum(v<0.05 for v in padj):,} meet FDR < 0.05. This count is descriptive, not a measure of model quality or causal evidence. Positive log2 fold changes indicate higher expression at 72h.

The PCA separates the time groups along PC1 but shows substantial variation among 72h samples. Inspect sample-level metadata and sensitivity to individual samples before interpreting pathways. No sample has been removed on the basis of its PCA position. Fold changes are unshrunk.

See `results/bulk/pca.pdf`, `MA.pdf`, `72h_vs_24h.csv` and `sessionInfo.txt`.

## Single-cell RNA-seq

"""
q=ROOT/'results/sc/qc_after.csv'
if q.exists():
 cells=list(csv.DictReader(q.read_text().splitlines()))
 s+=f"Initial QC retained {len(cells):,} cells across {len(set(r['sample_id'] for r in cells))} sample tags. Provisional marker-based annotation and QC sensitivity review are available in `SC_QC_ANNOTATION.md`. Formal doublet modeling, sensitivity reclustering and mouse-level condition analysis remain pending. See `results/sc/`.\n"
else:s+='Conversion and full baseline analysis are in progress. No biological single-cell findings are claimed yet.\n'
if (ROOT/'docs/SC_PSEUDOBULK.md').exists(): s='Latest analysis: [exploratory mouse-level pseudobulk](SC_PSEUDOBULK.md). Earlier baseline outputs are summarized below.\n\n'+s
(ROOT/'docs/RESULTS.md').write_text(s)
