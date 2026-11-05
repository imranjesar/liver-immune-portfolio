"""Report exploratory pseudobulk comparisons without causal or batch-adjusted claims."""
import csv,gzip,pathlib,collections,json
ROOT=pathlib.Path(__file__).resolve().parents[1];OUT=ROOT/'results/sc/pseudobulk'
def read(name):
 p=OUT/name
 with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return list(csv.DictReader(f))
def main():
 if not (OUT/'DE_COMPLETE').exists():raise RuntimeError('Incomplete pseudobulk analysis')
 summary=read('de_summary.csv');elig=read('population_eligibility.csv');selection=read('selection_summary.csv');samples=read('sample_cell_counts.csv');res=read('all_tested_populations_de.csv')
 n=sum(int(x['cells'])for x in selection if x['reason']=='included');sig=[r for r in res if r['padj_across_populations']!='NA' and float(r['padj_across_populations'])<.05]
 lines=['# Exploratory mouse-level pseudobulk — CCS1477 versus vehicle','',
 '**Interpretation limit:** condition and capture library are completely confounded. There is one capture library per condition, with three mouse sample tags in each. Mouse-level aggregation avoids treating cells as independent replicates; it does not create replicate libraries or identify a library-adjusted drug effect. All P values and adjusted P values below are conditional on the stated model, not confirmatory evidence of efficacy.','',
 'Dataset: independent reanalysis of public mouse liver CD45+ immune-cell counts from [GSE338157](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE338157). No unpublished inputs are used.','',
 '## Cell-selection freeze','',
 f'{n:,} cells across all nine mouse tags were assigned to eight prespecified high-confidence populations. Each retained cell must have no baseline or new-cluster hold, no rate-sensitive model call and no strong same-cell B/T marker conflict. All cells had already passed the ≥500-gene branch and consensus doublet exclusion. The selection is conservative and does not prove all retained cells are singlets.','',
 'The population map is `config/pseudobulk_populations.tsv`; the frozen per-cell decision and one hierarchical exclusion reason are in `results/sc/pseudobulk/cell_selection.csv.gz`. Earlier reasons take priority when a cell fails multiple rules. Moderate/low-confidence populations are outside the primary analysis, even if plausible biological labels exist. The two excluded sparse populations refer to the reviewed Ltf/Ngp-rich neutrophil subset and basophils; no claim is made that all neutrophils are too sparse.','',
 '| Selection outcome | Cells |','|---|---:|']
 for r in selection:lines.append(f"| {r['reason']} | {int(r['cells']):,} |")
 lines+=['','## Replicates and eligibility','',
 'Raw integer UMI counts were summed by population and mouse sample tag. The 72 aggregates preserve eight populations × nine mouse tags. Fibrosis samples are retained in the aggregate/provenance tables but are not tested in this primary contrast. No zero-expression sample is invented to replace an absent population.','',
 'A population is tested only if **all six** ACLF mice have at least 20 selected cells and 10,000 summed UMIs. There are no dropped outlier mice, no pooling of mouse tags and no cell-level treatment test. The source GEO design describes the tags as biological mice; the analysis cannot independently verify that description.','',
 '| Population | Tested | Minimum cells among six mice | Minimum UMIs | Reason |','|---|---|---:|---:|---|']
 for r in elig:lines.append(f"| {r['population']} | {r['tested']} | {r['min_cells']} | {int(float(r['min_total_umi'])):,} | {r['reason']} |")
 lines+=['','## Model and multiple testing','',
 '- DESeq2 1.48.2, design `~ condition`, unpaired three CCS1477 and three vehicle mouse tags; positive log2 fold change means higher in CCS1477.',
 '- Including capture together with condition is rank-deficient. This is explicitly checked and exported for every tested population; no unidentifiable batch term is fitted.',
 '- Genes require ≥10 raw counts in at least three of the six mice. Size factors use `poscounts`; coefficient estimates are unshrunk. No normalized values are passed as count inputs.',
 '- Wald testing uses alpha=0.05, standard independent filtering and Cook outlier handling. Automatic count replacement is disabled. Failed coefficient convergence, if present, is not treated as a successful test.',
 '- `padj` is the DESeq2 within-population adjusted value. `padj_across_populations` applies BH jointly to the nonmissing P values from all six tested populations in this single contrast. Gene–population pairs are distinct tests, not necessarily distinct genes.',
 '- Six leave-one-mouse-out fits use the same prefiltered gene set with two versus three mice. The exported direction counts and LFC ranges measure sensitivity to omission, not repeat significance or correction for shared-library effects. No LOO P value is used as an extra discovery test.','',
 '| Population | Prefiltered genes | Nonmissing P values | Within-population BH <0.05 | Across-population BH <0.05 | Same LFC direction in all six omission fits, among across-population discoveries |','|---|---:|---:|---:|---:|---:|']
 for r in summary:lines.append(f"| {r['population']} | {r['genes_prefiltered']} | {r['nonmissing_p']} | {r['within_population_FDR05']} | {r['across_population_FDR05']} | {r['significant_stable_direction_all6']} |")
 lines+=['',f'The model yields {len(sig):,} gene–population pairs at across-population BH <0.05, representing {len(set(r["gene"] for r in sig)):,} unique gene symbols. These totals are descriptive model outputs and are not a score for portfolio quality or evidence of causal treatment activity.','',
 '## Sample coverage of tested populations','',
 '| Population | VEH tag04 | VEH tag05 | VEH tag06 | CCS tag10 | CCS tag11 | CCS tag12 |','|---|---:|---:|---:|---:|---:|---:|']
 for pop in [r['population']for r in elig if r['tested']=='TRUE']:
  lookup={r['sample_id']:r['n_cells']for r in samples if r['population']==pop}
  lines.append('| '+pop+' | '+' | '.join(lookup[t]for t in ['tag04','tag05','tag06','tag10','tag11','tag12'])+' |')
 lines+=['','## Diagnostic interpretation','',
 'Volcano facets use separate axis ranges; fold changes remain unshrunk. Several top-ranked features recur across populations, including mt-Nd4l and Gm42418. Gm42418 is retained as an unresolved feature annotation; no function is assigned here. Ribosomal-protein features also appear among the NK/T-cell results. Antibody-secreting-cell rankings include immunoglobulin-variable genes, whose variation can reflect repertoire or state composition. A B-cell Ccl5 signal alone does not establish B-cell-specific chemokine regulation and merits ambient/mixed-cell review. None of these signals is automatically deleted from the statistical tables. Total mitochondrial UMI fractions do not rise uniformly across populations even though mt-Nd4l has a positive modeled coefficient in all six; one mitochondrial gene should not be treated as a complete mitochondrial program.','',
 'The feature-class audit reports raw UMI fractions for mitochondrial symbols (`^mt-`), an operational, non-exhaustive ribosomal-protein symbol set (`^Rp[sl][0-9]`), immunoglobulin-variable symbols (`^Ig[hkl]v`) and Gm42418 separately. Shared shifts can arise from technical differences or broad biological states. Because all populations in each condition share a capture library, cross-population recurrence and stable mouse-omission directions do not resolve that ambiguity.','',
 '![Feature-class audit](../results/sc/pseudobulk/feature_class_audit.png)','',
 '![Exploratory model plots](../results/sc/pseudobulk/exploratory_volcano_panels.png)','',
 '## Figures and tables','',
 'Each tested population directory contains sample-labelled PCA, VST sample distances, a heatmap of the top 20 model-ranked genes, raw-model results, normalized counts, size factors and six omission-fit LFC columns. Heatmap colors are row-scaled VST values, not fold changes; gene ranking is derived from these same data and is not an independent validation. Mice are labelled explicitly so within-group heterogeneity remains visible.','',
 'The complete joined result is `results/sc/pseudobulk/all_tested_populations_de.csv`. `gene_feature_flags.csv` labels the operational feature classes, and `shared_model_associations.csv` lists features associated in at least three populations. No pathway enrichment or target-mechanism claim is made at this stage.','',
 '## Reproduce and validation','',
 'After a successful `renv::restore()`, run `bash scripts/run_sc_pseudobulk.sh`. Its R launcher explicitly selects the restored project library. Current execution used the isolated matching library. Fresh-environment restoration remains a separate release check.','',
 'Aggregation checks conserve included cells and gene-wise raw UMI totals. Independent tests spot-check marker counts against the source Seurat count layer, verify the six biological sample units and contrast direction, reconcile selection and coverage, and check the across-population BH family and omission-fit accounting. Source/annotation/map hashes and package session records accompany the outputs.','',
 '## Next interpretation work','',
 'Review shared mitochondrial/ribosomal/Gm42418 signals and sample QC before interpreting pathway rankings. Inspect broad B/CD8 state composition and candidate gene expression in individual mice. If pathway analysis is pursued, use the actual tested-gene universe, retain the confounding caveat, and distinguish exploratory hypotheses from mechanisms. An independent study would be needed to test whether these patterns generalize.','',
 'Method reference: [DESeq2 documentation](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html), including raw-count input, model-rank limitations and sample-level analysis. The executed package is the locked 1.48.2 version; the release website may describe a newer version.','']
 (ROOT/'docs/SC_PSEUDOBULK.md').write_text('\n'.join(lines))
 print('Pseudobulk report written:',len(sig),'model-associated gene–population pairs')
if __name__=='__main__':main()
