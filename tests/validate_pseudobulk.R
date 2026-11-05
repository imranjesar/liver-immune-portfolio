suppressPackageStartupMessages({library(Seurat);library(DESeq2);library(jsonlite)})
p<-"results/sc/pseudobulk";cfg<-jsonlite::fromJSON("config/pseudobulk.json")
stopifnot(file.exists(file.path(p,"AGGREGATION_COMPLETE")),file.exists(file.path(p,"DE_COMPLETE")))
x<-readRDS(file.path(p,"aggregates.rds"));sel<-read.csv(gzfile(file.path(p,"cell_selection.csv.gz")),stringsAsFactors=FALSE)
stopifnot(ncol(x$counts)==72,all(x$counts>=0),all(x$counts==round(x$counts)),sum(x$samples$n_cells)==sum(sel$included),!anyDuplicated(sel$cell_id))
stopifnot(all(!sel$needs_further_review[sel$included]),all(!sel$strong_B_T_screen[sel$included]))
for(j in seq_len(nrow(x$samples))){z<-x$samples[j,];stopifnot(z$n_cells==sum(sel$included & sel$population==z$population & sel$sample_id==z$sample_id,na.rm=TRUE),z$total_umi==sum(x$counts[,j]))}
elig<-read.csv(file.path(p,"population_eligibility.csv"))
stopifnot(sum(elig$tested)==6,all(!elig$tested[elig$population%in%c("Ltf_Ngp_neutrophils","Basophils")]))
r<-read.csv(file.path(p,"all_tested_populations_de.csv"))
stopifnot(!anyDuplicated(paste(r$population,r$gene)),all(r$padj_across_populations[!is.na(r$padj_across_populations)]>=0),all(r$padj_across_populations[!is.na(r$padj_across_populations)]<=1))
expected<-p.adjust(r$pvalue,"BH");stopifnot(isTRUE(all.equal(r$padj_across_populations,expected,tolerance=1e-12)))
for(pop in elig$population[elig$tested]){
 d<-readRDS(file.path(p,pop,"deseq2.rds"));m<-as.data.frame(colData(d));audit<-read.csv(file.path(p,pop,"design_audit.csv"))
 stopifnot(ncol(d)==6,all(table(m$condition)==3),!anyDuplicated(m$sample_id),all(m$n_cells>=cfg$min_cells_each_mouse),all(m$total_umi>=cfg$min_umi_each_mouse),audit$condition_capture_confounded)
 stopifnot("condition_aclf_ccs1477_vs_aclf_vehicle"%in%resultsNames(d))
 z<-r[r$population==pop,];stopifnot(all(z$loo_valid_fits<=6),all(z$loo_same_direction_fits<=z$loo_valid_fits))
}
# Independent spot checks against raw Seurat counts guard against aggregation of normalized values.
s<-readRDS("results/sc/strict500_consensus/seurat.rds")
for(pop in c("B_cells","CD8_T_cells","NK_cells"))for(tag in c("tag04","tag10")){
 gene<-c(B_cells="Cd79a",CD8_T_cells="Cd3e",NK_cells="Ncr1")[[pop]]
 ids<-sel$cell_id[sel$included & !is.na(sel$population) & sel$population==pop & sel$sample_id==tag]
 actual<-sum(GetAssayData(s,layer="counts")[gene,ids,drop=FALSE])
 stopifnot(actual==x$counts[gene,paste(pop,tag,sep="__")])
}
cat("Raw-count aggregation, selection, replicate units, contrast, FDR family and leave-one-mouse-out accounting verified\n")
