suppressPackageStartupMessages(library(jsonlite))
h<-read.csv("results/sc/pathways/hallmark_all_results.csv")
stopifnot(all(is.finite(h$NES)),all(h$pval>=0&h$pval<=1),!anyDuplicated(paste(h$population,h$pathway,h$variant)))
for(v in unique(h$variant)){
 z<-h[h$variant==v,];stopifnot(isTRUE(all.equal(z$padj_across_populations,p.adjust(z$pval,"BH"),tolerance=1e-12)))
}
coverage<-read.csv("results/sc/pathways/gene_set_coverage.csv")
stopifnot(length(unique(coverage$pathway))==50)
params<-jsonlite::fromJSON("metadata/hallmark_source.json")
stopifnot(identical(digest::digest(file="data/annotations/mh.all.v2026.2.Mm.symbols.gmt",algo="sha256"),params$sha256))
for(pop in unique(h$population)){
 a<-read.csv(paste0("results/sc/pathways/",pop,"_all_valid_features_ranks.csv"))
 b<-read.csv(paste0("results/sc/pathways/",pop,"_exclude_flagged_features_ranks.csv"))
 stopifnot(!anyDuplicated(a$gene),!anyDuplicated(b$gene),all(b$gene%in%a$gene),!any(grepl("^mt-|^Rp[sl][0-9]|^Ig[hkl]v",b$gene)),!any(b$gene=="Gm42418"))
 stopifnot(all(b$wald_stat==a$wald_stat[match(b$gene,a$gene)]))
}
c<-read.csv(gzfile("figures/source_data/atlas_coordinates.csv.gz"))
k<-read.csv("figures/source_data/atlas_lineage_key.csv")
d<-read.csv("figures/source_data/atlas_marker_evidence.csv")
stopifnot(nrow(c)==44188,!anyDuplicated(c$cell_id),all(is.finite(c$UMAP1)),all(is.finite(c$UMAP2)),all(c$type_id%in%k$type_id),all(d$pct>=0 & d$pct<=100),all(is.finite(d$mean_log_z)))
pca<-read.csv("figures/source_data/pseudobulk_pca.csv");stopifnot(nrow(pca)==36,all(table(pca$population,pca$condition)==3))
e<-read.csv("figures/source_data/illustrative_ifn_gene_effects.csv");stopifnot(nrow(e)==36,length(unique(e$gene))==6,all(is.finite(e$lfcSE)))
cat("Hallmark ranks, global adjustment, atlas counts, sample units and gene-effect source data verified\n")
