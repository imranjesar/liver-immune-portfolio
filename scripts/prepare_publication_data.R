suppressPackageStartupMessages({library(Seurat);library(Matrix);library(jsonlite);library(DESeq2)})
out<-"figures/source_data";dir.create(out,recursive=TRUE,showWarnings=FALSE)
write_gz<-function(x,path){con<-gzfile(path,"wt");on.exit(close(con));write.csv(x,con,row.names=FALSE)}
s<-readRDS("results/sc/strict500_consensus/seurat.rds")
a<-read.csv(gzfile("results/sc/strict500_consensus/cell_annotations.csv.gz"),stringsAsFactors=FALSE)
a<-a[match(colnames(s),a$cell_id),];stopifnot(!anyNA(a$cell_id))
types<-c("B cells","T cells","NK/ILC-like cells","Monocytes","Macrophages","Neutrophils","pDC","cDC-like cells","Antibody-secreting cells","Eosinophil-like cells","Basophils","Low-complexity unresolved","Mixed-lineage unresolved")
displays<-c("B cells","T cells","NK/ILC-like","Monocytes","Macrophages","Neutrophils","pDC","cDC-like","Antibody-secreting","Eosinophil-like","Basophils","Low complexity","Mixed lineage")
info<-data.frame(broad_type=types,type_id=seq_along(types),display=displays)
write.csv(info,file.path(out,"atlas_lineage_key.csv"),row.names=FALSE)
e<-as.data.frame(Embeddings(s,"umap"));names(e)<-c("UMAP1","UMAP2")
e$cell_id<-rownames(e);e$broad_type<-a$broad_type;e$type_id<-match(e$broad_type,types)
stopifnot(!anyNA(e$type_id));write_gz(e,file.path(out,"atlas_coordinates.csv.gz"))
genes<-intersect(c("Cd79a","Cd3e","Ncr1","Csf1r","C1qa","Ly6g","Ltf","Siglech","Zbtb46","Jchain","Ccr3","Mcpt8"),rownames(s))
counts<-GetAssayData(s,layer="counts")[genes,,drop=FALSE];norm<-GetAssayData(s,layer="data")[genes,,drop=FALSE]
dots<-do.call(rbind,lapply(types,function(t){i<-which(a$broad_type==t);data.frame(broad_type=t,gene=genes,pct=100*Matrix::rowMeans(counts[,i,drop=FALSE]>0),mean_log=Matrix::rowMeans(norm[,i,drop=FALSE]),cells=length(i))}))
dots$mean_log_z<-ave(dots$mean_log,dots$gene,FUN=function(x)as.numeric(scale(x)))
write.csv(dots,file.path(out,"atlas_marker_evidence.csv"),row.names=FALSE)
rm(s,counts,norm);gc()
pop_order<-c("B_cells","CD4_T_cells","CD8_T_cells","NK_cells","pDC","Antibody_secreting_cells")
pc<-do.call(rbind,lapply(pop_order,function(p){x<-read.csv(file.path("results/sc/pseudobulk",p,"pca_coordinates.csv"));x$population<-p;x}))
write.csv(pc,file.path(out,"pseudobulk_pca.csv"),row.names=FALSE)
r<-read.csv("results/sc/pseudobulk/all_tested_populations_de.csv",stringsAsFactors=FALSE)
h<-read.csv("results/sc/pathways/hallmark_all_results.csv",stringsAsFactors=FALSE)
ifn<-h[h$variant=="all_valid_features"&h$pathway%in%c("HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE")&h$padj_across_populations<.05,]
le<-do.call(rbind,lapply(seq_len(nrow(ifn)),function(i)data.frame(gene=strsplit(ifn$leadingEdge[i],";",fixed=TRUE)[[1]],population=ifn$population[i])))
le_count<-table(unique(le)$gene)
sig<-r[!is.na(r$padj_across_populations)&r$padj_across_populations<.05,]
sig<-sig[!grepl("^mt-|^Rp[sl][0-9]|^Ig[hkl]v|^Gm[0-9]",sig$gene),]
elig<-names(le_count)[le_count>=3];sig<-sig[sig$gene%in%elig,]
group<-split(sig,sig$gene);tab<-do.call(rbind,lapply(group,function(z)data.frame(gene=z$gene[1],populations=length(unique(z$population)),median_global_q=median(z$padj_across_populations))))
tab<-tab[tab$populations>=3,];tab<-tab[order(tab$median_global_q,tab$gene),];chosen<-head(tab$gene,6);stopifnot(length(chosen)==6)
sel<-r[r$gene%in%chosen,];sel$illustration_order<-match(sel$gene,chosen)
write.csv(sel,file.path(out,"illustrative_ifn_gene_effects.csv"),row.names=FALSE)
write.csv(tab,file.path(out,"illustrative_gene_selection_candidates.csv"),row.names=FALSE)
jsonlite::write_json(list(rule="Top six by median global q among gene-population discoveries; gene significant in >=3 populations and IFN-alpha/gamma primary leading edge in >=3 populations; display exclusions: mt, Rp[sl][0-9], Ig variable, Gm-prefixed symbols",selected_genes=chosen,scope="Data-selected illustration, not independent validation or a prespecified gene panel"),file.path(out,"illustrative_gene_selection.json"),auto_unbox=TRUE,pretty=TRUE)
file.copy("results/sc/pathways/hallmark_all_results.csv",file.path(out,"hallmark_all_results.csv"),overwrite=TRUE)
file.copy("results/sc/pathways/rank_feature_sensitivity.csv",file.path(out,"hallmark_rank_sensitivity.csv"),overwrite=TRUE)
file.copy("results/sc/pseudobulk/feature_class_fractions.csv",file.path(out,"feature_class_fractions.csv"),overwrite=TRUE)
file.copy("results/sc/doublets/scenario_summary.csv",file.path(out,"doublet_scenarios.csv"),overwrite=TRUE)
file.copy("results/sc/pseudobulk/sample_cell_counts.csv",file.path(out,"pseudobulk_cell_coverage.csv"),overwrite=TRUE)
# Store the single curve chosen by the smallest primary global pathway q.
primary<-h[h$variant=="all_valid_features"&is.finite(h$padj_across_populations),];primary<-primary[order(primary$padj_across_populations,primary$population,primary$pathway),]
best<-primary[1,];write.csv(best,file.path(out,"illustrative_enrichment_choice.csv"),row.names=FALSE)
file.copy(file.path("results/sc/pathways",paste0(best$population,"_all_valid_features_ranks.csv")),file.path(out,"illustrative_enrichment_ranks.csv"),overwrite=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,"preparation_sessionInfo.txt"))

# Reconstruct variance fractions from the same model objects, without refitting DE.
ps<-c("B_cells","CD4_T_cells","CD8_T_cells","NK_cells","pDC","Antibody_secreting_cells")
pv<-do.call(rbind,lapply(ps,function(p){d<-readRDS(file.path("results/sc/pseudobulk",p,"deseq2.rds"));v<-varianceStabilizingTransformation(d,blind=TRUE);pc<-plotPCA(v,intgroup="condition",returnData=TRUE);data.frame(population=p,PC1_percent=100*attr(pc,"percentVar")[1],PC2_percent=100*attr(pc,"percentVar")[2])}))
write.csv(pv,file.path(out,"pca_variance.csv"),row.names=FALSE)
