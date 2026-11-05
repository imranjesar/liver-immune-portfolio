suppressPackageStartupMessages({library(Seurat);library(ggplot2)})
s<-readRDS("results/sc/seurat.rds")
a<-read.delim("config/celltypes.tsv",stringsAsFactors=FALSE,colClasses="character")
stopifnot(!anyDuplicated(a$cluster),setequal(a$cluster,as.character(unique(Idents(s)))))
idx<-match(as.character(s$seurat_clusters),a$cluster)
for(n in setdiff(names(a),"cluster"))s[[n]]<-a[[n]][idx]
s$annotation_status<-"provisional_manual_marker_review"
s$hold_from_condition_analysis<-grepl("^hold_",s$qc_status)
s$passes_500_gene_sensitivity<-s$nFeature_RNA>=500
m<-s[[]];m$cell_id<-rownames(m)
con<-gzfile("results/sc/review/cell_annotations.csv.gz","wt")
write.csv(m,con,row.names=FALSE);close(con)
summary<-merge(a,as.data.frame(table(as.character(s$seurat_clusters))),by.x="cluster",by.y="Var1");names(summary)[names(summary)=="Freq"]<-"cells"
write.csv(summary,"results/sc/review/annotation_summary.csv",row.names=FALSE)
per_sample<-as.data.frame(table(s$sample_id,s$broad_type));names(per_sample)<-c("sample_id","broad_type","cells")
per_sample$within_sample_pct<-100*per_sample$cells/as.numeric(table(s$sample_id)[per_sample$sample_id])
write.csv(per_sample,"results/sc/review/annotated_sample_composition.csv",row.names=FALSE)
p<-DimPlot(s,group.by="broad_type",shuffle=TRUE,seed=20261103,raster=TRUE)+ggtitle("Provisional lineage annotation; unresolved cells retained")
ggsave("results/sc/review/annotated_umap.png",p,width=12,height=7,dpi=160,bg="white")
p<-DimPlot(s,group.by="qc_status",shuffle=TRUE,seed=20261103,raster=TRUE)+ggtitle("Review status: holds are not confirmed doublet calls")
ggsave("results/sc/review/qc_status_umap.png",p,width=11,height=6,dpi=150,bg="white")
p<-ggplot(per_sample,aes(sample_id,within_sample_pct,fill=broad_type))+geom_col()+theme_bw()+labs(x="Mouse sample tag",y="Baseline-QC cells (%)",fill="Provisional lineage")
ggsave("results/sc/review/annotated_composition.png",p,width=12,height=7,dpi=150,bg="white")
# Reapply annotation to the baseline object with this script; no duplicated large RDS is required.
writeLines(capture.output(sessionInfo()),"results/sc/review/annotation_sessionInfo.txt")
