suppressPackageStartupMessages({library(Seurat);library(Matrix);library(ggplot2)})
s<-readRDS("results/sc/seurat.rds");m<-s[[]];m$cell_id<-rownames(m)
source("scripts/settings.R")
before<-read.csv("results/sc/qc_before.csv",row.names=1)
# Sensitivity only: no automatic rewriting of the object or its cluster assignments.
thresholds<-c(200,500,750)
r<-do.call(rbind,lapply(thresholds,function(th){d<-before;d$keep<-d$nFeature_RNA>=th&d$nFeature_RNA<=max_features&d$percent.mt<=max_mt
 x<-aggregate(keep~sample_id,d,function(z)c(retained=sum(z),total=length(z)));data.frame(min_features=th,sample_id=x$sample_id,retained=x$keep[,1],total=x$keep[,2])}))
r$retained_pct<-100*r$retained/r$total
write.csv(r,"results/sc/review/qc_threshold_sensitivity.csv",row.names=FALSE)
r2<-do.call(rbind,lapply(thresholds,function(th){x<-aggregate(nFeature_RNA~seurat_clusters,m,function(z)c(retained=sum(z>=th),total=length(z)));data.frame(min_features=th,cluster=x$seurat_clusters,retained=x$nFeature_RNA[,1],total=x$nFeature_RNA[,2])}))
write.csv(r2,"results/sc/review/cluster_threshold_sensitivity.csv",row.names=FALSE)
g<-c("Cd79a","Cd79b","Ms4a1","Cd3d","Cd3e","Trac","Csf1r","Adgre1","Fcgr1","Ly6g","Retnlg","S100a8","Ltf","Ngp","Camp")
x<-GetAssayData(s,layer="counts")[g,,drop=FALSE]
m$B_umi<-Matrix::colSums(x[c("Cd79a","Cd79b","Ms4a1"),]);m$T_umi<-Matrix::colSums(x[c("Cd3d","Cd3e","Trac"),]);m$granule_umi<-Matrix::colSums(x[c("S100a8","Ltf","Ngp","Camp","Retnlg","Ly6g"),])
m$B_two_markers_2umi<-Matrix::colSums(x[c("Cd79a","Cd79b","Ms4a1"),]>=2)>=2
m$T_two_markers_2umi<-Matrix::colSums(x[c("Cd3d","Cd3e","Trac"),]>=2)>=2
m$B_T_strong_screen<-m$B_two_markers_2umi&m$T_two_markers_2umi
out<-do.call(rbind,lapply(split(m,m$seurat_clusters),function(d)data.frame(cluster=as.character(d$seurat_clusters[1]),cells=nrow(d),median_B_umi=median(d$B_umi),median_T_umi=median(d$T_umi),median_granule_umi=median(d$granule_umi),strong_B_T_fraction=mean(d$B_T_strong_screen))))
write.csv(out,"results/sc/review/lineage_coexpression.csv",row.names=FALSE)
write.csv(m[,c("cell_id","sample_id","seurat_clusters","B_umi","T_umi","granule_umi","B_T_strong_screen")],"results/sc/review/cell_lineage_coexpression.csv",row.names=FALSE)
m$cluster<-m$seurat_clusters
p<-ggplot(m,aes(cluster,nFeature_RNA,fill=cluster))+geom_violin(scale="width",show.legend=FALSE)+geom_hline(yintercept=c(200,500,750),linetype="dashed")+scale_y_log10()+theme_bw()+labs(x="Baseline cluster",y="Detected genes (log scale)")
ggsave("results/sc/review/cluster_complexity.png",p,width=12,height=5,dpi=150,bg="white")
p<-ggplot(m,aes(B_umi+1,T_umi+1))+geom_point(size=.1,alpha=.15)+scale_x_log10()+scale_y_log10()+facet_wrap(~cluster,ncol=7)+theme_bw()+labs(x="B marker UMIs + 1",y="T marker UMIs + 1")
ggsave("results/sc/review/B_T_coexpression.png",p,width=14,height=9,dpi=150,bg="white")
