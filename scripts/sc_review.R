suppressPackageStartupMessages({library(Seurat);library(Matrix);library(ggplot2)})
set.seed(20261103)
dir.create("results/sc/review",showWarnings=FALSE)
s <- readRDS("results/sc/seurat.rds")
m <- s[[]]; m$cell_id <- rownames(m)
before <- read.csv("results/sc/qc_before.csv",row.names=1)
summary_qc <- function(d,group) do.call(rbind,lapply(split(d,d[[group]]),function(x)data.frame(group=x[[group]][1],cells=nrow(x),median_genes=median(x$nFeature_RNA),median_umi=median(x$nCount_RNA),median_mt=median(x$percent.mt),p95_genes=unname(quantile(x$nFeature_RNA,.95)),p95_umi=unname(quantile(x$nCount_RNA,.95)),p95_mt=unname(quantile(x$percent.mt,.95)))))
a<-summary_qc(before,"sample_id");b<-summary_qc(m,"sample_id");names(a)[-1]<-paste0(names(a)[-1],"_before");names(b)[-1]<-paste0(names(b)[-1],"_after")
q<-merge(a,b,by="group");q$retained_pct<-100*q$cells_after/q$cells_before
write.csv(q,"results/sc/review/sample_qc.csv",row.names=FALSE)
write.csv(summary_qc(m,"seurat_clusters"),"results/sc/review/cluster_qc.csv",row.names=FALSE)
panels<-list(B=c("Cd79a","Cd79b","Ms4a1","Cd19","Cd74","H2-Aa"),T=c("Cd3d","Cd3e","Trac","Cd4","Cd8a","Cd8b1","Tcf7","Foxp3","Il2ra","Rorc","Trdc","Trgc1"),NK=c("Ncr1","Klrd1","Nkg7","Gzmb","Xcl1"),Myeloid=c("Lyz2","Csf1r","Adgre1","Fcgr1","C1qa","C1qb","Ccr2","Ly6c2","S100a8","S100a9","Ly6g","Retnlg","Ltf","Ngp"),DC=c("Flt3","Zbtb46","Itgax","Xcr1","Clec9a","Clec10a","Cd209a","Siglech","Bst2","Tcf4"),Other=c("Jchain","Mzb1","Sdc1","Igha","Siglecf","Epx","Prg2","Ccr3","Mcpt8","Fcer1a","Ms4a2","Kit","Cpa3","Pf4","Ppbp","Hba-a1","Hbb-bs","Pecam1","Alb","Krt19","Mki67","Top2a","Isg15","Ifit3"))
genes<-intersect(unique(unlist(panels)),rownames(s))
counts<-GetAssayData(s,layer="counts")[genes,,drop=FALSE];norm<-GetAssayData(s,layer="data")[genes,,drop=FALSE]
levels_c<-levels(Idents(s));rows<-lapply(levels_c,function(cl){idx<-which(as.character(Idents(s))==cl);data.frame(cluster=cl,gene=genes,pct=100*Matrix::rowMeans(counts[,idx,drop=FALSE]>0),mean_log=Matrix::rowMeans(norm[,idx,drop=FALSE]))})
write.csv(do.call(rbind,rows),"results/sc/review/panel_evidence.csv",row.names=FALSE)
for(n in names(panels)){
 g<-intersect(panels[[n]],genes)
 p<-DotPlot(s,features=g)+RotatedAxis()+ggtitle(paste(n,"marker evidence"))
 ggsave(paste0("results/sc/review/markers_",n,".png"),p,width=max(8,length(g)*.48),height=9,dpi=140,bg="white")
}
# QC flags are screens, not doublet predictions or automatic exclusion criteria.
# High complexity: above the per-mouse median log10 UMI + 3 scaled MAD.
limits<-sapply(split(m,m$sample_id),function(x) median(log10(x$nCount_RNA+1))+3*mad(log10(x$nCount_RNA+1)))
m$high_complexity_screen<-log10(m$nCount_RNA+1)>limits[m$sample_id]
# Two independent lineage markers in each of two discordant programs.
detected<-function(gs)Matrix::colSums(counts[intersect(gs,rownames(counts)),,drop=FALSE]>0)>=2
T<-detected(c("Cd3d","Cd3e","Trac"));B<-detected(c("Cd79a","Cd79b","Ms4a1"));My<-detected(c("Csf1r","Adgre1","Fcgr1","Ly6g","Retnlg"))
m$discordant_lineage_screen<-(T&B)|(T&My)|(B&My)
write.csv(m[,c("cell_id","sample_id","condition","seurat_clusters","nCount_RNA","nFeature_RNA","percent.mt","high_complexity_screen","discordant_lineage_screen")],"results/sc/review/cell_qc_flags.csv",row.names=FALSE)
flag<-aggregate(cbind(high_complexity_screen,discordant_lineage_screen)~seurat_clusters,m,mean)
write.csv(flag,"results/sc/review/cluster_flags.csv",row.names=FALSE)
composition<-as.data.frame(table(m$sample_id,m$seurat_clusters));names(composition)<-c("sample_id","cluster","cells")
composition$within_sample_pct<-100*composition$cells/as.numeric(table(m$sample_id)[composition$sample_id])
write.csv(composition,"results/sc/review/sample_composition.csv",row.names=FALSE)
p<-ggplot(composition,aes(sample_id,within_sample_pct,fill=cluster))+geom_col()+labs(x="Mouse sample tag",y="Cells within sample (%)")+theme_bw()
ggsave("results/sc/review/sample_composition.png",p,width=10,height=6,dpi=150,bg="white")
writeLines(capture.output(sessionInfo()),"results/sc/review/sessionInfo.txt")
