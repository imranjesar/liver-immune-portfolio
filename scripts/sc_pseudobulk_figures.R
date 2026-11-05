suppressPackageStartupMessages({library(ggplot2);library(DESeq2)})
root<-"results/sc/pseudobulk";stopifnot(file.exists(file.path(root,"DE_COMPLETE")))
r<-read.csv(file.path(root,"all_tested_populations_de.csv"),stringsAsFactors=FALSE)
r$evidence<-ifelse(!is.na(r$padj_across_populations)&r$padj_across_populations<.05,"Across-population BH < 0.05","Other/untested")
r$minus_log10_global<- -log10(pmax(r$padj_across_populations,1e-300))
p<-ggplot(r,aes(log2FoldChange,minus_log10_global,color=evidence))+geom_point(size=.45,alpha=.6,na.rm=TRUE)+facet_wrap(~population,scales="free",ncol=3)+theme_bw()+theme(legend.position="bottom")+scale_color_manual(values=c("Across-population BH < 0.05"="#377eb8","Other/untested"="grey75"))+labs(x="Unshrunk log2 fold change (CCS1477 / vehicle)",y="−log10 across-population BH",color="Model-based exploratory evidence",title="Mouse-level pseudobulk — condition and library are confounded")
ggsave(file.path(root,"exploratory_volcano_panels.png"),p,width=14,height=9,dpi=160,bg="white")
summary<-read.csv(file.path(root,"de_summary.csv"))
p<-ggplot(summary,aes(reorder(population,across_population_FDR05),across_population_FDR05))+geom_col(fill="#377eb8")+coord_flip()+theme_bw()+labs(x=NULL,y="Gene–population tests with across-population BH < 0.05",title="Exploratory model counts; not measures of treatment efficacy")
ggsave(file.path(root,"model_evidence_counts.png"),p,width=10,height=6,dpi=150,bg="white")
for(pop in unique(r$population)){
 out<-file.path(root,pop);z<-r[r$population==pop & !is.na(r$padj_across_populations),]
 z<-z[order(z$padj_across_populations),];top<-head(z,20)
 write.csv(top,file.path(out,"top20_exploratory_genes.csv"),row.names=FALSE)
 d<-readRDS(file.path(out,"deseq2.rds"));v<-varianceStabilizingTransformation(d,blind=TRUE)
 mat<-assay(v)[top$gene,,drop=FALSE];mat<-t(scale(t(mat)));mat[!is.finite(mat)]<-0
 dimnames(mat)<-list(top$gene,as.character(colData(d)$sample_id))
 dat<-as.data.frame(as.table(mat));names(dat)<-c("gene","mouse","row_z")
 dat$gene<-factor(dat$gene,levels=rev(top$gene))
 p<-ggplot(dat,aes(mouse,gene,fill=row_z))+geom_tile()+theme_bw()+scale_fill_gradient2(low="#2166ac",mid="white",high="#b2182b",midpoint=0)+labs(x="Mouse tag (10–12 CCS1477; 04–06 vehicle)",y=NULL,fill="Row z score",title=paste(pop,"— top model-ranked genes"))
 ggsave(file.path(out,"top20_heatmap.png"),p,width=8,height=8,dpi=150,bg="white")
}
