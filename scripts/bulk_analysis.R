suppressPackageStartupMessages({library(DESeq2);library(ggplot2)})
set.seed(20261103)
dir.create("results/bulk",recursive=TRUE,showWarnings=FALSE)
x <- read.table(gzfile("data/raw/GSE247373_countsAll.txt.gz"),header=TRUE,check.names=FALSE)
# File-specific mapping is validated against the observed header in config/bulk_samples.tsv.
meta <- read.delim("config/bulk_samples.tsv",stringsAsFactors=FALSE)
stopifnot(!anyDuplicated(x[[1]]),setequal(names(x)[-1],meta$count_column))
counts <- as.matrix(x[,-1]);rownames(counts)<-x[[1]]
stopifnot(all(is.finite(counts)),all(counts>=0),all(counts==round(counts)))
meta <- meta[match(colnames(counts),meta$count_column),]
rownames(meta)<-meta$count_column
# Primary contrast stays within liver; spleen remains an exploratory tissue-confounded comparison.
keep <- meta$condition %in% c("liver_24h","liver_72h")
meta<-meta[keep,];counts<-counts[,keep]
meta$condition<-factor(meta$condition,levels=c("liver_24h","liver_72h"))
dds<-DESeqDataSetFromMatrix(counts,meta,design=~condition)
dds<-dds[rowSums(counts(dds)>=10)>=3,]
dds<-DESeq(dds,quiet=TRUE)
r<-results(dds,contrast=c("condition","liver_72h","liver_24h"),alpha=0.05)
write.csv(as.data.frame(r[order(r$padj),]),"results/bulk/72h_vs_24h.csv")
v<-vst(dds,blind=TRUE)
pca_plot <- plotPCA(v,intgroup="condition")
pdf("results/bulk/pca.pdf");print(pca_plot);dev.off()
ggsave("results/bulk/pca.png",pca_plot,width=7,height=5,dpi=160)
write.csv(plotPCA(v,intgroup="condition",returnData=TRUE),"results/bulk/pca_coordinates.csv",row.names=FALSE)
pdf("results/bulk/MA.pdf");plotMA(r,ylim=c(-5,5));dev.off()
write.csv(as.data.frame(colData(dds)),"results/bulk/samples_used.csv")
saveRDS(dds,"results/bulk/deseq2.rds")
writeLines(capture.output(sessionInfo()),"results/bulk/sessionInfo.txt")
