r <- read.csv("results/bulk/72h_vs_24h.csv",row.names=1)
stopifnot(nrow(r)>1000,all(r$padj[!is.na(r$padj)]>=0),all(r$padj[!is.na(r$padj)]<=1))
m <- read.csv("results/bulk/samples_used.csv")
stopifnot(sum(m$condition=="liver_24h")==4,sum(m$condition=="liver_72h")==5)
stopifnot(file.exists("results/sc/seurat.rds"))
{
  suppressPackageStartupMessages(library(Seurat))
  s <- readRDS("results/sc/seurat.rds")
  stopifnot(length(unique(s$sample_id))==9,!anyDuplicated(colnames(s)),all(is.finite(Embeddings(s,"umap"))))
  stopifnot(all(s$nFeature_RNA>=200),all(s$nFeature_RNA<=6000),all(s$percent.mt<=20))
}
cat("Available analysis result integrity checks passed\n")
