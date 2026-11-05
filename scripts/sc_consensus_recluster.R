source("scripts/sc_sensitivity_helpers.R")
stopifnot(file.exists("results/sc/doublets/COMPLETE"))
x<-load_strict_inputs();d<-read.csv("results/sc/doublets/cell_doublet_summary.csv",stringsAsFactors=FALSE)
i<-match(rownames(x$meta),d$cell_id);stopifnot(!anyNA(i),setequal(d$cell_id,rownames(x$meta)))
x$meta$doublet_votes<-d$doublet_votes[i];x$meta$rate_sensitive_call<-d$rate_sensitive_call[i]
keep<-!d$consensus_doublet[i]
stopifnot(length(unique(x$meta$sample_id[keep]))==9)
rerun_sc(x$counts[,keep,drop=FALSE],x$meta[keep,,drop=FALSE],"results/sc/strict500_consensus",markers=TRUE)
