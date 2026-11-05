source("scripts/sc_sensitivity_helpers.R")
x<-load_strict_inputs()
stopifnot(ncol(x$counts)==46954,length(unique(x$meta$sample_id))==9)
rerun_sc(x$counts,x$meta,"results/sc/strict500",markers=FALSE)
