# Run a script with the project-local library populated by renv::restore().
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)stop("Usage: Rscript scripts/run_r.R scripts/analysis.R")
if(!requireNamespace("renv",quietly=TRUE))stop("Install renv and restore renv.lock first")
renv::load(project=getwd())
source(args[[1]],echo=FALSE)
