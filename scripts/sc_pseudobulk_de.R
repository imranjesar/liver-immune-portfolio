suppressPackageStartupMessages({library(DESeq2);library(ggplot2);library(jsonlite)})
cfg<-jsonlite::fromJSON("config/pseudobulk.json")
root<-"results/sc/pseudobulk";stopifnot(file.exists(file.path(root,"AGGREGATION_COMPLETE")))
if(file.exists(file.path(root,"DE_COMPLETE")))unlink(file.path(root,"DE_COMPLETE"))
x<-readRDS(file.path(root,"aggregates.rds"));set.seed(cfg$seed)
all_results<-list();coverage<-list()
for(pop in x$populations$population){
 out<-file.path(root,pop);dir.create(out,showWarnings=FALSE)
 m<-x$samples[x$samples$population==pop & x$samples$condition%in%c(cfg$reference,cfg$treatment),]
 m<-m[order(m$condition,m$sample_id),];rownames(m)<-m$aggregate_id
 stopifnot(nrow(m)==6,all(table(m$condition)==3),!anyDuplicated(m$sample_id))
 eligible<-all(m$n_cells>=cfg$min_cells_each_mouse & m$total_umi>=cfg$min_umi_each_mouse)
 coverage[[pop]]<-data.frame(population=pop,tested=eligible,min_cells=min(m$n_cells),min_total_umi=min(m$total_umi),reason=ifelse(eligible,"six_mice_pass","insufficient_cells_or_umi_in_at_least_one_mouse"))
 write.csv(m,file.path(out,"comparison_samples.csv"),row.names=FALSE)
 if(!eligible){cat("Skipping sparse population",pop,"\n");next}
 cat("Fitting",pop,"\n")
 cts<-x$counts[,m$aggregate_id,drop=FALSE]
 keep<-rowSums(cts>=cfg$min_gene_count)>=cfg$min_gene_samples
 cts<-cts[keep,,drop=FALSE];stopifnot(nrow(cts)>100)
 m$condition<-factor(m$condition,levels=c(cfg$reference,cfg$treatment))
 # Condition and capture are completely confounded: no estimable capture-adjusted contrast.
 design_check<-model.matrix(~condition+capture,m)
 write.csv(data.frame(columns=ncol(design_check),rank=qr(design_check)$rank,condition_capture_confounded=qr(design_check)$rank<ncol(design_check)),file.path(out,"design_audit.csv"),row.names=FALSE)
 stopifnot(qr(design_check)$rank<ncol(design_check))
 dds<-DESeqDataSetFromMatrix(cts,m,design=~condition)
 dds<-DESeq(dds,sfType=cfg$size_factors,minReplicatesForReplace=Inf,quiet=TRUE)
 contrast<-c("condition",cfg$treatment,cfg$reference)
 res<-results(dds,contrast=contrast,alpha=cfg$alpha)
 result<-data.frame(population=pop,gene=rownames(res),as.data.frame(res),row.names=NULL)
 result$wald_converged<-mcols(dds)$betaConv
 result$max_cooks<-apply(assays(dds)[["cooks"]],1,max,na.rm=TRUE)
 # Failed coefficient convergence is not treated as an ordinary successful test.
 result$pvalue[!result$wald_converged]<-NA_real_;result$padj[!result$wald_converged]<-NA_real_
 write.csv(result,file.path(out,"primary_de.csv"),row.names=FALSE)
 write.csv(data.frame(gene=rownames(dds),counts(dds,normalized=TRUE),check.names=FALSE),file.path(out,"normalized_counts.csv"),row.names=FALSE)
 write.csv(data.frame(sample_id=m$sample_id,aggregate_id=m$aggregate_id,n_cells=m$n_cells,total_umi=m$total_umi,size_factor=sizeFactors(dds)),file.path(out,"size_factors.csv"),row.names=FALSE)
 saveRDS(dds,file.path(out,"deseq2.rds"))
 # Blind VST is for sample QC only; all testing above uses integer raw sums.
 vsd<-varianceStabilizingTransformation(dds,blind=TRUE)
 pc<-plotPCA(vsd,intgroup="condition",returnData=TRUE)
 pc$sample_id<-m$sample_id[match(pc$name,m$aggregate_id)]
 variance<-round(100*attr(pc,"percentVar"),1)
 p<-ggplot(pc,aes(PC1,PC2,color=condition,label=sample_id))+geom_point(size=3)+geom_text(vjust=-.6,show.legend=FALSE)+theme_bw()+labs(title=paste(pop,"— mouse-level PCA"),x=paste0("PC1 (",variance[1],"%)"),y=paste0("PC2 (",variance[2],"%)"))
 ggsave(file.path(out,"pca.png"),p,width=8,height=5,dpi=150,bg="white")
 write.csv(pc,file.path(out,"pca_coordinates.csv"),row.names=FALSE)
 distmat<-as.matrix(dist(t(assay(vsd))));dimnames(distmat)<-list(m$sample_id,m$sample_id)
 distdf<-as.data.frame(as.table(distmat));names(distdf)<-c("mouse_1","mouse_2","distance")
 p<-ggplot(distdf,aes(mouse_1,mouse_2,fill=distance))+geom_tile()+coord_equal()+theme_bw()+labs(title=paste(pop,"VST sample distances"),x=NULL,y=NULL)
 ggsave(file.path(out,"sample_distances.png"),p,width=6,height=5,dpi=150,bg="white")
 if(isTRUE(cfg$leave_one_mouse_out)){
  loo<-matrix(NA_real_,nrow=nrow(cts),ncol=6,dimnames=list(rownames(cts),m$sample_id))
  for(j in seq_len(nrow(m))){
   cat("  Leaving out",m$sample_id[j],"\n")
   sub<-DESeqDataSetFromMatrix(cts[,-j,drop=FALSE],m[-j,,drop=FALSE],design=~condition)
   sub<-DESeq(sub,sfType=cfg$size_factors,minReplicatesForReplace=Inf,quiet=TRUE)
   sr<-results(sub,contrast=contrast,alpha=cfg$alpha)
   good<-mcols(sub)$betaConv
   loo[good,j]<-sr$log2FoldChange[good]
  }
  write.csv(data.frame(gene=rownames(loo),loo,check.names=FALSE),file.path(out,"leave_one_mouse_out_lfc.csv"),row.names=FALSE)
  base<-result$log2FoldChange[match(rownames(loo),result$gene)]
  result$loo_valid_fits<-rowSums(is.finite(loo))
  result$loo_same_direction_fits<-rowSums(sign(loo)==sign(base),na.rm=TRUE)
  result$loo_lfc_min<-apply(loo,1,function(z)if(any(is.finite(z)))min(z,na.rm=TRUE)else NA_real_)
  result$loo_lfc_max<-apply(loo,1,function(z)if(any(is.finite(z)))max(z,na.rm=TRUE)else NA_real_)
 }
 all_results[[pop]]<-result
 write.csv(result,file.path(out,"primary_de_with_loo.csv"),row.names=FALSE)
 writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
}
write.csv(do.call(rbind,coverage),file.path(root,"population_eligibility.csv"),row.names=FALSE)
combined<-do.call(rbind,all_results)
combined$padj_across_populations<-p.adjust(combined$pvalue,method="BH")
write.csv(combined,file.path(root,"all_tested_populations_de.csv"),row.names=FALSE)
summary<-do.call(rbind,lapply(split(combined,combined$population),function(z)data.frame(population=z$population[1],genes_prefiltered=nrow(z),nonmissing_p=sum(!is.na(z$pvalue)),within_population_FDR05=sum(z$padj<cfg$alpha,na.rm=TRUE),across_population_FDR05=sum(z$padj_across_populations<cfg$alpha,na.rm=TRUE),significant_stable_direction_all6=sum(z$padj_across_populations<cfg$alpha & z$loo_valid_fits==6 & z$loo_same_direction_fits==6,na.rm=TRUE))))
write.csv(summary,file.path(root,"de_summary.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(root,"de_sessionInfo.txt"))
writeLines("COMPLETE",file.path(root,"DE_COMPLETE"))
