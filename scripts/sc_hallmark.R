suppressPackageStartupMessages({library(fgsea);library(jsonlite)})
set.seed(20261103)
out<-"results/sc/pathways";dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(out,"COMPLETE")))unlink(file.path(out,"COMPLETE"))
gmt<-"data/annotations/mh.all.v2026.2.Mm.symbols.gmt"
source_info<-jsonlite::fromJSON("metadata/hallmark_source.json")
stopifnot(identical(digest::digest(file=gmt,algo="sha256"),source_info$sha256))
sets<-fgsea::gmtPathways(gmt);stopifnot(length(sets)==50)
r<-read.csv("results/sc/pseudobulk/all_tested_populations_de.csv",stringsAsFactors=FALSE)
all<-list();coverage<-list()
for(pop in sort(unique(r$population)))for(variant in c("all_valid_features","exclude_flagged_features")){
 x<-r[r$population==pop & is.finite(r$stat) & is.finite(r$pvalue) & r$wald_converged,]
 if(variant=="exclude_flagged_features")x<-x[!grepl("^mt-|^Rp[sl][0-9]|^Ig[hkl]v",x$gene)&x$gene!="Gm42418",]
 stopifnot(!anyDuplicated(x$gene))
 x<-x[order(-x$stat,x$gene),];ranks<-setNames(x$stat,x$gene)
 membership<-data.frame(pathway=names(sets),overlap=vapply(sets,function(g)sum(g%in%names(ranks)),integer(1)))
 membership$population<-pop;membership$variant<-variant
 coverage[[paste(pop,variant)]]<-membership
 set.seed(20261103)
 z<-fgseaMultilevel(pathways=sets,stats=ranks,minSize=15,maxSize=500,eps=0,nPermSimple=10000,nproc=1)
 z<-as.data.frame(z);z$population<-pop;z$variant<-variant
 z$leadingEdge<-vapply(z$leadingEdge,paste,collapse=";",FUN.VALUE=character(1))
 all[[paste(pop,variant)]]<-z
 write.csv(data.frame(gene=names(ranks),wald_stat=as.numeric(ranks)),file.path(out,paste0(pop,"_",variant,"_ranks.csv")),row.names=FALSE)
}
z<-do.call(rbind,all)
z$padj_across_populations<-NA_real_
for(v in unique(z$variant)){i<-which(z$variant==v);z$padj_across_populations[i]<-p.adjust(z$pval[i],"BH")}
write.csv(z,file.path(out,"hallmark_all_results.csv"),row.names=FALSE)
write.csv(do.call(rbind,coverage),file.path(out,"gene_set_coverage.csv"),row.names=FALSE)
a<-z[z$variant=="all_valid_features",];b<-z[z$variant=="exclude_flagged_features",]
comp<-merge(a[,c("population","pathway","NES","padj_across_populations")],b[,c("population","pathway","NES","padj_across_populations")],by=c("population","pathway"),suffixes=c("_primary","_sensitivity"))
comp$same_direction<-sign(comp$NES_primary)==sign(comp$NES_sensitivity)
comp$both_q05<-comp$padj_across_populations_primary<.05 & comp$padj_across_populations_sensitivity<.05
write.csv(comp,file.path(out,"rank_feature_sensitivity.csv"),row.names=FALSE)
jsonlite::write_json(list(gene_sets=source_info,ranking="signed DESeq2 Wald statistic; finite p and statistic; converged model",min_size=15,max_size=500,nPermSimple=10000,eps=0,seed=20261103,fgsea=as.character(packageVersion("fgsea")),multiple_testing="BH across all populations separately within each of two rank variants",interpretation="Rank-feature sensitivity is not a batch correction or a DE model refit"),file.path(out,"parameters.json"),auto_unbox=TRUE,pretty=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
writeLines("COMPLETE",file.path(out,"COMPLETE"))
