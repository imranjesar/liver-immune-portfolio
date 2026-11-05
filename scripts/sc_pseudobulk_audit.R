suppressPackageStartupMessages(library(ggplot2))
p<-"results/sc/pseudobulk";x<-readRDS(file.path(p,"aggregates.rds"))
genes<-rownames(x$counts)
sets<-list(mitochondrial=grepl("^mt-",genes),ribosomal_protein=grepl("^Rp[sl][0-9]",genes),immunoglobulin_variable=grepl("^Ig[hkl]v",genes),Gm42418_unresolved=genes=="Gm42418")
a<-do.call(rbind,lapply(names(sets),function(n)data.frame(x$samples,feature_class=n,umi_fraction=colSums(x$counts[sets[[n]],,drop=FALSE])/pmax(colSums(x$counts),1))))
write.csv(a,file.path(p,"feature_class_fractions.csv"),row.names=FALSE)
keep<-a$condition%in%c("aclf_vehicle","aclf_ccs1477") & a$population%in%read.csv(file.path(p,"population_eligibility.csv"))$population[read.csv(file.path(p,"population_eligibility.csv"))$tested]
z<-a[keep,];z$condition<-factor(z$condition,levels=c("aclf_vehicle","aclf_ccs1477"))
p1<-ggplot(z,aes(sample_id,100*umi_fraction,fill=condition))+geom_col()+facet_grid(feature_class~population,scales="free_y")+theme_bw()+theme(axis.text.x=element_text(angle=90,vjust=.5,hjust=1))+labs(x="Mouse tag",y="Raw pseudobulk UMIs (%)",fill="Condition",title="Shared feature-class shifts require technical and biological review")
ggsave(file.path(p,"feature_class_audit.png"),p1,width=16,height=10,dpi=150,bg="white")
r<-read.csv(file.path(p,"all_tested_populations_de.csv"));r$feature_class<-ifelse(grepl("^mt-",r$gene),"mitochondrial",ifelse(grepl("^Rp[sl][0-9]",r$gene),"ribosomal_protein",ifelse(grepl("^Ig[hkl]v",r$gene),"immunoglobulin_variable",ifelse(r$gene=="Gm42418","Gm42418_unresolved","other"))))
write.csv(r[,c("population","gene","feature_class","log2FoldChange","padj_across_populations")],file.path(p,"gene_feature_flags.csv"),row.names=FALSE)
sig<-r[!is.na(r$padj_across_populations)&r$padj_across_populations<.05,]
shared<-aggregate(population~gene,sig,function(v)length(unique(v)));names(shared)[2]<-"significant_populations"
shared<-shared[shared$significant_populations>=3,]
shared<-shared[order(-shared$significant_populations,shared$gene),]
write.csv(shared,file.path(p,"shared_model_associations.csv"),row.names=FALSE)
