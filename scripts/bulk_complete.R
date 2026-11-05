suppressPackageStartupMessages({library(DESeq2);library(AnnotationDbi);library(org.Mm.eg.db);library(fgsea);library(jsonlite)})
set.seed(20261103)
out <- 'results/bulk/curated';dir.create(out,recursive=TRUE,showWarnings=FALSE)
unlink(file.path(out,'COMPLETE'))
write_table<-function(x,n)write.csv(x,file.path(out,paste0(n,'.csv')),row.names=FALSE)
x<-read.table(gzfile('data/raw/GSE247373_countsAll.txt.gz'),header=TRUE,check.names=FALSE)
cts<-as.matrix(x[,-1]);rownames(cts)<-x[[1]]
m<-read.delim('config/bulk_samples.tsv');stopifnot(setequal(colnames(cts),m$count_column),!anyDuplicated(rownames(cts)),all(cts>=0),all(cts==round(cts)))
m<-m[match(colnames(cts),m$count_column),];rownames(m)<-m$count_column
# Preserve original versioned IDs; never merge distinct count rows by symbol.
ids<-sub('\\.[0-9]+$','',rownames(cts))
keys_valid<-intersect(unique(ids),keys(org.Mm.eg.db,keytype='ENSEMBL'))
map<-AnnotationDbi::select(org.Mm.eg.db,keys=keys_valid,keytype='ENSEMBL',columns='SYMBOL')
map<-unique(map[!is.na(map$SYMBOL),]);write_table(map,'annotation_all_mappings')
maps<-split(map$SYMBOL,map$ENSEMBL)
ann<-data.frame(gene_id=rownames(cts),ensembl=ids,symbol=vapply(ids,function(k){v<-unique(maps[[k]]);if(length(v)==1)v else NA_character_},character(1)),symbol_count=vapply(ids,function(k)length(unique(maps[[k]])),integer(1)))
write_table(ann,'gene_annotation')
# All 12 libraries are descriptive QC; only nine liver libraries enter inference.
qc<-cbind(m,total_gene_counts=colSums(cts),genes_ge10=colSums(cts>=10),zero_fraction=colMeans(cts==0))
mt<-!is.na(ann$symbol)&grepl('^mt-',ann$symbol)
qc$mitochondrial_fraction<-colSums(cts[mt,,drop=FALSE])/colSums(cts)
write_table(qc,'library_qc')
keep<-m$condition!='spleen_0h';meta<-m[keep,];meta$condition<-factor(meta$condition,c('liver_24h','liver_72h'))
y<-cts[,keep];filter<-rowSums(y>=10)>=3
d<-DESeqDataSetFromMatrix(y[filter,],meta,~condition);d<-DESeq(d,quiet=TRUE,minReplicatesForReplace=Inf)
r<-as.data.frame(results(d,contrast=c('condition','liver_72h','liver_24h'),alpha=.05))
r<-cbind(ann[match(rownames(r),ann$gene_id),],r);r$converged<-mcols(d)$betaConv
r$max_cooks<-apply(assays(d)[['cooks']],1,max,na.rm=TRUE)
# Leave-one-library-out checks use the fixed primary gene universe. They are not new discovery tests.
loo<-list()
for(i in seq_len(ncol(d))){
 message('Leave-one-library-out ',i,'/',ncol(d))
 z<-DESeqDataSetFromMatrix(y[filter,-i,drop=FALSE],meta[-i,,drop=FALSE],~condition)
 z<-DESeq(z,quiet=TRUE,minReplicatesForReplace=Inf)
 q<-as.data.frame(results(z,contrast=c('condition','liver_72h','liver_24h'),alpha=.05))
 loo[[i]]<-data.frame(gene_id=rownames(q),omitted=colnames(d)[i],log2FoldChange=q$log2FoldChange,converged=mcols(z)$betaConv)
}
l<-do.call(rbind,loo);write_table(l,'leave_one_out')
lm<-matrix(l$log2FoldChange,nrow=nrow(r));cm<-matrix(l$converged,nrow=nrow(r))
r$loo_all_converged<-rowSums(!cm)==0
r$loo_same_direction<-rowSums(sign(lm)==sign(r$log2FoldChange),na.rm=TRUE)==ncol(d)
r$loo_min_lfc<-apply(lm,1,min);r$loo_max_lfc<-apply(lm,1,max)
write_table(r,'differential_expression')
qc2<-cbind(meta,size_factor=sizeFactors(d),max_cooks=apply(assays(d)[['cooks']],2,max,na.rm=TRUE));write_table(qc2,'model_samples')
v<-vst(d,blind=TRUE);pc<-plotPCA(v,intgroup='condition',returnData=TRUE)
pv<-attr(pc,'percentVar');write_table(pc,'pca');write_table(data.frame(component=1:2,variance=pv),'pca_variance')
write_table(cbind(sample=colnames(v),as.data.frame(as.matrix(dist(t(assay(v)))))), 'sample_distances')
# Symbol-level ranks: exclude ambiguous mappings; among duplicate symbols retain the
# highest primary baseMean (then original ID), rather than selecting by significance.
rankx<-r[is.finite(r$stat)&is.finite(r$pvalue)&r$converged&!is.na(r$symbol),]
rankx<-rankx[order(-rankx$baseMean,rankx$gene_id),];rankx<-rankx[!duplicated(rankx$symbol),]
write_table(rankx,'rank_feature_selection')
sets<-gmtPathways('data/annotations/mh.all.v2026.2.Mm.symbols.gmt')
src<-fromJSON('metadata/hallmark_source.json');stopifnot(digest::digest(file='data/annotations/mh.all.v2026.2.Mm.symbols.gmt',algo='sha256')==src$sha256)
path<-list()
for(variant in c('all_valid_features','exclude_flagged_features')){
 a<-rankx
 if(variant=='exclude_flagged_features')a<-a[!grepl('^mt-|^Rp[sl][0-9]|^Ig[hkl]v',a$symbol)&a$symbol!='Gm42418',]
 a<-a[order(-a$stat,a$symbol),];ranks<-setNames(a$stat,a$symbol)
 write_table(data.frame(symbol=names(ranks),stat=as.numeric(ranks)),paste0('ranks_',variant))
 set.seed(20261103);p<-as.data.frame(fgseaMultilevel(sets,ranks,minSize=15,maxSize=500,eps=0,nPermSimple=10000,BPPARAM=BiocParallel::SerialParam()))
 p$leadingEdge<-vapply(p$leadingEdge,paste,collapse=';',FUN.VALUE=character(1));p$variant<-variant;path[[variant]]<-p
}
p<-do.call(rbind,path);write_table(p,'hallmark')
# Descriptive cross-study context: contrasts and cell populations are different.
s<-read.csv('results/sc/pathways/hallmark_all_results.csv');s<-s[s$variant=='all_valid_features',c('population','pathway','NES','padj_across_populations')]
b<-p[p$variant=='all_valid_features',c('pathway','NES','padj')]
context<-merge(s,b,by='pathway',suffixes=c('_sc','_bulk'));context$same_direction<-sign(context$NES_sc)==sign(context$NES_bulk)
write_table(context,'cross_study_context')
saveRDS(d,file.path(out,'deseq2.rds'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
write_json(list(contrast='liver 72 h minus liver 24 h',sample_units='GEO libraries; individual-mouse provenance/pooling not established from retrieved metadata',n_24h=4,n_72h=5,seed=20261103,gene_filter='count >=10 in >=3 of nine libraries',replacement='disabled',annotation=as.character(packageVersion('org.Mm.eg.db')),rank_duplicate_policy='highest baseMean then original ID; one unambiguous symbol per feature',fdr='DE gene BH; Hallmark BH separately per rank variant',pathway_min=15,pathway_max=500),file.path(out,'parameters.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines('COMPLETE',file.path(out,'COMPLETE'))
print(c(genes=nrow(r),de=sum(r$padj<.05,na.rm=TRUE),stable=sum(r$padj<.05&r$loo_same_direction&r$loo_all_converged,na.rm=TRUE),hallmark=sum(p$variant=='all_valid_features'&p$padj<.05)))
