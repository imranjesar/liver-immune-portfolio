suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ggrepel);library(ragg);library(svglite);library(fgsea)})
set.seed(20261103)
src<-"figures/source_data";pub<-"figures/publication";preview<-"figures/previews"
dir.create(pub,recursive=TRUE,showWarnings=FALSE);dir.create(preview,recursive=TRUE,showWarnings=FALSE)
read<-function(n)read.csv(if(grepl("[.]gz$",n))gzfile(file.path(src,n))else file.path(src,n),stringsAsFactors=FALSE)
font<-"Liberation Sans"
theme_pub<-function(){theme_classic(base_size=7.5,base_family=font)+theme(axis.title=element_text(size=7.5),axis.text=element_text(size=6.7,color="black"),axis.line=element_line(linewidth=.25),axis.ticks=element_line(linewidth=.25),axis.ticks.length=grid::unit(1,"mm"),plot.title=element_text(size=8.5,face="bold",hjust=0),plot.subtitle=element_text(size=7),plot.tag=element_text(size=10,face="bold"),strip.background=element_blank(),strip.text=element_text(size=7.2,face="bold"),legend.title=element_text(size=7),legend.text=element_text(size=6.5),legend.key.height=grid::unit(2.8,"mm"),legend.key.width=grid::unit(3.5,"mm"),legend.spacing=grid::unit(1.5,"mm"),plot.margin=margin(3,3,3,3),panel.spacing=grid::unit(3.5,"mm"),plot.background=element_rect(fill="white",color=NA))}
theme_set(theme_pub())
condition_colors<-c(aclf_vehicle="#777777",aclf_ccs1477="#0072B2")
condition_labels<-c(aclf_vehicle="Vehicle",aclf_ccs1477="CCS1477")
poporder<-c("B_cells","CD4_T_cells","CD8_T_cells","NK_cells","pDC","Antibody_secreting_cells")
popnames<-c(B_cells="B cells",CD4_T_cells="CD4 T",CD8_T_cells="CD8 T",NK_cells="NK",pDC="pDC",Antibody_secreting_cells="ASC")
export<-function(p,name,w,h){
 ggsave(file.path(pub,paste0(name,".pdf")),p,device=grDevices::cairo_pdf,width=w,height=h,units="mm",family=font,bg="white",limitsize=FALSE)
 ggsave(file.path(pub,paste0(name,".svg")),p,device=svglite::svglite,width=w,height=h,units="mm",bg="white",limitsize=FALSE)
 ggsave(file.path(pub,paste0(name,".png")),p,device=ragg::agg_png,width=w,height=h,units="mm",dpi=600,bg="white",limitsize=FALSE)
 ggsave(file.path(preview,paste0(name,".png")),p,device=ragg::agg_png,width=w,height=h,units="mm",dpi=240,bg="white",limitsize=FALSE)
}
# Figure 1: direct numeric labels connect the atlas to the marker-validation key.
key<-read("atlas_lineage_key.csv");coords<-read("atlas_coordinates.csv.gz");dots<-read("atlas_marker_evidence.csv")
palette<-setNames(c("#D89000","#CC79A7","#7B61A8","#0072B2","#009E73","#56B4E9","#9B6D37","#648C16","#C44E52","#E69F00","#8C564B","#A8A8A8","#D1D1D1"),key$broad_type)
coords<-coords[sample(nrow(coords)),]
centers<-do.call(rbind,lapply(split(coords,coords$type_id),function(z){cx<-median(z$UMAP1);cy<-median(z$UMAP2);z[which.min((z$UMAP1-cx)^2+(z$UMAP2-cy)^2),]}))
a<-ggplot(coords,aes(UMAP1,UMAP2,color=broad_type))+geom_point(size=.12,alpha=.75,stroke=0)+scale_color_manual(values=palette,guide="none")+geom_label_repel(data=centers,aes(label=type_id),color="black",fill="white",size=2.3,label.size=.15,label.padding=grid::unit(.6,"mm"),box.padding=.25,point.padding=.1,min.segment.length=0,segment.color="grey50",segment.size=.2,seed=20261103,max.overlaps=Inf)+coord_equal()+labs(x="UMAP 1",y="UMAP 2",title="Immune-cell landscape")
key$label<-paste0(key$type_id,"  ",key$display)
dots$lineage<-factor(key$label[match(dots$broad_type,key$broad_type)],levels=rev(key$label))
dots$gene<-factor(dots$gene,levels=unique(dots$gene))
b<-ggplot(dots,aes(gene,lineage))+geom_point(aes(size=pct,fill=mean_log_z),shape=21,color="grey60",stroke=.15)+scale_size_continuous(range=c(.1,2),limits=c(0,100),breaks=c(0,50,100),name="Detected (%)")+scale_fill_gradient(low="#EEEEEE",high="#2166AC",limits=c(-1.5,1.5),oob=scales::squish,breaks=c(-1,0,1),name="Mean expression\n(gene-wise z)")+labs(x=NULL,y=NULL,title="Lineage-marker support")+theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1,face="italic",size=6.3),axis.text.y=element_text(size=6.2),legend.position="bottom",legend.box="horizontal",legend.key.width=grid::unit(3,"mm"),legend.title=element_text(size=6.2),legend.text=element_text(size=6),axis.ticks.y=element_blank())+guides(fill=guide_colorbar(display="rectangles",barwidth=grid::unit(13,"mm"),barheight=grid::unit(2,"mm")),size=guide_legend(override.aes=list(fill="grey40")))
f1<-(a+b)+plot_layout(widths=c(1,1.12))+plot_annotation(tag_levels="A")
export(f1,"Figure_1_cell_atlas",180,113)
# Figure 2: each point is one mouse, not one cell.
pc<-read("pseudobulk_pca.csv");pv<-read("pca_variance.csv")
pc$population<-factor(pc$population,levels=poporder)
strips<-setNames(paste0(popnames[pv$population],"\nPC1 ",round(pv$PC1_percent),"%, PC2 ",round(pv$PC2_percent),"%"),pv$population)
pc$mouse<-sub("tag","",pc$sample_id)
a<-ggplot(pc,aes(PC1,PC2,color=condition))+geom_point(size=1.6)+geom_text_repel(aes(label=mouse),size=2,seed=20261103,box.padding=.2,point.padding=.1,min.segment.length=Inf,max.overlaps=Inf,show.legend=FALSE)+facet_wrap(~population,ncol=3,scales="free",labeller=as_labeller(strips))+scale_color_manual(values=condition_colors,labels=condition_labels,name=NULL)+labs(title="Mouse-level pseudobulk",x="PC1",y="PC2")+theme(legend.position="top",legend.justification="left")+scale_x_continuous(expand=expansion(mult=.22))+scale_y_continuous(expand=expansion(mult=.22))
effects<-read("illustrative_ifn_gene_effects.csv")
gene_order<-unique(effects$gene[order(effects$illustration_order)])
effects$gene<-factor(effects$gene,levels=gene_order);effects$population<-factor(popnames[effects$population],levels=rev(unname(popnames[poporder])))
effects$supported<-!is.na(effects$padj_across_populations)&effects$padj_across_populations<.05
b<-ggplot(effects,aes(log2FoldChange,population))+geom_vline(xintercept=0,linewidth=.25,color="grey65",linetype="dashed")+geom_segment(aes(x=log2FoldChange-1.96*lfcSE,xend=log2FoldChange+1.96*lfcSE,yend=population),linewidth=.35,color="grey35")+geom_point(aes(fill=supported),shape=21,size=1.6,stroke=.3,color="#174A7E")+scale_fill_manual(values=c(`TRUE`="#0072B2",`FALSE`="white"),labels=c(`TRUE`="q < 0.05",`FALSE`="q >= 0.05"),name="Across-population BH")+facet_wrap(~gene,ncol=3)+labs(title="Illustrative interferon-leading-edge genes",x="log2 fold change (CCS1477 / vehicle)",y=NULL)+theme(strip.text=element_text(face="italic",size=7.5),legend.position="bottom",legend.justification="left",axis.text.y=element_text(size=6.2))
f2<-(a/b)+plot_layout(heights=c(1.1,1))+plot_annotation(tag_levels="A")
export(f2,"Figure_2_mouse_level_expression",180,180)
# Figure 3: show every Hallmark set associated in any primary population (no top-N trimming).
h<-read("hallmark_all_results.csv");primary<-h[h$variant=="all_valid_features",];comp<-read("hallmark_rank_sensitivity.csv")
shown<-unique(primary$pathway[!is.na(primary$padj_across_populations)&primary$padj_across_populations<.05])
stopifnot(length(shown)>0)
order_paths<-shown[order(vapply(shown,function(p)min(primary$padj_across_populations[primary$pathway==p],na.rm=TRUE),numeric(1)),shown)]
pretty<-function(x){
 y<-tolower(gsub("_"," ",sub("HALLMARK_","",x)));y<-paste0(toupper(substr(y,1,1)),substr(y,2,nchar(y)))
 replacements<-c("Myc"="MYC"," v1"=" V1"," v2"=" V2","Il2 stat5"="IL2-STAT5","Mtorc1"="mTORC1","E2f"="E2F","Tnfa signaling via nfkb"="TNF-alpha signaling via NF-kappaB","Il6 jak stat3"="IL6-JAK-STAT3","G2m"="G2M","P53"="p53","Dna"="DNA","Kras"="KRAS")
 for(n in names(replacements))y<-gsub(n,replacements[[n]],y,fixed=TRUE)
 y
}
heat<-primary[primary$pathway%in%shown,];heat$pathway<-factor(heat$pathway,levels=rev(order_paths));heat$population<-factor(heat$population,levels=poporder)
heat$primary_q05<-heat$padj_across_populations<.05
ci<-match(paste(heat$population,heat$pathway),paste(comp$population,comp$pathway))
heat$stable<-comp$both_q05[ci]&comp$same_direction[ci]
write.csv(heat,file.path(src,"displayed_hallmark_heatmap.csv"),row.names=FALSE)
a<-ggplot(heat,aes(population,pathway))+geom_tile(aes(fill=NES),color="white",linewidth=.25)+geom_point(data=heat[heat$primary_q05,],shape=21,size=.7,fill="white",stroke=.2,color="black")+geom_point(data=heat[heat$stable,],shape=16,size=.65,color="black")+scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",midpoint=0,limits=c(-3,3),oob=scales::squish,breaks=c(-2,0,2),name="NES")+scale_x_discrete(labels=popnames)+scale_y_discrete(labels=pretty)+labs(x=NULL,y=NULL,title="Condition-associated Hallmark programs")+theme(axis.text.y=element_text(size=6.4),axis.text.x=element_text(size=7),axis.line=element_blank(),axis.ticks=element_blank(),legend.key.height=grid::unit(5,"mm"))
a<-a+guides(fill=guide_colorbar(display="rectangles"))
comp$population<-factor(comp$population,levels=poporder)
comp$robust<-comp$both_q05&comp$same_direction
popcols<-setNames(c("#D89000","#009E73","#0072B2","#7B61A8","#9B6D37","#C44E52"),poporder)
b<-ggplot(comp,aes(NES_primary,NES_sensitivity,color=population))+geom_abline(slope=1,intercept=0,linewidth=.3,color="grey65")+geom_hline(yintercept=0,linewidth=.2,color="grey85")+geom_vline(xintercept=0,linewidth=.2,color="grey85")+geom_point(aes(shape=robust),size=1.25,alpha=.7,stroke=.4)+scale_shape_manual(values=c(`TRUE`=16,`FALSE`=1),guide="none")+scale_color_manual(values=popcols,labels=popnames,name=NULL)+coord_equal()+labs(x="Primary NES",y="Feature-filtered NES",title="Rank-feature sensitivity")+theme(legend.position="bottom",legend.key.width=grid::unit(2.5,"mm"),legend.text=element_text(size=6))+guides(color=guide_legend(nrow=2,override.aes=list(alpha=1,size=1.5)))
choice<-read("illustrative_enrichment_choice.csv");rankdf<-read("illustrative_enrichment_ranks.csv");ranks<-setNames(rankdf$wald_stat,rankdf$gene)
sets<-fgsea::gmtPathways("data/annotations/mh.all.v2026.2.Mm.symbols.gmt")
c<-fgsea::plotEnrichment(sets[[choice$pathway]],ranks)+theme_pub()+labs(title=paste0(popnames[choice$population],": ",pretty(choice$pathway)),x="Gene rank (signed Wald statistic)",y="Running enrichment score")+theme(plot.title=element_text(size=7.4))
c<-c+labs(subtitle=paste0("NES = ",formatC(choice$NES,digits=2,format="f"),"; global q = ",formatC(choice$padj_across_populations,digits=1,format="e")))+theme(plot.subtitle=element_text(size=6.5))
c$layers[[1]]$aes_params$colour<-"#0072B2"
c$layers[[1]]$aes_params$linewidth<-.5
for(i in c(3,4)){c$layers[[i]]$aes_params$colour<-"grey55";c$layers[[i]]$aes_params$linewidth<-.2}
f3<-(a/(b+c))+plot_layout(heights=c(1.5,1))+plot_annotation(tag_levels="A")
export(f3,"Figure_3_hallmark_programs",180,178)
# Supplement: per-capture models and technical-feature audit, outside the main biological figures.
rates<-read("doublet_scenarios.csv");rates$predicted_percent<-100*rates$predicted_doublets/rates$cells
capcols<-c(GSM9868022="#777777",GSM9868023="#D55E00",GSM9868024="#0072B2")
a<-ggplot(rates,aes(100*assumed_rate,predicted_percent,color=capture))+geom_line(linewidth=.45)+geom_point(size=1.6)+scale_color_manual(values=capcols,labels=c("Fibrosis capture","Vehicle capture","CCS1477 capture"),name=NULL)+labs(x="Assumed residual doublet rate (%)",y="Predicted doublets (%)",title="Capture-wise doublet-rate sensitivity")+theme(legend.position="bottom")
frac<-read("feature_class_fractions.csv");frac<-frac[frac$condition%in%names(condition_colors)&frac$population%in%poporder,]
frac$population<-factor(frac$population,levels=poporder)
frac$condition<-factor(frac$condition,levels=names(condition_colors))
make_fraction<-function(type,title){z<-frac[frac$feature_class==type,];ggplot(z,aes(population,100*umi_fraction,color=condition))+geom_point(position=position_dodge(width=.45),size=1.3,alpha=.85)+scale_color_manual(values=condition_colors,labels=condition_labels,name=NULL)+scale_x_discrete(labels=popnames)+labs(x=NULL,y="Raw pseudobulk UMIs (%)",title=title)+theme(axis.text.x=element_text(angle=35,hjust=1),legend.position="bottom")}
b<-make_fraction("mitochondrial","Mitochondrial UMI fraction")
c<-make_fraction("ribosomal_protein","Ribosomal-protein UMI fraction")
d<-make_fraction("Gm42418_unresolved","Gm42418 UMI fraction")
fS<-(a+b)/(c+d)+plot_annotation(tag_levels="A")
export(fS,"Supplementary_Figure_S1_technical_audit",180,143)
spec<-data.frame(figure=c("Figure_1_cell_atlas","Figure_2_mouse_level_expression","Figure_3_hallmark_programs","Supplementary_Figure_S1_technical_audit"),width_mm=180,height_mm=c(113,180,178,143),png_dpi=600,font=font)
if(file.exists(file.path(pub,"export_specifications.csv"))){oldspec<-read.csv(file.path(pub,"export_specifications.csv"));spec<-rbind(spec,oldspec[!oldspec$figure %in% spec$figure,])}
write.csv(spec,file.path(pub,"export_specifications.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(pub,"sessionInfo.txt"))
