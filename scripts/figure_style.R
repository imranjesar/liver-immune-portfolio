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
