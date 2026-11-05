suppressPackageStartupMessages(library(ggplot2))
dir.create("results/sc/doublets/figures",showWarnings=FALSE)
d<-read.csv("results/sc/doublets/cell_doublet_summary.csv")
r<-read.csv("results/sc/doublets/scenario_summary.csv")
r$predicted_pct<-100*r$predicted_doublets/r$cells
p<-ggplot(r,aes(100*assumed_rate,predicted_pct,color=capture))+geom_line()+geom_point(size=2)+theme_bw()+labs(x="Assumed residual doublet rate (%)",y="Predicted doublets (%)",color="Capture library")
ggsave("results/sc/doublets/figures/rate_scenarios.png",p,width=8,height=5,dpi=150,bg="white")
x<-aggregate(consensus_doublet~baseline_cluster,d,mean)
x$baseline_cluster<-factor(x$baseline_cluster,levels=0:25)
p<-ggplot(x,aes(baseline_cluster,100*consensus_doublet))+geom_col(fill="#377eb8")+theme_bw()+labs(x="Baseline cluster (original numbering)",y="Consensus predicted doublets (%)")
ggsave("results/sc/doublets/figures/consensus_by_baseline_cluster.png",p,width=11,height=5,dpi=150,bg="white")
b<-read.csv("results/sc/qc_after.csv");out<-list()
for(tag in sort(unique(b$sample_id))){
 nbase<-sum(b$sample_id==tag);strict<-sum(d$sample_id==tag);calls<-sum(d$sample_id==tag & d$consensus_doublet)
 out[[tag]]<-data.frame(sample_id=tag,status=c("Retained after consensus filter","Consensus predicted doublet","Below 500 genes"),cells=c(strict-calls,calls,nbase-strict))
}
x<-do.call(rbind,out)
p<-ggplot(x,aes(sample_id,cells,fill=status))+geom_col()+theme_bw()+labs(x="Mouse tag",y="Cells from baseline QC",fill="Disposition")
ggsave("results/sc/doublets/figures/sample_retention.png",p,width=11,height=6,dpi=150,bg="white")
