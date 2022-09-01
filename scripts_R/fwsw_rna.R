# read in tables
mo <- read.table('tf/fsm_heatmaps/regions/fw_ps_cats.bed',header=F,sep="\t")
# name columns
names(mo)<-c("chr","start","stop","geneid","tss_exp","strand","act_trip_d","firstwave","ps_affect","category")
# add balancer column for heatmap.2
mo$balance <- as.numeric(1)
# make heatmap tables
mo_aff_hm <- mo[,c("ps_affect","balance")]
# plot heatmap
library(gplots)
hmcols<-colorRampPalette(c("#330066","white"))(50)
pdf('tf/fsm_heatmaps/mo_aff_rnaFC_heatmap.pdf',7,7)
heatmap.2(as.matrix(-mo_aff_hm), Rowv=FALSE, Colv=FALSE, col = hmcols, scale="none", margins=c(6,10), trace="none", dendrogram="none", density.info="none", cexCol=1, breaks = seq(-9,0, length.out = 51))
dev.off()
