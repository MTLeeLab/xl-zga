# read in tables
L_tfs <- read.table('enh/tfs/L_tf_exp.txt',header=F,sep='\t')
S_tfs <- read.table('enh/tfs/S_tf_exp.txt',header=F,sep='\t')
# add subgenome column
L_tfs$subg <- "L"
S_tfs$subg <- "S"
# annotate maternal genes
tfs_total <- rbind(L_tfs, S_tfs)
tfs_total$mat <- ifelse(tfs_total$V3 >= 1, "X",".")
# name columns
names(tfs_total)[1]<-"Geneid"
names(tfs_total)[2]<-"Name"
names(tfs_total)[3]<-"st5_avg"
names(tfs_total)[4]<-"GO_terms"
# make beeswarms
library(beeswarm)
pdf('enh/tfs/LS_bee.pdf',7,7)
beeswarm(log2(st5_avg)~subg, data=subset(tfs_total, mat == "X"), ylim = c(0,10), spacing=2, method="swarm", main = "log2 maternal TF expression", ylab = "log2 TPM", xlab = "subgenome", pch =21, cex = 0.5, col = c("#FF3030", "#1E90FF"), bg = c("#FF303015","#1E90FF15"))
dev.off()
write.table(tfs_total, 'enh/tfs/tfs_go.txt', sep="\t", quote=F, row.names=F, col.names=T)
