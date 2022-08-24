# make beeswarm plots
# read in centered files
fw_conf_chx_L <- read.table('../prom/heatmaps/regions/ls_fw_conf_chx_Lprom.bed', header=F, sep='\t')
fw_conf_chx_S <- read.table('../prom/heatmaps/regions/ls_fw_conf_chx_Sprom.bed', header=F, sep='\t')
# read in upstream files
fw_conf_chxUS_L <- read.table('../prom/data/ls_fw_conf_chx_Lprom_k27US.bed', header=F, sep='\t')
fw_conf_chxUS_S <- read.table('../prom/data/ls_fw_conf_chx_Sprom_k27US.bed', header=F, sep='\t')
# name columns
names(fw_conf_chx_L)<-c("chr.L","start.L","stop.L","geneid.L","name","strand.L","act_cat","rna_fc_chx","s8_cov.L","s9_cov.L","k27_cov.L")
names(fw_conf_chx_S)<-c("chr.S","start.S","stop.S","geneid.S","name","strand.S","act_cat","rna_fc_chx","s8_cov.S","s9_cov.S","k27_cov.S")
names(fw_conf_chxUS_L)<-c("chr.L","start.L","stop.L","geneid.L","name","strand.L","act_cat","rna_fc_chx","k27US_cov.L")
names(fw_conf_chxUS_S)<-c("chr.S","start.S","stop.S","geneid.S","name","strand.S","act_cat","rna_fc_chx","k27US_cov.S")
# calculate rpkm
s8pm <- 1367050/1000000
s9pm <- 21944333/1000000
k27pm <- 23598517/1000000
fw_conf_chx_L$s8_rpkm.L <- (fw_conf_chx_L$s8_cov.L/s8pm)/((fw_conf_chx_L$stop.L-fw_conf_chx_L$start.L)/1000)
fw_conf_chx_S$s8_rpkm.S <- (fw_conf_chx_S$s8_cov.S/s8pm)/((fw_conf_chx_S$stop.S-fw_conf_chx_S$start.S)/1000)
fw_conf_chx_L$s9_rpkm.L <- (fw_conf_chx_L$s9_cov.L/s9pm)/((fw_conf_chx_L$stop.L-fw_conf_chx_L$start.L)/1000)
fw_conf_chx_S$s9_rpkm.S <- (fw_conf_chx_S$s9_cov.S/s9pm)/((fw_conf_chx_S$stop.S-fw_conf_chx_S$start.S)/1000)
fw_conf_chx_L$k27_rpkm.L <- (fw_conf_chx_L$k27_cov.L/k27pm)/((fw_conf_chx_L$stop.L-fw_conf_chx_L$start.L)/1000)
fw_conf_chx_S$k27_rpkm.S <- (fw_conf_chx_S$k27_cov.S/k27pm)/((fw_conf_chx_S$stop.S-fw_conf_chx_S$start.S)/1000)
fw_conf_chxUS_L$k27_rpkm.L <- (fw_conf_chxUS_L$k27US_cov.L/k27pm)/((fw_conf_chxUS_L$stop.L-fw_conf_chxUS_L$start.L)/1000)
fw_conf_chxUS_S$k27_rpkm.S <- (fw_conf_chxUS_S$k27US_cov.S/k27pm)/((fw_conf_chxUS_S$stop.S-fw_conf_chxUS_S$start.S)/1000)
# join S coverage to L homeolog
fw_conf_chx_L$s8_rpkm.S <- fw_conf_chx_S$s8_rpkm.S
fw_conf_chx_L$s9_rpkm.S <- fw_conf_chx_S$s9_rpkm.S
fw_conf_chx_L$k27_rpkm.S <- fw_conf_chx_S$k27_rpkm.S
fw_conf_chxUS_L$k27_rpkm.S <- fw_conf_chxUS_S$k27_rpkm.S
# calculate L:S ratio
fw_conf_chx_L$Ratio_s8 <- log2(fw_conf_chx_L$s8_rpkm.L+0.02)-log2(fw_conf_chx_L$s8_rpkm.S+0.02)
fw_conf_chx_L$Ratio_s9 <- log2(fw_conf_chx_L$s9_rpkm.L+0.02)-log2(fw_conf_chx_L$s9_rpkm.S+0.02)
fw_conf_chx_L$Ratio_k27 <- log2(fw_conf_chx_L$k27_rpkm.L+0.02)-log2(fw_conf_chx_L$k27_rpkm.S+0.02)
fw_conf_chxUS_L$Ratio_k27US <- log2(fw_conf_chxUS_L$k27_rpkm.L+0.02)-log2(fw_conf_chxUS_L$k27_rpkm.S+0.02)
# plot with beeswarm and boxplot
# s9 k4me3 with stats
library(beeswarm)
pdf(file='../prom/pdfs/s9_3bees_de.pdf',height=7,width=7)
beeswarm(Ratio_s9~act_cat, data=fw_conf_chx_L, main = "log2 L:S ratio of activated s9 H3K4me3 promoter coverage", ylab = "log2 L:S coverage", xlab = "Category", pch=21, cex = 0.25, col = c("#330066","#FF3030", "#1E90FF"), bg = c("#33006615","#FF303015","#1E90FF15"), ylim = c(-6.5,6.5))
boxplot(Ratio_s9~act_cat, data=fw_conf_chx_L, main = "log2 L:S ratio of activated s9 H3K4me3 promoter coverage", ylab = "log2 L:S coverage", xlab = "Category", outline=FALSE, col = "#FFFFFF00", lwd=1, add=TRUE)
dev.off()
L <- subset(fw_conf_chx_L, act_cat == "Lon")
B <- subset(fw_conf_chx_L, act_cat == "Both")
S <- subset(fw_conf_chx_L, act_cat == "Son")
# count rows for stats
nrow(L)
nrow(S)
nrow(B)
# stats
t.test(L$s9_rpkm.L, L$s9_rpkm.S, paired=T)
t.test(L$s9_rpkm.L, L$s9_rpkm.S, paired=T)$p.value
t.test(S$s9_rpkm.L, S$s9_rpkm.S, paired=T)
t.test(S$s9_rpkm.L, S$s9_rpkm.S, paired=T)$p.value
t.test(B$s9_rpkm.L, B$s9_rpkm.S, paired=T)
t.test(B$s9_rpkm.L, B$s9_rpkm.S, paired=T)$p.value
# s8 k4me3 with stats
pdf(file='../prom/pdfs/s8_3bees_de.pdf',height=7,width=7)
beeswarm(Ratio_s8~act_cat, data=fw_conf_chx_L, main = "log2 L:S ratio of activated s8 H3K4me3 promoter coverage", ylab = "log2 L:S coverage", xlab = "Category", pch=21, cex = 0.25, col = c("#330066", "#FF3030", "#1E90FF"), bg = c("#33006615","#FF303015","#1E90FF15"), ylim = c(-6.5,6.5))
boxplot(Ratio_s8~act_cat, data=fw_conf_chx_L, main = "log2 L:S ratio of activated s8 H3K4me3 promoter coverage", ylab = "log2 L:S coverage", xlab = "Category", outline=FALSE, col = "#FFFFFF00", lwd=1, add=TRUE)
dev.off()
L <- subset(fw_conf_chx_L, act_cat == "Lon")
B <- subset(fw_conf_chx_L, act_cat == "Both")
S <- subset(fw_conf_chx_L, act_cat == "Son")
t.test(L$s8_rpkm.L, L$s8_rpkm.S, paired=T)
t.test(L$s8_rpkm.L, L$s8_rpkm.S, paired=T)$p.value
t.test(S$s8_rpkm.L, S$s8_rpkm.S, paired=T)
t.test(S$s8_rpkm.L, S$s8_rpkm.S, paired=T)$p.value
t.test(B$s8_rpkm.L, B$s8_rpkm.S, paired=T)
t.test(B$s8_rpkm.L, B$s8_rpkm.S, paired=T)$p.value
# s8 k27ac with stats
pdf(file='../prom/pdfs/k27_3bees_de.pdf',height=7,width=7)
beeswarm(Ratio_k27~act_cat, data=fw_conf_chx_L, main = "log2 L:S ratio of activated s8 H3K27ac promoter coverage", ylab = "log2 L:S coverage", xlab = "Category", pch=21, cex = 0.25, col = c("#330066", "#FF3030", "#1E90FF"), bg = c("#33006615","#FF303015","#1E90FF15"), ylim = c(-6.5,6.5))
boxplot(Ratio_k27~act_cat, data=fw_conf_chx_L, main = "log2 L:S ratio of activated s8 H3K27ac promoter coverage", ylab = "log2 L:S coverage", xlab = "Category", outline=FALSE, col = "#FFFFFF00", lwd=1, add=TRUE)
dev.off()
L <- subset(fw_conf_chx_L, act_cat == "Lon")
B <- subset(fw_conf_chx_L, act_cat == "Both")
S <- subset(fw_conf_chx_L, act_cat == "Son")
# stat test
t.test(L$k27_rpkm.L, L$k27_rpkm.S, paired=T)
t.test(L$k27_rpkm.L, L$k27_rpkm.S, paired=T)$p.value
t.test(S$k27_rpkm.L, S$k27_rpkm.S, paired=T)
t.test(S$k27_rpkm.L, S$k27_rpkm.S, paired=T)$p.value
t.test(B$k27_rpkm.L, B$k27_rpkm.S, paired=T)
t.test(B$k27_rpkm.L, B$k27_rpkm.S, paired=T)$p.value
# US k27ac with stats
pdf(file='../prom/pdfs/k27US_3bees_de.pdf',height=7,width=7)
beeswarm(Ratio_k27US~act_cat, data=fw_conf_chxUS_L, main = "log2 L:S ratio of activated s8 H3K27ac upstream coverage", ylab = "log2 L:S coverage", xlab = "Category", pch=21, cex = 0.25, col = c("#330066", "#FF3030", "#1E90FF"), bg = c("#33006615","#FF303015","#1E90FF15"), ylim = c(-6.5,6.5))
boxplot(Ratio_k27US~act_cat, data=fw_conf_chxUS_L, main = "log2 L:S ratio of activated s8 H3K27ac upstream coverage", ylab = "log2 L:S coverage", xlab = "Category", outline=FALSE, col = "#FFFFFF00", lwd=1, add=TRUE)
dev.off()
L <- subset(fw_conf_chxUS_L, act_cat == "Lon")
B <- subset(fw_conf_chxUS_L, act_cat == "Both")
S <- subset(fw_conf_chxUS_L, act_cat == "Son")
# count rows for stats
nrow(L)
nrow(S)
nrow(B)
# stat test
t.test(L$k27_rpkm.L, L$k27_rpkm.S, paired=T)
t.test(L$k27_rpkm.L, L$k27_rpkm.S, paired=T)$p.value
t.test(S$k27_rpkm.L, S$k27_rpkm.S, paired=T)
t.test(S$k27_rpkm.L, S$k27_rpkm.S, paired=T)$p.value
t.test(B$k27_rpkm.L, B$k27_rpkm.S, paired=T)
t.test(B$k27_rpkm.L, B$k27_rpkm.S, paired=T)$p.value
# write out tables
write.table(fw_conf_chx_L, '../prom/data/tss_centered_bee_coverage.txt', sep="\t", quote=F, col.names=T, row.names=F)
write.table(fw_conf_chxUS_L, '../prom/data/k27_US_bee_coverage.txt', sep="\t", quote=F, col.names=T, row.names=F)
