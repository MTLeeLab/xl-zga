# sliding window in R
# read in headless matrices
LpromK4 <- read.table('../prom/heatmaps/R_matrix/ls_fw_conf_chx_L_s9k4_lfc.txt',header=F,sep="\t")
SpromK4 <- read.table('../prom/heatmaps/R_matrix/ls_fw_conf_chx_S_s9k4_lfc.txt',header=F,sep="\t")
LpromK27 <- read.table('../prom/heatmaps/R_matrix/ls_fw_conf_chx_L_s8k27_lfc.txt',header=F,sep="\t")
SpromK27 <- read.table('../prom/heatmaps/R_matrix/ls_fw_conf_chx_S_s8k27_lfc.txt',header=F,sep="\t")
# name lfc column
names(LpromK4)[77]<-"rna_lfc"
names(SpromK4)[77]<-"rna_lfc"
names(LpromK27)[77]<-"rna_lfc"
names(SpromK27)[77]<-"rna_lfc"
# sum total signal per row
library(zoo)
LpromK4$count <- rowSums(LpromK4[,7:76])
SpromK4$count <- rowSums(SpromK4[,7:76])
LpromK27$count <- rowSums(LpromK27[,7:76])
SpromK27$count <- rowSums(SpromK27[,7:76])
# calculate sliding window
l4 <- rollapply(LpromK4$count, width=110, FUN=mean)
s4 <- rollapply(SpromK4$count, width=110, FUN=mean)
l4_2 <- as.data.frame(l4)
s4_2 <- as.data.frame(s4)
l4_2$bin <- 1:nrow(l4_2)
s4_2$bin <- 1:nrow(s4_2)
l27 <- rollapply(LpromK27$count, width=110, FUN=mean)
s27  <- rollapply(SpromK27$count, width=110, FUN=mean)
l27_2 <- as.data.frame(l27)
s27_2 <- as.data.frame(s27)
l27_2$bin <- 1:nrow(l27_2)
s27_2$bin <- 1:nrow(s27_2)
# plot meta plots
pdf(file="../prom/pdfs/k4_fw_meta.pdf",height=7,width=7)
plot(l4_2$bin, l4_2$l4, type='l', ylim = c(0,75), col = "red")
points(s4_2$bin, s4_2$s4, type='l', ylim = c(0,75), col = "blue")
dev.off()
pdf(file="../prom/pdfs/k27_fw_meta.pdf",height=7,width=7)
plot(l27_2$bin, l27_2$l27, type='l', ylim = c(0,20), col = "red")
points(s27_2$bin, s27_2$s27, type='l', ylim = c(0,20), col = "blue")
dev.off()
# stat test
# k4me3 count
cor.test(LpromK4$rna_lfc, log2(LpromK4$count+1)-log2(SpromK4$count+1), alternative="two.sided",method="pearson")
cor.test(LpromK4$rna_lfc, log2(LpromK4$count+1)-log2(SpromK4$count+1), alternative="two.sided",method="pearson")$p.value
# k27ac count
cor.test(LpromK27$rna_lfc, log2(LpromK27$count+1)-log2(SpromK27$count+1), alternative="two.sided",method="pearson")
cor.test(LpromK27$rna_lfc, log2(LpromK27$count+1)-log2(SpromK27$count+1), alternative="two.sided",method="pearson")$p.value
