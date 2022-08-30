# load in counts table
enh_counts <- read.table('enh/distance/count_enh/total_enh_counts.txt', header=F, sep="\t")
# name columns
names(enh_counts)<-c("name","chr.L","start.L","stop.L","geneid.L","strand.L","chr.S","start.S","stop.S","geneid.S","strand.S","act_cat","LS_rna_chx_FC","both_L_length","Lon_L_length","both_S_length","Son_S_length")
# slideing window of total enh length
library(zoo)
wide=144
bl <- rollapply(enh_counts$both_L_length, width=wide, FUN=sum)
bs <- rollapply(enh_counts$both_S_length, width=wide, FUN=sum)
l <- rollapply(enh_counts$Lon_L_length, width=wide, FUN=sum)
s <- rollapply(enh_counts$Son_S_length, width=wide, FUN=sum)
bl2 <- as.data.frame(bl)
bs2 <- as.data.frame(bs)
l2 <- as.data.frame(l)
s2 <- as.data.frame(s)
bl2$norm <- (bl2$bl+1)/(wide*50000)
bs2$norm <- (bs2$bs+1)/(wide*50000)
l2$norm <- (l2$l+1)/(wide*50000)
s2$norm <- (s2$s+1)/(wide*50000)
bl2$bin <- 1:nrow(bl2)
bs2$bin <- 1:nrow(bs2)
l2$bin <- 1:nrow(l2)
s2$bin <- 1:nrow(s2)
# plot meta lines
pdf('enh/distance/count_enh/con_enh_lin.pdf',7,7)
plot(bs2$bin, bs2$norm*10000, type='l', ylim=c(0,250), col = "blue", main = "con enhancer density")
points(bl2$bin, bl2$norm*10000, type='l', ylim=c(0,250), col = "red")
dev.off()
pdf('enh/distance/count_enh/diff_enh_lin.pdf',7,7)
plot(l2$bin, l2$norm*10000, type='l', ylim=c(0,250), col = "maroon", main = "diff enhancer density")
points(s2$bin, s2$norm*10000, type='l', ylim=c(0,250), col = "navy")
dev.off()
# stat test
# differential enh density
cor.test(enh_counts$LS_rna_chx_FC, log2(enh_counts$Lon_L_length+250)-log2(enh_counts$Son_S_length+250), alternative="two.sided",method="pearson")
cor.test(enh_counts$LS_rna_chx_FC, log2(enh_counts$Lon_L_length+250)-log2(enh_counts$Son_S_length+250), alternative="two.sided",method="pearson")$p.value
# conserved enh density
cor.test(enh_counts$LS_rna_chx_FC, log2(enh_counts$both_L_length+250)-log2(enh_counts$both_S_length+250), alternative="two.sided",method="pearson")
cor.test(enh_counts$LS_rna_chx_FC, log2(enh_counts$both_L_length+250)-log2(enh_counts$both_S_length+250), alternative="two.sided",method="pearson")$p.value
