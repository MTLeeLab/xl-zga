# read in tables
bl <- read.table('enh/quant/regions/both_on_L.rpm.bed', header=F, sep="\t")
bs <- read.table('enh/quant/regions/both_on_S.rpm.bed', header=F, sep="\t")
ll <- read.table('enh/quant/regions/L_on_L.rpm.bed', header=F, sep="\t")
ls <- read.table('enh/quant/regions/L_on_S.rpm.bed', header=F, sep="\t")
sl <- read.table('enh/quant/regions/S_on_L.rpm.bed', header=F, sep="\t")
ss <- read.table('enh/quant/regions/S_on_S.rpm.bed', header=F, sep="\t")
# name columns
names(bl) <- c("chr","start","stop","enhid","pairid","k27_raw","atac_raw","k27_rpm","atac_rpm")
names(bs) <- c("chr","start","stop","enhid","pairid","k27_raw","atac_raw","k27_rpm","atac_rpm")
names(ll) <- c("chr","start","stop","enhid","pairid","k27_raw","atac_raw","k27_rpm","atac_rpm")
names(ls) <- c("chr","start","stop","enhid","pairid","k27_raw","atac_raw","k27_rpm","atac_rpm")
names(sl) <- c("chr","start","stop","enhid","pairid","k27_raw","atac_raw","k27_rpm","atac_rpm")
names(ss) <- c("chr","start","stop","enhid","pairid","k27_raw","atac_raw","k27_rpm","atac_rpm")
# merge LS together
both <- merge(bl, bs, by="pairid", suffix=c(".L",".S"))
lon <- merge(ll, ls, by="pairid", suffix=c(".L",".S"))
son <- merge(sl, ss, by="pairid", suffix=c(".L",".S"))
# boxplot
pdf("enh/quant/k27atac_ls_quant.pdf",7,7)
par(pty="s")
par(mfrow=c(1,2))
boxplot(log2(lon$k27_rpm.L+0.04), log2(lon$k27_rpm.S+0.04), log2(both$k27_rpm.L+0.04), log2(both$k27_rpm.S+0.04), log2(son$k27_rpm.L+0.04), log2(son$k27_rpm.S+0.04), names=c("lonl","lons","bothl","boths","sonl","sons"), cex=0.5, main = "LvS k27ac coverage", ylab = "log2 rpm", las=2)
boxplot(log2(lon$atac_rpm.L+0.07), log2(lon$atac_rpm.S+0.07), log2(both$atac_rpm.L+0.07), log2(both$atac_rpm.S+0.07), log2(son$atac_rpm.L+0.07), log2(son$atac_rpm.S+0.07), names=c("lonl","lons","bothl","boths","sonl","sons"), cex=0.5, main = "LvS atac coverage", ylab = "log2 rpm", las=2)
dev.off()
# stat test
# k27 - L;LS;S
t.test(log2(lon$k27_rpm.L+0.04), log2(lon$k27_rpm.S+0.04), paired=T, mu=0)$p.value
t.test(log2(both$k27_rpm.L+0.04), log2(both$k27_rpm.S+0.04), paired=T, mu=0)$p.value
t.test(log2(son$k27_rpm.L+0.04), log2(son$k27_rpm.S+0.04), paired=T, mu=0)$p.value
# atac - L;LS;S
t.test(log2(lon$atac_rpm.L+0.07), log2(lon$atac_rpm.S+0.07), paired=T, mu=0)$p.value
t.test(log2(both$atac_rpm.L+0.07), log2(both$atac_rpm.S+0.07), paired=T, mu=0)$p.value
t.test(log2(son$atac_rpm.L+0.07), log2(son$atac_rpm.S+0.07), paired=T, mu=0)$p.value
