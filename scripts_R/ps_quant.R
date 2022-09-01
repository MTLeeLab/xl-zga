# read it tables
bl <- read.table('tf/quant/regions/both_on_L.rpm.bed', header=F, sep="\t")
bs <- read.table('tf/quant/regions/both_on_S.rpm.bed', header=F, sep="\t")
ll <- read.table('tf/quant/regions/L_on_L.rpm.bed', header=F, sep="\t")
ls <- read.table('tf/quant/regions/L_on_S.rpm.bed', header=F, sep="\t")
sl <- read.table('tf/quant/regions/S_on_L.rpm.bed', header=F, sep="\t")
ss <- read.table('tf/quant/regions/S_on_S.rpm.bed', header=F, sep="\t")
# name columns
names(bl) <- c("chr","start","stop","enhid","pairid","pou_raw","sox_raw","pou_rpm","sox_rpm")
names(bs) <- c("chr","start","stop","enhid","pairid","pou_raw","sox_raw","pou_rpm","sox_rpm")
names(ll) <- c("chr","start","stop","enhid","pairid","pou_raw","sox_raw","pou_rpm","sox_rpm")
names(ls) <- c("chr","start","stop","enhid","pairid","pou_raw","sox_raw","pou_rpm","sox_rpm")
names(sl) <- c("chr","start","stop","enhid","pairid","pou_raw","sox_raw","pou_rpm","sox_rpm")
names(ss) <- c("chr","start","stop","enhid","pairid","pou_raw","sox_raw","pou_rpm","sox_rpm")
# merge LS together
both <- merge(bl, bs, by="pairid", suffix=c(".L",".S"))
lon <- merge(ll, ls, by="pairid", suffix=c(".L",".S"))
son <- merge(sl, ss, by="pairid", suffix=c(".L",".S"))
# boxplot
pdf("tf/quant/tf_ls_quant.pdf",7,7)
par(pty="s")
par(mfrow=c(1,2))
boxplot(log2(lon$pou_rpm.L+0.08), log2(lon$pou_rpm.S+0.08), log2(both$pou_rpm.L+0.08), log2(both$pou_rpm.S+0.08), log2(son$pou_rpm.L+0.08), log2(son$pou_rpm.S+0.08), names=c("lonl","lons","bothl","boths","sonl","sons"), cex=0.5, main = "LvS Pou5f3.3 coverage", ylab = "log2 rpm", las=2)
boxplot(log2(lon$sox_rpm.L+0.05), log2(lon$sox_rpm.S+0.05), log2(both$sox_rpm.L+0.05), log2(both$sox_rpm.S+0.05), log2(son$sox_rpm.L+0.05), log2(son$sox_rpm.S+0.05), names=c("lonl","lons","bothl","boths","sonl","sons"), cex=0.5, main = "LvS Sox3 coverage", ylab = "log2 rpm", las=2)
dev.off()
# t-test LS coverage for each factor - log2 with 1 read smoothed
# pou
t.test(log2(lon$pou_rpm.L+0.08), log2(lon$pou_rpm.S+0.08), paired=T, mu=0)$p.value
t.test(log2(both$pou_rpm.L+0.08), log2(both$pou_rpm.S+0.08), paired=T, mu=0)$p.value
t.test(log2(son$pou_rpm.L+0.08), log2(son$pou_rpm.S+0.08), paired=T, mu=0)$p.value
# sox
t.test(log2(lon$sox_rpm.L+0.05), log2(lon$sox_rpm.S+0.05), paired=T, mu=0)$p.value
t.test(log2(both$sox_rpm.L+0.05), log2(both$sox_rpm.S+0.05), paired=T, mu=0)$p.value
t.test(log2(son$sox_rpm.L+0.05), log2(son$sox_rpm.S+0.05), paired=T, mu=0)$p.value

