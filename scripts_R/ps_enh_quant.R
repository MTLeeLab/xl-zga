# read in tables
bl <- read.table("quant/enh_count/LpromLS_pscount.bed", header=F, sep="\t")
bs <- read.table("quant/enh_count/SpromLS_pscount.bed", header=F, sep="\t")
ll <- read.table("quant/enh_count/LpromL_pscount.bed", header=F, sep="\t")
sl <- read.table("quant/enh_count/LpromS_pscount.bed", header=F, sep="\t")
ls <- read.table("quant/enh_count/SpromL_pscount.bed", header=F, sep="\t")
ss <- read.table("quant/enh_count/SpromS_pscount.bed", header=F, sep="\t")
# name columns
names(bl) <- c("chr","start","stop","geneid","name","strand","act","firstwave_conf","act_cat","ps_fc","PSaff","cons_count","diff_count")
names(bs) <- c("chr","start","stop","geneid","name","strand","act","firstwave_conf","act_cat","ps_fc","PSaff","cons_count","diff_count")
names(ll) <- c("chr","start","stop","geneid","name","strand","act","firstwave_conf","act_cat","ps_fc","PSaff","cons_count","diff_count")
names(ls) <- c("chr","start","stop","geneid","name","strand","act","firstwave_conf","act_cat","ps_fc","PSaff","cons_count","diff_count")
names(sl) <- c("chr","start","stop","geneid","name","strand","act","firstwave_conf","act_cat","ps_fc","PSaff","cons_count","diff_count")
names(ss) <- c("chr","start","stop","geneid","name","strand","act","firstwave_conf","act_cat","ps_fc","PSaff","cons_count","diff_count")
# merge LS together
both <- merge(bl,bs,by="name",suffix=c(".L",".S"))
lon <- merge(ll,ls,by="name",suffix=c(".L",".S"))
son <- merge(sl,ss,by="name",suffix=c(".L",".S"))
# boxplot
pdf("quant/enh_count/LS_condiff_count.pdf",7,7)
par(pty="s")
par(mfrow=c(1,2))
boxplot(lon$cons_count.L, lon$cons_count.S, both$cons_count.L, both$cons_count.S, son$cons_count.L, son$cons_count.S, names=c("lonl","lons","bothl","boths","sonl","sons"), cex=0.5, main = "LvS con binding counts", ylab = "enh counts", las=2)
boxplot(lon$diff_count.L, lon$diff_count.S, both$diff_count.L, both$diff_count.S, son$diff_count.L, son$diff_count.S, names=c("lonl","lons","bothl","boths","sonl","sons"), cex=0.5, main = "LvS diff binding counts", ylab = "enh counts", las=2)
dev.off()
# stat test
# paired t - cons;diff - L;LS;S
# cons
t.test(lon$cons_count.L, lon$cons_count.S, paired=T, mu=0)$p.value
t.test(both$cons_count.L, both$cons_count.S, paired=T, mu=0)$p.value
t.test(son$cons_count.L, son$cons_count.S, paired=T, mu=0)$p.value
# diff
t.test(lon$diff_count.L, lon$diff_count.S, paired=T, mu=0)$p.value
t.test(both$diff_count.L, both$diff_count.S, paired=T, mu=0)$p.value
t.test(son$diff_count.L, son$diff_count.S, paired=T, mu=0)$p.value
