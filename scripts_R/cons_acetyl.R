# read in table
xl_xt_zf <- read.table('hom/xl_xt_zf_enh_full.txt', header=F, sep="\t")
# name columns
names(xl_xt_zf)<-c("pairid","chr.L","start.L","stop.L","enhid.L","chr.S","start.S","stop.S","enhid.S","trop_coord","trop_acetyl","zf_coord","zf_acetyl")
# count total number of enhancers
nrow(subset(xl_xt_zf, grepl("LS-",xl_xt_zf$pairid)))
nrow(subset(xl_xt_zf, grepl("L-",xl_xt_zf$pairid)))
nrow(subset(xl_xt_zf, grepl("S-NA",xl_xt_zf$pairid)))
nrow(subset(xl_xt_zf, grepl("l-",xl_xt_zf$pairid)))
nrow(subset(xl_xt_zf, grepl("s-",xl_xt_zf$pairid)))
# count number not lifted to trop
nrow(subset(xl_xt_zf, grepl("LS-",xl_xt_zf$pairid) & is.na(xl_xt_zf$trop_coord)))
nrow(subset(xl_xt_zf, grepl("L-",xl_xt_zf$pairid) & is.na(xl_xt_zf$trop_coord)))
nrow(subset(xl_xt_zf, grepl("S-NA",xl_xt_zf$pairid) & is.na(xl_xt_zf$trop_coord)))
nrow(subset(xl_xt_zf, grepl("l-",xl_xt_zf$pairid) & is.na(xl_xt_zf$trop_coord)))
nrow(subset(xl_xt_zf, grepl("s-",xl_xt_zf$pairid) & is.na(xl_xt_zf$trop_coord)))
# count number lifted but not acetylated
nrow(subset(xl_xt_zf, grepl("LS-",xl_xt_zf$pairid) & xl_xt_zf$trop_acetyl == "."))
nrow(subset(xl_xt_zf, grepl("L-",xl_xt_zf$pairid) & xl_xt_zf$trop_acetyl == "."))
nrow(subset(xl_xt_zf, grepl("S-NA",xl_xt_zf$pairid) & xl_xt_zf$trop_acetyl == "."))
nrow(subset(xl_xt_zf, grepl("l-",xl_xt_zf$pairid) & xl_xt_zf$trop_acetyl == "."))
nrow(subset(xl_xt_zf, grepl("s-",xl_xt_zf$pairid) & xl_xt_zf$trop_acetyl == "."))
# count acetylated trop regions
nrow(subset(xl_xt_zf, grepl("LS-",xl_xt_zf$pairid) & xl_xt_zf$trop_acetyl == "X"))
nrow(subset(xl_xt_zf, grepl("L-",xl_xt_zf$pairid) & xl_xt_zf$trop_acetyl == "X"))
nrow(subset(xl_xt_zf, grepl("S-NA",xl_xt_zf$pairid) & xl_xt_zf$trop_acetyl == "X"))
nrow(subset(xl_xt_zf, grepl("l-",xl_xt_zf$pairid) & xl_xt_zf$trop_acetyl == "X"))
nrow(subset(xl_xt_zf, grepl("s-",xl_xt_zf$pairid) & xl_xt_zf$trop_acetyl == "X"))
# put counts into a data table for ploting
trop_bars  <- as.data.frame(matrix(c("LS","L","S","l","s",6586,5183,3394,2974,1096,763,909,628,1181,594,2058,2566,1661,1100,269,3765,1708,1105,693,233),5))
names(trop_bars)<-c("cat","total","non_lift","lift","acetyl")
trop_bars$total <- as.numeric(trop_bars$total)
trop_bars$non_lift <- as.numeric(trop_bars$non_lift)
trop_bars$lift <- as.numeric(trop_bars$lift)
trop_bars$acetyl <- as.numeric(trop_bars$acetyl)
trop_bars$non_acetyl <- trop_bars$non_lift+trop_bars$lift
#
# in zebrafish
# count number not lifted to zf
nrow(subset(xl_xt_zf, grepl("LS-",xl_xt_zf$pairid) & is.na(xl_xt_zf$zf_coord)))
nrow(subset(xl_xt_zf, grepl("L-",xl_xt_zf$pairid) & is.na(xl_xt_zf$zf_coord)))
nrow(subset(xl_xt_zf, grepl("S-NA",xl_xt_zf$pairid) & is.na(xl_xt_zf$zf_coord)))
nrow(subset(xl_xt_zf, grepl("l-",xl_xt_zf$pairid) & is.na(xl_xt_zf$zf_coord)))
nrow(subset(xl_xt_zf, grepl("s-",xl_xt_zf$pairid) & is.na(xl_xt_zf$zf_coord)))
# count number lifted but not acetylated
nrow(subset(xl_xt_zf, grepl("LS-",xl_xt_zf$pairid) & xl_xt_zf$zf_acetyl == "."))
nrow(subset(xl_xt_zf, grepl("L-",xl_xt_zf$pairid) & xl_xt_zf$zf_acetyl == "."))
nrow(subset(xl_xt_zf, grepl("S-NA",xl_xt_zf$pairid) & xl_xt_zf$zf_acetyl == "."))
nrow(subset(xl_xt_zf, grepl("l-",xl_xt_zf$pairid) & xl_xt_zf$zf_acetyl == "."))
nrow(subset(xl_xt_zf, grepl("s-",xl_xt_zf$pairid) & xl_xt_zf$zf_acetyl == "."))
# count number acetylated
nrow(subset(xl_xt_zf, grepl("LS-",xl_xt_zf$pairid) & xl_xt_zf$zf_acetyl == "X"))
nrow(subset(xl_xt_zf, grepl("L-",xl_xt_zf$pairid) & xl_xt_zf$zf_acetyl == "X"))
nrow(subset(xl_xt_zf, grepl("S-NA",xl_xt_zf$pairid) & xl_xt_zf$zf_acetyl == "X"))
nrow(subset(xl_xt_zf, grepl("l-",xl_xt_zf$pairid) & xl_xt_zf$zf_acetyl == "X"))
nrow(subset(xl_xt_zf, grepl("s-",xl_xt_zf$pairid) & xl_xt_zf$zf_acetyl == "X"))
# put counts into a data table to plot
zf_bars  <- as.data.frame(matrix(c("LS","L","S","l","s",6586,5183,3394,2974,1096,4679,4238,2715,2537,927,1042,684,502,281,105,865,261,177,156,64),5))
names(zf_bars)<-c("cat","total","non_lift","lift","acetyl")
zf_bars$total <- as.numeric(zf_bars$total)
zf_bars$non_lift <- as.numeric(zf_bars$non_lift)
zf_bars$lift <- as.numeric(zf_bars$lift)
zf_bars$acetyl <- as.numeric(zf_bars$acetyl)
# do distal enhancers only
# trop
nrow(subset(xl_xt_zf, grepl("LS-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid)))
nrow(subset(xl_xt_zf, grepl("L-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid)))
nrow(subset(xl_xt_zf, grepl("S-NA",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid)))
nrow(subset(xl_xt_zf, grepl("l-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid)))
nrow(subset(xl_xt_zf, grepl("s-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid)))
#
nrow(subset(xl_xt_zf, grepl("LS-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & is.na(xl_xt_zf$trop_coord)))
nrow(subset(xl_xt_zf, grepl("L-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & is.na(xl_xt_zf$trop_coord)))
nrow(subset(xl_xt_zf, grepl("S-NA",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & is.na(xl_xt_zf$trop_coord)))
nrow(subset(xl_xt_zf, grepl("l-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & is.na(xl_xt_zf$trop_coord)))
nrow(subset(xl_xt_zf, grepl("s-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & is.na(xl_xt_zf$trop_coord)))
#
nrow(subset(xl_xt_zf, grepl("LS-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & xl_xt_zf$trop_acetyl == "."))
nrow(subset(xl_xt_zf, grepl("L-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & xl_xt_zf$trop_acetyl == "."))
nrow(subset(xl_xt_zf, grepl("S-NA",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & xl_xt_zf$trop_acetyl == "."))
nrow(subset(xl_xt_zf, grepl("l-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & xl_xt_zf$trop_acetyl == "."))
nrow(subset(xl_xt_zf, grepl("s-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & xl_xt_zf$trop_acetyl == "."))
#
nrow(subset(xl_xt_zf, grepl("LS-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & xl_xt_zf$trop_acetyl == "X"))
nrow(subset(xl_xt_zf, grepl("L-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & xl_xt_zf$trop_acetyl == "X"))
nrow(subset(xl_xt_zf, grepl("S-NA",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & xl_xt_zf$trop_acetyl == "X"))
nrow(subset(xl_xt_zf, grepl("l-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & xl_xt_zf$trop_acetyl == "X"))
nrow(subset(xl_xt_zf, grepl("s-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & xl_xt_zf$trop_acetyl == "X"))
#
#
distal_trop_bars  <- as.data.frame(matrix(c("LS","L","S","l","s",1312,2561,1776,999,452,151,526,347,518,310,430,1362,977,324,85,731,673,452,157,57),5))
names(distal_trop_bars)<-c("cat","total","non_lift","lift","acetyl")
distal_trop_bars$total <- as.numeric(distal_trop_bars$total)
distal_trop_bars$non_lift <- as.numeric(distal_trop_bars$non_lift)
distal_trop_bars$lift <- as.numeric(distal_trop_bars$lift)
distal_trop_bars$acetyl <- as.numeric(distal_trop_bars$acetyl)
# zebrafish
nrow(subset(xl_xt_zf, grepl("LS-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & is.na(xl_xt_zf$zf_coord)))
nrow(subset(xl_xt_zf, grepl("L-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & is.na(xl_xt_zf$zf_coord)))
nrow(subset(xl_xt_zf, grepl("S-NA",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & is.na(xl_xt_zf$zf_coord)))
nrow(subset(xl_xt_zf, grepl("l-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & is.na(xl_xt_zf$zf_coord)))
nrow(subset(xl_xt_zf, grepl("s-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & is.na(xl_xt_zf$zf_coord)))
#
nrow(subset(xl_xt_zf, grepl("LS-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & xl_xt_zf$zf_acetyl == "."))
nrow(subset(xl_xt_zf, grepl("L-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & xl_xt_zf$zf_acetyl == "."))
nrow(subset(xl_xt_zf, grepl("S-NA",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & xl_xt_zf$zf_acetyl == "."))
nrow(subset(xl_xt_zf, grepl("l-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & xl_xt_zf$zf_acetyl == "."))
nrow(subset(xl_xt_zf, grepl("s-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & xl_xt_zf$zf_acetyl == "."))
#
nrow(subset(xl_xt_zf, grepl("LS-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & xl_xt_zf$zf_acetyl == "X"))
nrow(subset(xl_xt_zf, grepl("L-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & xl_xt_zf$zf_acetyl == "X"))
nrow(subset(xl_xt_zf, grepl("S-NA",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & xl_xt_zf$zf_acetyl == "X"))
nrow(subset(xl_xt_zf, grepl("l-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & xl_xt_zf$zf_acetyl == "X"))
nrow(subset(xl_xt_zf, grepl("s-",xl_xt_zf$pairid) & grepl("distal",xl_xt_zf$pairid) & xl_xt_zf$zf_acetyl == "X"))
#
distal_zf_bars  <- as.data.frame(matrix(c("LS","L","S","l","s",1312,2561,1776,999,452,960,2176,1482,885,416,276,354,271,111,32,76,31,23,3,4),5))
names(distal_zf_bars)<-c("cat","total","non_lift","lift","acetyl")
distal_zf_bars$total <- as.numeric(distal_zf_bars$total)
distal_zf_bars$non_lift <- as.numeric(distal_zf_bars$non_lift)
distal_zf_bars$lift <- as.numeric(distal_zf_bars$lift)
distal_zf_bars$acetyl <- as.numeric(distal_zf_bars$acetyl)
# plot proportions in stacked barplot
pdf("hom/enh_cons_barplots.pdf",7,7)
par(pty="s")
par(mfrow=c(2,2))
barplot(as.matrix(t((trop_bars[,c(3:5)]/trop_bars$total)*100)), beside=F, names.arg=trop_bars$cat, main="trop regulatory conservation", col=c("#CC6600","#00CC00","#6666FF"))
barplot(as.matrix(t((zf_bars[,c(3:5)]/zf_bars$total)*100)), beside=F, names.arg=trop_bars$cat, main="ZF regulatory conservation", col=c("#CC6600","#00CC00","#6666FF"))
barplot(as.matrix(t((distal_trop_bars[,c(3:5)]/distal_trop_bars$total)*100)), beside=F, names.arg=trop_bars$cat, main="trop enh conservation", col=c("#CC6600","#00CC00","#6666FF"))
barplot(as.matrix(t((distal_zf_bars[,c(3:5)]/distal_zf_bars$total)*100)), beside=F, names.arg=trop_bars$cat, main="ZF enh conservation", col=c("#CC6600","#00CC00","#6666FF"))
dev.off()
# make tables and export
reg_counts <- merge(trop_bars, zf_bars[,c(1,3:5)], by = "cat", suffix=c(".tro",".zf"))
enh_counts <- merge(distal_trop_bars, distal_zf_bars[,c(1,3:5)], by = "cat", suffix=c(".tro",".zf"))
write.table(reg_counts,"hom/reg_counts.txt",sep="\t",quote=F,col.names=T,row.names=F)
write.table(enh_counts,"hom/enh_counts.txt",sep="\t",quote=F,col.names=T,row.names=F)
