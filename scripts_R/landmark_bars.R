# read tables
both <- read.table('enh/landmark/match/both_enh_match.txt', header=F, sep="\t")
lon <- read.table('enh/landmark/match/L_enh_match.txt', header=F, sep="\t")
son <- read.table('enh/landmark/match/S_enh_match.txt', header=F, sep="\t")
# name columns
names(both)<-c("pairid","L_enh","S_enh","US_match","DS_match")
names(lon)<-c("pairid","L_enh","S_enh","US_match","DS_match")
names(son)<-c("pairid","L_enh","S_enh","US_match","DS_match")
# calculate proportions
bars<-c(nrow(subset(both, US_match != 0 & DS_match != 0))/nrow(both), nrow(subset(lon, US_match != 0 & DS_match != 0))/nrow(lon), nrow(subset(son, US_match != 0 & DS_match != 0))/nrow(son))
# barplot
pdf("enh/landmark/match_barplot.pdf",7,7)
barplot(bars*100, beside=T, names=c("LS","L","S"), las=2, ylim=c(0,100))
dev.off()
