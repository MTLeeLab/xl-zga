# read in homeolog rna-seq table
homeo_lfc <- read.table('../files/anno/homeolog_lfc_all.txt', header=T, sep="\t")
g <- homeo_lfc
g[is.na(g)] <- 0


# Make table of homeolog pairs where one or both are first-wave activated
gsw <- subset(g, (firstwave_conf.L == "X" | firstwave_conf.S == "X") & homeo_status == "LS")
fw_conf_chx <- gsw[,c("family", "Geneid.L", "Geneid.S", "act_trip_d.L", "firstwave_conf.L", "act_trip_d.S", "firstwave_conf.S")]


# calculate L vs S log fold difference of DESeq log fold differences
# for exons and introns, then select the higher

max_fold_diff = function(e_lfc_l, e_lfc_s, i_lfc_l, i_lfc_s) {
	
	#Zero out decreases, since we're calculating activation
	e_lfc_l[e_lfc_l < 0] <- 0
	e_lfc_s[e_lfc_s < 0] <- 0
	i_lfc_l[i_lfc_l < 0] <- 0
	i_lfc_s[i_lfc_s < 0] <- 0

	# Pick exons or introns to work with based on magnitude of lfc
	e_lfc_max <- ifelse(e_lfc_l > e_lfc_s, e_lfc_l, e_lfc_s)	
	i_lfc_max <- ifelse(i_lfc_l > i_lfc_s, i_lfc_l, i_lfc_s)	

	lfc.L <- ifelse(e_lfc_max > i_lfc_max, e_lfc_l, i_lfc_l)
	lfc.S <- ifelse(e_lfc_max > i_lfc_max, e_lfc_s, i_lfc_s)
	lfc_feature <- ifelse(e_lfc_max > i_lfc_max, 'e', 'i')

	data.frame(lfc.L, lfc.S, lfc_L_vs_S = lfc.L - lfc.S, lfc_feature)
}


act_table <- max_fold_diff(gsw$s9.5_dmso_vs_trip_d.log2FoldChange.e.L, gsw$s9.5_dmso_vs_trip_d.log2FoldChange.e.S, gsw$s9.5_dmso_vs_trip_d.log2FoldChange.i.L, gsw$s9.5_dmso_vs_trip_d.log2FoldChange.i.S) 
names(act_table) <- paste(names(act_table), '_wt', sep='')

fw_table <- max_fold_diff(gsw$s9.5_chx_vs_trip_d.log2FoldChange.e.L, gsw$s9.5_chx_vs_trip_d.log2FoldChange.e.S, gsw$s9.5_chx_vs_trip_d.log2FoldChange.i.L, gsw$s9.5_chx_vs_trip_d.log2FoldChange.i.S) 
names(fw_table) <- paste(names(fw_table), '_chx', sep='')

fw_conf_chx <- cbind(fw_conf_chx, act_table, fw_table)
fw_conf_chx <- fw_conf_chx[order(-fw_conf_chx$lfc_L_vs_S_chx),]



# Add code indicating first wave activation category
# If both activated in first wave, Both
# If activated in first wave in L and not activated at all in S, or vice versa: Lon, Son
# Otherwise "unk" for gray area
# 

fw_conf_chx$fw_cat <- with(fw_conf_chx,
				ifelse(firstwave_conf.L == "X" & act_trip_d.S == "..", "Lon",
					ifelse(firstwave_conf.L == "X" & firstwave_conf.S == "X", "Both",
						ifelse(firstwave_conf.S == "X" & act_trip_d.L == "..", "Son",
							"unk"))))



###
# Biplot of L vs S lfc
###
# make french flag biplot for L and S activation
# both from wildtype and firstwave foldchange

# Thresholding the color scale
fc <- ifelse(fw_conf_chx$lfc_L_vs_S_chx < -5, -5, 
			ifelse(fw_conf_chx$lfc_L_vs_S_chx > 5, 5,
				fw_conf_chx$lfc_L_vs_S_chx))

fc2 <- floor((60/11)*(fc+6))
colorRampPalette(c("firebrick1","#CCCCCC","dodgerblue"))(5)


# output should be this, use these hex values below
# [1] "#FF3030" "#E57E7E" "#CCCCCC" "#75AEE5" "#1E90FF"


RPBcol2<-rev(c(colorRampPalette(c("#FF3030","#E57E7E"))(10),colorRampPalette(c("#E57E7E","#CCCCCC"))(20),colorRampPalette(c("#CCCCCC","#75AEE5"))(20),colorRampPalette(c("#75AEE5","#1E90FF"))(10)))

library(scales)

# plot
pdf(file='../rna/pdfs/rna_fc_LS_biplot.pdf',width=7,heigh=7)
par(las = 1, pty = 's', mfrow=c(1,2))
with(fw_conf_chx, plot(lfc.S_chx, lfc.L_chx, xlim=c(-1,15), ylim=c(-1,15), pch = 21, cex = 0.75, col = RPBcol2[fc2], bg = alpha(RPBcol2[fc2], 0.2), axes = FALSE, main = 'First wave activation'))
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)

with(fw_conf_chx, plot(lfc.S_wt, lfc.L_wt, xlim=c(-1,15), ylim=c(-1,15), pch = 21, cex = 0.75, col = RPBcol2[fc2], bg = alpha(RPBcol2[fc2], 0.2), axes = FALSE, main = 'All activation'))
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
dev.off()


###
# Parallel heatmaps of L and S fold change
###

library(gplots)
hmcols <- colorRampPalette(c("white","#330066"))(50)

#Filter out gene pairs where one gene is on a scaffold
genechrs <- read.table('../files/anno/xl_gene_chrs.txt')  #51065
genechrs <- subset(genechrs, !grepl('chrUn', V1))  #45130

fw_conf_chx_chr <- subset(fw_conf_chx, Geneid.L %in% genechrs$V2 & Geneid.S %in% genechrs$V2)  #1213 / 1359

#2-fold boundaries for the heatmap
num_l <- sum(fw_conf_chx_chr$lfc_L_vs_S_chx > 1)
num_ls <- sum(fw_conf_chx_chr$lfc_L_vs_S_chx > -1)

pdf(file='../rna/pdfs/rna_fc_LShm_chr.pdf',width=5,heigh=7)
heatmap.2(as.matrix(fw_conf_chx_chr[,c("lfc.L_chx","lfc.S_chx")]), Rowv=FALSE, Colv=FALSE, col = hmcols, scale="none", margins=c(6,10), trace="none", dendrogram="none", density.info="none", cexCol=1, breaks = seq(0,12, length.out = 51), rowsep=c(0,num_l, num_ls), sepcolor='red')
dev.off()




# write out table
write.table(fw_conf_chx_chr, '../rna/ls_fw_conf_chx_chr.txt', sep="\t", quote=F, col.names=T, row.names=F)
