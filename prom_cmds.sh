# commands for generateing k4me3/k27ac heatmaps and beeswarms over promoters
# plot CR heatmaps over TSS
# place s8&s9 k4me3, s8 k27ac, and s8&s9 no antibody bigwigs into files/bigwigs
# place k4me3 and k27ac beds into files/beds
# use chr size and selected tss file from tss_cmds.sh
# use ls_fw_conf_chx_chr.txt from rna_lfc.sh

TAB=$'\t'

K4_S8_BW=files/bigwigs/s8_k4me3_hq.bw
K4_S9_BW=files/bigwigs/s9_k4me3_hq.bw
S9_NA=files/bigwigs/s9_noab_hq.bw
K27_S8_BW=files/bigwigs/s8_k27ac_hq.bw
S8_NA=files/bigwigs/s8_noab_hq.bw
SIZE=files/anno/chr_sizes.txt
K4_S9_COV=files/beds/s9_k4me3_hq.bed
K4_S8_COV=files/beds/s8_k4me3_hq.bed
K27_S8_COV=files/beds/s8_k27ac_hq.bed
SELECT=tss/TSS_select.bed

:<<TSSCOV

mkdir -p prom/heatmaps
mkdir -p prom/heatmaps/matrix
mkdir -p prom/heatmaps/pdfs
mkdir -p prom/heatmaps/regions
mkdir -p prom/heatmaps/R_matrix

# plot k4me3 and k27ac over all TSS
# calculate coverage for sorting order
cut -f1-6 $SELECT | bedtools slop -b 100 -i stdin -g $SIZE | sort -k1,1 -k2,2n | bedtools coverage -sorted -counts -a stdin -b $K4_S8_COV | bedtools coverage -sorted -counts -a stdin -b $K4_S9_COV | bedtools coverage -sorted -counts -a stdin -b $K27_S8_COV | awk -v OFS="\t" -v FS="\t" '{k4pm = 1367050/1000000; k4pm2 = 21944333/1000000; k27pm = 23598517/1000000} {print $0,($7/k4pm)/1,($8/k4pm2)/1,($9/k27pm)/1}' | awk -v OFS="\t" -v FS="\t" '{print $0,($10+$11+$12)}' | sort -k13,13nr > tss/TSS_select_crSort.bed

# computematrix over sorted select tss - sort by k27ac
WD=prom/heatmaps

computeMatrix reference-point --referencePoint center -a 2000 -b 2000 -R $SELECT -S $K27_S8_BW --binSize 50 --sortRegions descend --outFileSortedRegions $WD/regions/prom_total_s8.k27sort.bed --missingDataAsZero -o $WD/matrix/prom_total_s8k27.s8sort.matrix -p8
computeMatrix reference-point --referencePoint center -a 2000 -b 2000 -R $WD/regions/prom_total_s8.k27sort.bed -S $K4_S9_BW --binSize 50 --sortRegions keep --missingDataAsZero -o $WD/matrix/prom_total_s9k4.s8sort.matrix -p8
computeMatrix reference-point --referencePoint center -a 2000 -b 2000 -R $WD/regions/prom_total_s8.k27sort.bed -S $S9_NA --binSize 50 --sortRegions keep --missingDataAsZero -o $WD/matrix/prom_total_s9NA.s8sort.matrix -p8
computeMatrix reference-point --referencePoint center -a 2000 -b 2000 -R $WD/regions/prom_total_s8.k27sort.bed -S $K4_S8_BW --binSize 50 --sortRegions keep --missingDataAsZero -o $WD/matrix/prom_total_s8k4.s8sort.matrix -p8
computeMatrix reference-point --referencePoint center -a 2000 -b 2000 -R $WD/regions/prom_total_s8.k27sort.bed -S $S8_NA --binSize 50 --sortRegions keep --missingDataAsZero -o $WD/matrix/prom_total_s8NA.s8sort.matrix -p8

# plot heatmaps over select TSS - keep sort order from previous section
plotHeatmap -m $WD/matrix/prom_total_s8k4.s8sort.matrix -o $WD/pdfs/prom_total_s8k4.s8sort.pdf --sortRegions keep --dpi 1000 --heatmapWidth 2 --heatmapHeight 7 --xAxisLabel '' --regionsLabel 'all promoters' --samplesLabel 'stage 8 HQ' --colorList white,sandybrown,saddlebrown
plotHeatmap -m $WD/matrix/prom_total_s9k4.s8sort.matrix -o $WD/pdfs/prom_total_s9k4.s8sort.pdf --sortRegions keep --dpi 1000 --heatmapWidth 2 --heatmapHeight 7 --xAxisLabel '' --regionsLabel 'all promoters' --samplesLabel 'stage 9 HQ' --colorList white,sandybrown,saddlebrown --zMin 0.0 --zMax 1.2
plotHeatmap -m $WD/matrix/prom_total_s9NA.s8sort.matrix -o $WD/pdfs/prom_total_s9NA.s8sort.pdf --sortRegions keep --dpi 1000 --heatmapWidth 2 --heatmapHeight 7 --xAxisLabel '' --regionsLabel 'all promoters' --samplesLabel 'stage 9 NA' --colorList white,sandybrown,saddlebrown --zMin 0.0 --zMax 1.2
plotHeatmap -m $WD/matrix/prom_total_s8k27.s8sort.matrix -o $WD/pdfs/prom_total_s8k27.s8sort.pdf --sortRegions keep --dpi 1000 --heatmapWidth 2 --heatmapHeight 7 --xAxisLabel '' --regionsLabel 'all promoters' --samplesLabel 'stage 8 HQ' --colorList white,darkmagenta --zMin 0.0 --zMax 0.6
plotHeatmap -m $WD/matrix/prom_total_s8NA.s8sort.matrix -o $WD/pdfs/prom_total_s8NA.s8sort.pdf --sortRegions keep --dpi 1000 --heatmapWidth 2 --heatmapHeight 7 --xAxisLabel '' --regionsLabel 'all promoters' --samplesLabel 'stage 8 NA' --colorList white,darkmagenta --zMin 0.0 --zMax 0.6

TSSCOV


:<<FWCOV
# plot coverage over firstwave activated genes
# add TSS to activated genes

# Prepare BED files for deeptools and downstream analyses
WD2=prom/data
mkdir -p $WD2

cut -f1-6 $SELECT | bedtools slop -b 500 -i stdin -g $SIZE | sort -k1,1 -k2,2n | bedtools coverage -counts -sorted -a stdin -b $K4_S8_COV | bedtools coverage -counts -sorted -a stdin -b $K4_S9_COV | bedtools coverage -counts -sorted -a stdin -b $K27_S8_COV | cut -f1-4,6-9 > $WD2/prom500_histCov.bed

# Join coverage based on L id, then by S id; sort by lfc L vs S (column 14)
LSFILE=rna/ls_fw_conf_chx_ps.txt_tmp
join -1 2 -2 4 -t $'\t' <(sed '1d' $LSFILE | sort -k2,2 -t "${TAB}") <(sort -k4,4 -t "${TAB}" $WD2/prom500_histCov.bed) | join -1 3 -2 4 -t $'\t' <(sort -k3,3 -t "${TAB}" -) <(sort -k4,4 -t "${TAB}" $WD2/prom500_histCov.bed) | sort -t "${TAB}" -k14,14nr > $WD2/ls_fw_conf_chx_PROM.txt
## columns are in a strange order -- column 1 is the S gene, column 2 is the L gene

# make +/-500 tss BED files for two sets of promoters (L and S) - including scaffolds
# chr.L, start, end, genename, geneID, strand, fw_cat, lfc L vs S, (coverage L)
awk -v FS="\t" -v OFS="\t" '{print $20,$21,$22,$2,$3,$23,$16,$17,$24,$25,$26}' $WD2/ls_fw_conf_chx_PROM.txt > $WD/regions/ls_fw_conf_chx_Lprom.bed

# chr.S, start, end, genename, geneID, strand, fw_cat, lfc L vs S, (coverage S)
awk -v FS="\t" -v OFS="\t" '{print $27,$28,$29,$1,$3,$30,$16,$17,$31,$32,$33}' $WD2/ls_fw_conf_chx_PROM.txt > $WD/regions/ls_fw_conf_chx_Sprom.bed


# make +/-500 tss BED for enhancer analysis. Code excluded gene pairs where one is on a scaffold, but that is redundant if the original LSFILE has already been filtered
# chr, start, end, genename, geneID, strand, fw_cat, lfc L vs S
grep -v chrUn $WD2/ls_fw_conf_chx_PROM.txt | awk -v FS="\t" -v OFS="\t" '{print $20,$21,$22,$2,$3,$23,$16,$17,$24,$25,$26}' > $WD/regions/ls_fw_conf_chx_Lprom_noscaff.bed
grep -v chrUn $WD2/ls_fw_conf_chx_PROM.txt | awk -v FS="\t" -v OFS="\t" '{print $27,$28,$29,$1,$3,$30,$16,$17,$31,$32,$33}' > $WD/regions/ls_fw_conf_chx_Sprom_noscaff.bed

# computeMatrix over firstwave genes in L:S lfc order - leave out scaffold genes to facilitate homeolog comparisons
computeMatrix reference-point --referencePoint center -a 3000 -b 500 -R $WD/regions/ls_fw_conf_chx_Lprom_noscaff.bed $WD/regions/ls_fw_conf_chx_Sprom_noscaff.bed -S $K4_S8_BW --binSize 50 --sortRegions keep --missingDataAsZero -o $WD/matrix/ls_fw_conf_chx_k4s8.matrix -p6
computeMatrix reference-point --referencePoint center -a 3000 -b 500 -R $WD/regions/ls_fw_conf_chx_Lprom_noscaff.bed $WD/regions/ls_fw_conf_chx_Sprom_noscaff.bed -S $K4_S9_BW --binSize 50 --sortRegions keep --missingDataAsZero -o $WD/matrix/ls_fw_conf_chx_k4s9.matrix -p6
computeMatrix reference-point --referencePoint center -a 500 -b 3000 -R $WD/regions/ls_fw_conf_chx_Lprom_noscaff.bed $WD/regions/ls_fw_conf_chx_Sprom_noscaff.bed -S $K27_S8_BW --binSize 50 --sortRegions keep --missingDataAsZero -o $WD/matrix/ls_fw_conf_chx_k27s8.matrix -p6

# plotHeatmaps
plotHeatmap -m $WD/matrix/ls_fw_conf_chx_k4s9.matrix -o $WD/pdfs/ls_fw_conf_chx_s9k4.pdf --sortRegions keep --dpi 1000 --heatmapWidth 2 --heatmapHeight 7 --xAxisLabel '' --regionsLabel 'L proms' 'S proms' --samplesLabel 'stage 9 HQ' --colorList white,sandybrown,saddlebrown --zMax 1.2 --zMin 0
plotHeatmap -m $WD/matrix/ls_fw_conf_chx_k4s8.matrix -o $WD/pdfs/ls_fw_conf_chx_s8k4.pdf --sortRegions keep --dpi 1000 --heatmapWidth 2 --heatmapHeight 7 --xAxisLabel '' --regionsLabel 'L proms' 'S proms' --samplesLabel 'stage 8 HQ' --colorList white,sandybrown,saddlebrown
plotHeatmap -m $WD/matrix/ls_fw_conf_chx_k27s8.matrix -o $WD/pdfs/ls_fw_conf_chx_s8k27.pdf --sortRegions keep --dpi 1000 --heatmapWidth 2 --heatmapHeight 7 --xAxisLabel '' --regionsLabel 'L proms' 'S proms' --samplesLabel 'stage 8 HQ' --colorList white,darkmagenta --zMax 0.6 --zMin 0

# make meta plots for s9k4 and s8k27 over firstwave selected TSS
# k4 goes into gene body; k27 goes upstream
# make output directory

mkdir -p prom/pdfs

# compute matrix over L and S individually
computeMatrix reference-point --referencePoint center -a 3000 -b 500 -R $WD/regions/ls_fw_conf_chx_Lprom_noscaff.bed -S $K4_S9_BW --binSize 50 --sortRegions keep --missingDataAsZero -o $WD/R_matrix/ls_fw_conf_chx_L_s9k4.matrix -p6
computeMatrix reference-point --referencePoint center -a 3000 -b 500 -R $WD/regions/ls_fw_conf_chx_Sprom_noscaff.bed -S $K4_S9_BW --binSize 50 --sortRegions keep --missingDataAsZero -o $WD/R_matrix/ls_fw_conf_chx_S_s9k4.matrix -p6
computeMatrix reference-point --referencePoint center -a 500 -b 3000 -R $WD/regions/ls_fw_conf_chx_Lprom_noscaff.bed -S $K27_S8_BW --binSize 50 --sortRegions keep --missingDataAsZero -o $WD/R_matrix/ls_fw_conf_chx_L_s8k27.matrix -p6
computeMatrix reference-point --referencePoint center -a 500 -b 3000 -R $WD/regions/ls_fw_conf_chx_Sprom_noscaff.bed -S $K27_S8_BW --binSize 50 --sortRegions keep --missingDataAsZero -o $WD/R_matrix/ls_fw_conf_chx_S_s8k27.matrix -p6

# unzip and remove the headerline from the matrix output
gunzip -c $WD/R_matrix/ls_fw_conf_chx_L_s9k4.matrix | sed '1d' > $WD/R_matrix/ls_fw_conf_chx_L_s9k4_headless.txt
gunzip -c $WD/R_matrix/ls_fw_conf_chx_S_s9k4.matrix | sed '1d' > $WD/R_matrix/ls_fw_conf_chx_S_s9k4_headless.txt
gunzip -c $WD/R_matrix/ls_fw_conf_chx_L_s8k27.matrix | sed '1d' > $WD/R_matrix/ls_fw_conf_chx_L_s8k27_headless.txt
gunzip -c $WD/R_matrix/ls_fw_conf_chx_S_s8k27.matrix | sed '1d' > $WD/R_matrix/ls_fw_conf_chx_S_s8k27_headless.txt

# add rna lfc to headless matrix - make sure both files are sorted by L:S rna lfc!
# column 11 for total, column 8 for noscaff
paste -d "${TAB}" $WD/R_matrix/ls_fw_conf_chx_L_s9k4_headless.txt <(cut -f8 $WD/regions/ls_fw_conf_chx_Lprom_noscaff.bed) > $WD/R_matrix/ls_fw_conf_chx_L_s9k4_lfc.txt
paste -d "${TAB}" $WD/R_matrix/ls_fw_conf_chx_S_s9k4_headless.txt <(cut -f8 $WD/regions/ls_fw_conf_chx_Sprom_noscaff.bed) > $WD/R_matrix/ls_fw_conf_chx_S_s9k4_lfc.txt
paste -d "${TAB}" $WD/R_matrix/ls_fw_conf_chx_L_s8k27_headless.txt <(cut -f8 $WD/regions/ls_fw_conf_chx_Lprom_noscaff.bed) > $WD/R_matrix/ls_fw_conf_chx_L_s8k27_lfc.txt
paste -d "${TAB}" $WD/R_matrix/ls_fw_conf_chx_S_s8k27_headless.txt <(cut -f8 $WD/regions/ls_fw_conf_chx_Sprom_noscaff.bed) > $WD/R_matrix/ls_fw_conf_chx_S_s8k27_lfc.txt


# call R script to output metaplots
Rscript scripts_R/LS_cr_metaplot_cmds.R
# P-values for cor.test on L vs S s9 K4me3 and s8 K27ac
#[1] 4.395795e-67
#[1] 2.562815e-26

# Old p values that included the scaffolds
#[1] 2.532302e-56
#[1] 6.932618e-23


FWCOV

:<<BEESWARM
WD=prom/heatmaps
WD2=prom/data

# make beeswarm of k4 and k27 coverage around TSS
# k4 regions generated previously; make upstream regions for k27ac
# flank L and S promoters upt to 3kb upstream

cut -f1-8 $WD/regions/ls_fw_conf_chx_Lprom_noscaff.bed | bedtools flank -l 2500 -r 0 -s -i - -g $SIZE | sort -k1,1 -k2,2n | bedtools coverage -sorted -counts -a stdin -b $K27_S8_COV | sort -k8,8nr > $WD2/ls_fw_conf_chx_Lprom_k27US.bed
cut -f1-8 $WD/regions/ls_fw_conf_chx_Sprom_noscaff.bed | bedtools flank -l 2500 -r 0 -s -i - -g $SIZE | sort -k1,1 -k2,2n | bedtools coverage -sorted -counts -a stdin -b $K27_S8_COV | sort -k8,8nr > $WD2/ls_fw_conf_chx_Sprom_k27US.bed


# call Rscript for beeswarms
Rscript scripts_R/beeswarm_rpkm.R

#P-values for paired t-tests on Lon, Son, Both on

# K4me3 s9
#[1] 1.153935e-25
#[1] 1.94822e-10
#[1] 0.6844665

# K4me3 s8
#[1] 5.822327e-09
#[1] 0.01291845
#[1] 0.3914565

# K27ac s8 centered tss
#[1] 1.546704e-25
#[1] 3.247503e-08
#[1] 0.0003751275

# K27ac s8 upstream
#[1] 2.38626e-16
#[1] 0.001273631
#[1] 0.003133992


### Original p values including scaffolds
# K4me3 s9
#[1] 1.414862e-29
#[1] 2.818353e-12
#[1] 0.8764236

# K4me3 s8
#[1] 1.387972e-10
#[1] 0.003713669
#[1] 0.1601666

# K27ac s8 centered tss
#[1] 1.14734e-29
#[1] 4.948292e-09
#[1] 0.001887095

# K27ac s8 upstream
#[1] 9.233843e-18
#[1] 0.0004424943
#[1] 0.003309835

BEESWARM

:<<MOTIF

# motif calling
WD3=prom/homer
mkdir -p $WD3
mkdir -p $WD3/fasta
mkdir -p $WD3/regions
mkdir -p $WD3/Lact
mkdir -p $WD3/Sact
# pull out 500bp around TSS
# chr,start,stop,geneid,gene name,strand
awk -v OFS="\t" -v FS="\t" 'midpt=int(($21+$22)/2) {if($16=="Lon")print $20,midpt,midpt+1,$2,$3,$23}' $WD2/ls_fw_conf_chx_PROM.txt | bedtools slop -b 250 -i stdin -g $SIZE > $WD3/regions/Lact_L_500bp.bed
awk -v OFS="\t" -v FS="\t" 'midpt=int(($21+$22)/2) {if($16=="Son")print $20,midpt,midpt+1,$2,$3,$23}' $WD2/ls_fw_conf_chx_PROM.txt | bedtools slop -b 250 -i stdin -g $SIZE > $WD3/regions/Sact_L_500bp.bed
awk -v OFS="\t" -v FS="\t" 'midpt=int(($28+$29)/2) {if($16=="Lon")print $27,midpt,midpt+1,$1,$3,$30}' $WD2/ls_fw_conf_chx_PROM.txt | bedtools slop -b 250 -i stdin -g $SIZE > $WD3/regions/Lact_S_500bp.bed
awk -v OFS="\t" -v FS="\t" 'midpt=int(($28+$29)/2) {if($16=="Son")print $27,midpt,midpt+1,$1,$3,$30}' $WD2/ls_fw_conf_chx_PROM.txt | bedtools slop -b 250 -i stdin -g $SIZE > $WD3/regions/Sact_S_500bp.bed
# make fasta files
GENOME=files/anno/genome.fa
bedtools getfasta -fo $WD3/fasta/Lact_L.fa -fi $GENOME -bed $WD3/regions/Lact_L_500bp.bed -nameOnly
bedtools getfasta -fo $WD3/fasta/Lact_S.fa -fi $GENOME -bed $WD3/regions/Lact_S_500bp.bed -nameOnly
bedtools getfasta -fo $WD3/fasta/Sact_L.fa -fi $GENOME -bed $WD3/regions/Sact_L_500bp.bed -nameOnly
bedtools getfasta -fo $WD3/fasta/Sact_S.fa -fi $GENOME -bed $WD3/regions/Sact_S_500bp.bed -nameOnly
# find Motifs
MOTIFS=files/data/vert_known.motifs
findMotifs.pl $WD3/fasta/Lact_L.fa fasta $WD3/Lact/ -fasta $WD3/fasta/Lact_S.fa -mknown $MOTIFS -p 6
findMotifs.pl $WD3/fasta/Sact_S.fa fasta $WD3/Sact/ -fasta $WD3/fasta/Sact_L.fa -mknown $MOTIFS -p 6

MOTIF

exit
