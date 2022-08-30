# figure 3 commands

TAB=$'\t'

GENOME=files/anno/genome.fa
CHRSIZES=files/anno/chr_sizes.txt
CHAIN=files/chains/xlL_xlS_med.liftOver.chain
ATAC_BED=files/beds/s9_atac_hq.bed
K27_BED=files/beds/s8_k27ac_hq.bed
K27AC_BW=files/beds/s8_k27ac_hq.bw
NOAB_BW=files/bigwigs/pool_noab_hq.bw
HIGH=files/data/distal_k27high_chr.bed
LOW=files/data/distal_k27low_chr.bed
ATAC_BW=files/bigwigs/s9_atac_hq.bw
ATAC8_BW=files/bigwigs/s8_atac_hq.bw
LPROM=prom/heatmaps/regions/ls_fw_conf_chx_Lprom_noscaff.bed
SPROM=prom/heatmaps/regions/ls_fw_conf_chx_Sprom_noscaff.bed
GO=files/anno/xb_go_terms.txt
HOMEO=files/data/homeolog_pairs_expanded.txt


:<<"HEATMAPS"
# make directories
mkdir -p enh/heatmap
mkdir -p enh/heatmap/regions
mkdir -p enh/heatmap/matrix
mkdir -p enh/heatmap/pdfs
mkdir -p enh/heatmap/bigwigs
#
# pull out enhancer pairs
# sort by k27ac
grep "LS-" lift/combined/all_enh_homeologs.txt | grep distal | sort -k1,1 | cut -f1-5 | awk -v OFS="\t" -v FS="\t" '{print $2,$3,$4,$5,$1}' | bedtools slop -l 0 -r 1 -i stdin -g $CHRSIZES > enh/heatmap/regions/both_on_L.bed
grep "LS-" lift/combined/all_enh_homeologs.txt | grep distal | sort -k1,1 | cut -f1,6-9 | awk -v OFS="\t" -v FS="\t" '{print $2,$3,$4,$5,$1}' | bedtools slop -l 0 -r 1 -i stdin -g $CHRSIZES > enh/heatmap/regions/both_on_S.bed
grep "L-" lift/combined/all_enh_homeologs.txt | grep distal | sort -k1,1 | cut -f1-5 | awk -v OFS="\t" -v FS="\t" '{print $2,$3,$4,$5,$1}' | bedtools slop -l 0 -r 1 -i stdin -g $CHRSIZES > enh/heatmap/regions/L_on_L.bed
grep "L-" lift/combined/all_enh_homeologs.txt | grep distal | sort -k1,1 | cut -f1,6-9 | awk -v OFS="\t" -v FS="\t" '{print $2,$3,$4,$5,$1}' | bedtools slop -l 0 -r 1 -i stdin -g $CHRSIZES > enh/heatmap/regions/L_on_S.bed
grep "S-NA" lift/combined/all_enh_homeologs.txt | grep distal | sort -k1,1 | cut -f1-5 | awk -v OFS="\t" -v FS="\t" '{print $2,$3,$4,$5,$1}' | bedtools slop -l 0 -r 1 -i stdin -g $CHRSIZES > enh/heatmap/regions/S_on_L.bed
grep "S-NA" lift/combined/all_enh_homeologs.txt | grep distal | sort -k1,1 | cut -f1,6-9 | awk -v OFS="\t" -v FS="\t" '{print $2,$3,$4,$5,$1}' | bedtools slop -l 0 -r 1 -i stdin -g $CHRSIZES > enh/heatmap/regions/S_on_S.bed
#
# sort by s8 k27ac signal
join -j5 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,2.1,2.2,2.3,2.4,2.5,2.6 <(bedtools coverage -counts -a enh/heatmap/regions/both_on_L.bed -b $K27_BED | sort -k5,5) <(bedtools coverage -counts -a enh/heatmap/regions/both_on_S.bed -b $K27_BED | sort -k5,5) | awk -v OFS="\t" -v FS="\t" '{print $1,$2,$3,$4,$5,$7,$8,$9,$10,$11,int(($6+$12)/2)}' | sort -k11,11nr | cut -f1-5,11 > enh/heatmap/regions/both_on_L.k27Sort.bed
join -j5 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,2.1,2.2,2.3,2.4,2.5,2.6 <(bedtools coverage -counts -a enh/heatmap/regions/both_on_L.bed -b $K27_BED | sort -k5,5) <(bedtools coverage -counts -a enh/heatmap/regions/both_on_S.bed -b $K27_BED | sort -k5,5) | awk -v OFS="\t" -v FS="\t" '{print $1,$2,$3,$4,$5,$7,$8,$9,$10,$11,int(($6+$12)/2)}' | sort -k11,11nr | cut -f6-11 > enh/heatmap/regions/both_on_S.k27Sort.bed
join -j5 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,2.1,2.2,2.3,2.4,2.5 <(bedtools coverage -counts -a enh/heatmap/regions/L_on_L.bed -b $K27_BED | sort -k5,5) <(sort -k5,5 enh/heatmap/regions/L_on_S.bed) | sort -k6,6nr | cut -f1-6 > enh/heatmap/regions/L_on_L.k27Sort.bed
join -j5 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,2.1,2.2,2.3,2.4,2.5 <(bedtools coverage -counts -a enh/heatmap/regions/L_on_L.bed -b $K27_BED | sort -k5,5) <(sort -k5,5 enh/heatmap/regions/L_on_S.bed) | sort -k6,6nr | cut -f6-11 | awk -v OFS="\t" -v FS="\t" '{print $2,$3,$4,$5,$6,$1}' > enh/heatmap/regions/L_on_S.k27Sort.bed
join -j5 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,2.1,2.2,2.3,2.4,2.5,2.6 <(sort -k5,5 enh/heatmap/regions/S_on_L.bed) <(bedtools coverage -counts -a enh/heatmap/regions/S_on_S.bed -b $K27_BED | sort -k5,5) | sort -k11,11nr | cut -f1-5,11 > enh/heatmap/regions/S_on_L.k27Sort.bed
join -j5 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,2.1,2.2,2.3,2.4,2.5,2.6 <(sort -k5,5 enh/heatmap/regions/S_on_L.bed) <(bedtools coverage -counts -a enh/heatmap/regions/S_on_S.bed -b $K27_BED | sort -k5,5) | sort -k11,11nr | cut -f6-11  > enh/heatmap/regions/S_on_S.k27Sort.bed
#
# calculate enrichment bigwig
bigwigCompare --bigwig1 $K27AC_BW --bigwig2 $NOAB_BW --operation ratio --pseudocount 0.1 --binSize 200 -o enh/heatmap/bigwigs/k27ac_st8_enrich_200bp.bw -of bigwig -p 8
bigwigCompare --bigwig1 $K27AC_BW --bigwig2 $NOAB_BW --operation ratio --pseudocount 0.1 --binSize 50 -o enh/heatmap/bigwigs/k27ac_st8_enrich_50bp.bw -of bigwig -p 6 --skipZeroOverZero
K27AC_FC=enh/heatmap/bigwigs/k27ac_st8_enrich_200bp.bw
#
# plot heatmap with s9 atac signal and s8 k27ac enrichment
LONL=enh/heatmap/regions/L_on_L.k27Sort.bed
LONS=enh/heatmap/regions/L_on_S.k27Sort.bed
BOTHL=enh/heatmap/regions/both_on_L.k27Sort.bed
BOTHS=enh/heatmap/regions/both_on_S.k27Sort.bed
SONL=enh/heatmap/regions/S_on_L.k27Sort.bed
SONS=enh/heatmap/regions/S_on_S.k27Sort.bed
computeMatrix reference-point --referencePoint center -a 5000 -b 5000 -R $HIGH $LOW -S $K27AC_FC --binSize 100 --sortRegions keep --missingDataAsZero -o enh/heatmap/matrix/s8_k27ac_hq_ratio.k27sort.matrix -p8
computeMatrix reference-point --referencePoint center -a 5000 -b 5000 -R $LONL $LONS $BOTHL $BOTHS $SONL $SONS -S $K27AC_FC --binSize 100 --sortRegions keep --missingDataAsZero -o enh/heatmap/matrix/LvS_s8_k27ac_hq_ratio.k27sort.matrix -p8
plotHeatmap -m enh/heatmap/matrix/s8_k27ac_hq_ratio.k27sort.matrix -o enh/heatmap/pdfs/s8_k27ac_hq_ratio.k27sort.pdf --sortRegions keep --dpi 1000 --heatmapWidth 2 --heatmapHeight 7 --xAxisLabel '' --regionsLabel 'k27ac high' 'k27ac low' --samplesLabel 'stage8 K27AC HQ ratio' --colorList white,darkmagenta --zMin 1 --zMax 5
plotHeatmap -m enh/heatmap/matrix/LvS_s8_k27ac_hq_ratio.k27sort.matrix -o enh/heatmap/pdfs/LvS_s8_k27ac_hq_ratio.k27sort.pdf --sortRegions keep --dpi 1000 --heatmapWidth 2 --heatmapHeight 7 --xAxisLabel '' --regionsLabel 'Lon L' 'Lon S' 'both L' 'both S' 'Son L' 'Son S' --samplesLabel 'stage8 K27AC HQ ratio' --colorList white,darkmagenta --zMin 1 --zMax 5
#
computeMatrix reference-point --referencePoint center -a 5000 -b 5000 -R $HIGH $LOW -S $ATAC_BW --binSize 100 --sortRegions keep --missingDataAsZero -o enh/heatmap/matrix/s9_open_hq.k27sort.matrix -p8
computeMatrix reference-point --referencePoint center -a 5000 -b 5000 -R $LONL $LONS $BOTHL $BOTHS $SONL $SONS -S $ATAC_BW --binSize 100 --sortRegions keep --missingDataAsZero -o heatmap/matrix/LvS_s9_open_hq.k27sort.matrix -p8
plotHeatmap -m enh/heatmap/matrix/s9_open_hq.k27sort.matrix -o enh/heatmap/pdfs/s9_open_hq.k27sort.pdf --sortRegions keep --dpi 1000 --heatmapWidth 2 --heatmapHeight 7 --xAxisLabel '' --regionsLabel 'k27ac high' 'k27ac low' --samplesLabel 'stage9 ATAC HQ' --colorList white,olive --zMax 0.5
plotHeatmap -m enh/heatmap/matrix/LvS_s9_open_hq.k27sort.matrix -o enh/heatmap/pdfs/LvS_s9_open_hq.k27sort.pdf --sortRegions keep --dpi 1000 --heatmapWidth 2 --heatmapHeight 7 --xAxisLabel '' --regionsLabel 'Lon L' 'Lon S' 'both L' 'both S' 'Son L' 'Son S' --samplesLabel 'stage9 ATAC HQ' --colorList white,olive --zMax 0.5
#
computeMatrix reference-point --referencePoint center -a 5000 -b 5000 -R $HIGH $LOW -S $ATAC8_BW --binSize 100 --sortRegions keep --missingDataAsZero -o heatmap/matrix/s8_open_hq.k27sort.matrix -p8
plotHeatmap -m enh/heatmap/matrix/s8_open_hq.k27sort.matrix -o enh/heatmap/pdfs/s8_open_hq.k27sort.pdf --sortRegions keep --dpi 1000 --heatmapWidth 2 --heatmapHeight 7 --xAxisLabel '' --regionsLabel 'k27ac high' 'k27ac low' --samplesLabel 'stage8 ATAC HQ' --colorList white,olive --zMax 0.5
#
# quantify LS differences in coverage
mkdir -p enh/quant
mkdir -p enh/quant/regions
# calculate k27 and atac coverage
bedtools coverage -counts -a $BOTHL -b $K27_BED | bedtools coverage -counts -a stdin -b $ATAC_BED | sort -k5,5 > enh/quant/regions/both_on_L.k27atac.bed
bedtools coverage -counts -a $BOTHS -b $K27_BED | bedtools coverage -counts -a stdin -b $ATAC_BED | sort -k5,5 > enh/quant/regions/both_on_S.k27atac.bed
bedtools coverage -counts -a $LONL -b $K27_BED | bedtools coverage -counts -a stdin -b $ATAC_BED | sort -k5,5 > enh/quant/regions/L_on_L.k27atac.bed
awk -v OFS="\t" -v FS="\t" 'midpt=int(($3+$2)/2) {print $1,midpt,midpt+1,$4,$5}' $LONS | bedtools slop -b 250 -i stdin -g $CHRSIZES | bedtools coverage -counts -a stdin -b $K27_BED | bedtools coverage -counts -a stdin -b $ATAC_BED | sort -k5,5 > enh/quant/regions/L_on_S.k27atac.bed
bedtools coverage -counts -a $SONS -b $K27_BED | bedtools coverage -counts -a stdin -b $ATAC_BED | sort -k5,5 > enh/quant/regions/S_on_S.k27atac.bed
awk -v OFS="\t" -v FS="\t" 'midpt=int(($3+$2)/2) {print $1,midpt,midpt+1,$4,$5}' $SONL | bedtools slop -b 250 -i stdin -g $CHRSIZES | bedtools coverage -counts -a stdin -b $K27_BED | bedtools coverage -counts -a stdin -b $ATAC_BED | sort -k5,5 > enh/quant/regions/S_on_L.k27atac.bed
# normalize to RPM
awk -v OFS="\t" -v FS="\t" 'k27pm=1000000/23598517, atacpm=1000000/14663728 {print $0,$6*k27pm,$7*atacpm}' enh/quant/regions/both_on_L.k27atac.bed > enh/quant/regions/both_on_L.rpm.bed
awk -v OFS="\t" -v FS="\t" 'k27pm=1000000/23598517, atacpm=1000000/14663728 {print $0,$6*k27pm,$7*atacpm}' enh/quant/regions/both_on_S.k27atac.bed > enh/quant/regions/both_on_S.rpm.bed
awk -v OFS="\t" -v FS="\t" 'k27pm=1000000/23598517, atacpm=1000000/14663728 {print $0,$6*k27pm,$7*atacpm}' enh/quant/regions/L_on_L.k27atac.bed > enh/quant/regions/L_on_L.rpm.bed
awk -v OFS="\t" -v FS="\t" 'k27pm=1000000/23598517, atacpm=1000000/14663728 {print $0,$6*k27pm,$7*atacpm}' enh/quant/regions/L_on_S.k27atac.bed > enh/quant/regions/L_on_S.rpm.bed
awk -v OFS="\t" -v FS="\t" 'k27pm=1000000/23598517, atacpm=1000000/14663728 {print $0,$6*k27pm,$7*atacpm}' enh/quant/regions/S_on_L.k27atac.bed > enh/quant/regions/S_on_L.rpm.bed
awk -v OFS="\t" -v FS="\t" 'k27pm=1000000/23598517, atacpm=1000000/14663728 {print $0,$6*k27pm,$7*atacpm}' enh/quant/regions/S_on_S.k27atac.bed > enh/quant/regions/S_on_S.rpm.bed
# plot distributions and stat test
# t-test of log2 rpm - k27;atac - L;LS;S
Rscript k27atac_quant.R

HEATMAPS

:<<"SUPP"

# plot heatmaps for  st10 and st12 atac and k27ac data
# sorted by s8 k27ac
K27AC_S10_BW=files/bigwigs/s10.5_k27ac_hq.bw
S10_BW=files/bigwigs/s10.5_atac_hq.bw
S12_BW=files/bigwigs/s12_atac_hq.bw
computeMatrix reference-point --referencePoint center -a 5000 -b 5000 -R $LONL $LONS $BOTHL $BOTHS $SONL $SONS -S $K27AC_S10_BW --binSize 100 --sortRegions keep --outFileSortedRegions enh/heatmap/regions/enh_s10_k27acsort.bed --missingDataAsZero -o enh/heatmap/matrix/s10_k27ac.s10k27acsort.matrix -p8
computeMatrix reference-point --referencePoint center -a 5000 -b 5000 -R enh/heatmap/regions/enh_s10_k27acsort.bed -S $S10_BW --binSize 100 --sortRegions keep --missingDataAsZero -o enh/heatmap/matrix/s10_pool.s10k27sort.matrix -p8
computeMatrix reference-point --referencePoint center -a 5000 -b 5000 -R enh/heatmap/regions/enh_s10_k27acsort.bed -S $S12_BW --binSize 100 --sortRegions keep --missingDataAsZero -o enh/heatmap/matrix/s12_pool.s10k27sort.matrix -p8
#
plotHeatmap -m enh/heatmap/matrix/s10_k27ac.s10k27acsort.matrix -o enh/heatmap/pdfs/s10_k27ac.s10k27acsort.pdf --sortRegions keep --dpi 1000 --heatmapWidth 2 --heatmapHeight 7 --xAxisLabel '' --regionsLabel 'Lon L' 'Lon S' 'both L' 'both S' 'Son L' 'Son S' --samplesLabel 'stage 10 K27ac HQ' --colorList white,darkmagenta
plotHeatmap -m enh/heatmap/matrix/s10_pool.s10k27sort.matrix -o enh/heatmap/pdfs/s10_pool.s10k27acsort.pdf --sortRegions keep --dpi 1000 --heatmapWidth 2 --heatmapHeight 7 --xAxisLabel '' --regionsLabel 'Lon L' 'Lon S' 'both L' 'both S' 'Son L' 'Son S' --samplesLabel 'stage 10 ATAC HQ' --colorList white,olive --zMax 4.5
plotHeatmap -m enh/heatmap/matrix/s12_pool.s10k27sort.matrix -o enh/heatmap/pdfs/s12_pool.s10k27acsort.pdf --sortRegions keep --dpi 1000 --heatmapWidth 2 --heatmapHeight 7 --xAxisLabel '' --regionsLabel 'Lon L' 'Lon S' 'both L' 'both S' 'Son L' 'Son S' --samplesLabel 'stage 12 ATAC HQ' --colorList white,olive --zMax 4.5

SUPP


:<<"DISTANCE"

# calculate enhancer density around activated homeologs
# make directories and move active proms from prom
mkdir -p enh/distance/bigwigs
mkdir -p enh/distance/count_enh
mkdir -p enh/distance/matrix
mkdir -p enh/distance/matrix_R
mkdir -p enh/distance/pdfs
mkdir -p enh/distance/regions

#
# make 500bp centered prom regions
awk -v OFS="\t" -v FS="\t" "{print $1,$2,$3,$10,$4,$5,$6,$7,$8,$9}" tss/TSS_nonzero_anno > tss/TSS_nonzero_anno.bed
TSS=tss/TSS_nonzero_anno.bed

bedtools slop -b 500 -i $TSS -g $CHRSIZES > enh/distance/regions/tss_subs.bed

#
# pull out enhancers
grep "LS-" lift/combined/all_enh_homeologs.txt | grep distal | cut -f2-5 | sort -k1,1 -k2,2n > enh/distance/regions/both_L_enh.bed
grep "LS-" lift/combined/all_enh_homeologs.txt | grep distal | cut -f6-9 | sort -k1,1 -k2,2n > enh/distance/regions/both_S_enh.bed
grep "L-" lift/combined/all_enh_homeologs.txt | grep distal | cut -f2-5 | sort -k1,1 -k2,2n > enh/distance/regions/Lon_L_enh.bed
grep "S-NA" lift/combined/all_enh_homeologs.txt | grep distal | cut -f6-9 | sort -k1,1 -k2,2n > enh/distance/regions/Son_S_enh.bed

#
# identify regions of 2-fold enriched k27ac
bigWigToBedGraph enh/heatmap/bigwigs/k27ac_st8_ratio.bw enh/heatmap/bigwigs/k27ac_st8_ratio.bedGraph
awk -v OFS="\t" -v FS="\t" '{if($4 >= 2) print $0}' enh/heatmap/bigwigs/k27ac_st8_ratio.bedGraph | bedtools merge -c 4 -o mean -i stdin > enh/k27ac_st8_2fold_peaks.bed

K27_2F=enh/k27ac_st8_2fold_peaks.bed

#
# intersect to find which enhancers have at least 2-fold k27ac enrichment
bedtools intersect -wa -wb -a $K27_2F -b enh/distance/regions/both_L_enh.bed | cut -f1-3,8 | cat enh/distance/regions/both_L_enh.bed - | sort -k4,4 | bedtools groupby -i stdin -g 4 -c 1,2,3 -o distinct,min,max | awk -v OFS="\t" -v FS="\t" '{print $2,$3,$4,$1}' | sort -k1,1 -k2,2n > enh/distance/regions/both_L_k27.bed
bedtools intersect -wa -wb -a $K27_2F -b enh/distance/regions/both_S_enh.bed | cut -f1-3,8 | cat enh/distance/regions/both_S_enh.bed - | sort -k4,4 | bedtools groupby -i stdin -g 4 -c 1,2,3 -o distinct,min,max | awk -v OFS="\t" -v FS="\t" '{print $2,$3,$4,$1}' | sort -k1,1 -k2,2n > enh/distance/regions/both_S_k27.bed
bedtools intersect -wa -wb -a $K27_2F -b enh/distance/regions/Lon_L_enh.bed | cut -f1-3,8 | cat enh/distance/regions/Lon_L_enh.bed - | sort -k4,4 | bedtools groupby -i stdin -g 4 -c 1,2,3 -o distinct,min,max | awk -v OFS="\t" -v FS="\t" '{print $2,$3,$4,$1}' | sort -k1,1 -k2,2n > enh/distance/regions/Lon_L_k27.bed
bedtools intersect -wa -wb -a $K27_2F -b enh/distance/regions/Son_S_enh.bed | cut -f1-3,8 | cat enh/distance/regions/Son_S_enh.bed - | sort -k4,4 | bedtools groupby -i stdin -g 4 -c 1,2,3 -o distinct,min,max | awk -v OFS="\t" -v FS="\t" '{print $2,$3,$4,$1}' | sort -k1,1 -k2,2n > enh/distance/regions/Son_S_k27.bed

#
# subtract out promoters and make faux-bedgraph
bedtools merge -i enh/distance/regions/both_L_k27.bed | bedtools subtract -a stdin -b enh/distance/regions/tss_subs.bed | awk -v OFS="\t" -v FS="\t" '{print $0,"1"}' > enh/distance/regions/both_L_prebw.bed
bedtools merge -i enh/distance/regions/both_S_k27.bed | bedtools subtract -a stdin -b enh/distance/regions/tss_subs.bed | awk -v OFS="\t" -v FS="\t" '{print $0,"1"}' > enh/distance/regions/both_S_prebw.bed
bedtools merge -i enh/distance/regions/Lon_L_k27.bed | bedtools subtract -a stdin -b enh/distance/regions/tss_subs.bed | awk -v OFS="\t" -v FS="\t" '{print $0,"1"}' > enh/distance/regions/Lon_L_prebw.bed
bedtools merge -i enh/distance/regions/Son_S_k27.bed | bedtools subtract -a stdin -b enh/distance/regions/tss_subs.bed | awk -v OFS="\t" -v FS="\t" '{print $0,"1"}' > enh/distance/regions/Son_S_prebw.bed

#
# make bigwigs
for i in enh/distance/regions/*_prebw.bed; do wigToBigWig $i $CHRSIZES enh/distance/bigwigs/`basename $i .bed`.bw; done


#
# both and diff separate and plotted together
# plot separately
# L regions
computeMatrix reference-point --referencePoint center -a 25000 -b 25000 -R $LPROM -S enh/distance/bigwigs/both_L_prebw.bw --averageTypeBins max --binSize 500 --sortRegions keep --missingDataAsZero -o enh/distance/matrix/enh_count_L_both.matrix -p6
computeMatrix reference-point --referencePoint center -a 25000 -b 25000 -R $LPROM -S enh/distance/bigwigs/Lon_L_prebw.bw --averageTypeBins max --binSize 500 --sortRegions keep --missingDataAsZero -o enh/distance/matrix/enh_count_L_L.matrix -p6

# S regions
computeMatrix reference-point --referencePoint center -a 25000 -b 25000 -R $SPROM -S enh/distance/bigwigs/both_S_prebw.bw --averageTypeBins max --binSize 500 --sortRegions keep --missingDataAsZero -o enh/distance/matrix/enh_count_S_both.matrix -p6
computeMatrix reference-point --referencePoint center -a 25000 -b 25000 -R $SPROM -S enh/distance/bigwigs/Son_S_prebw.bw --averageTypeBins max --binSize 500 --sortRegions keep --missingDataAsZero -o enh/distance/matrix/enh_count_S_S.matrix -p6

# plotheatmaps
# conserved
plotHeatmap -m  enh/distance/matrix/enh_count_L_both.matrix -o enh/distance/pdfs/enh_count_L_both.pdf --sortRegions keep --dpi 1000 --heatmapWidth 2 --heatmapHeight 7 --xAxisLabel '' --regionsLabel 'L proms' --samplesLabel 'con enh' --colorList white,black --zMax 1 --zMin 0
plotHeatmap -m  enh/distance/matrix/enh_count_S_both.matrix -o enh/distance/pdfs/enh_count_S_both.pdf --sortRegions keep --dpi 1000 --heatmapWidth 2 --heatmapHeight 7 --xAxisLabel '' --regionsLabel 'S proms' --samplesLabel 'con enh' --colorList white,black --zMax 1 --zMin 0

# differential
plotHeatmap -m  enh/distance/matrix/enh_count_L_L.matrix -o enh/distance/pdfs/enh_count_L_L.pdf --sortRegions keep --dpi 1000 --heatmapWidth 2 --heatmapHeight 7 --xAxisLabel '' --regionsLabel 'L proms' --samplesLabel 'diff enh' --colorList white,black --zMax 1 --zMin 0
plotHeatmap -m  enh/distance/matrix/enh_count_S_S.matrix -o enh/distance/pdfs/enh_count_S_S.pdf --sortRegions keep --dpi 1000 --heatmapWidth 2 --heatmapHeight 7 --xAxisLabel '' --regionsLabel 'S proms' --samplesLabel 'diff enh' --colorList white,black --zMax 1 --zMin 0

#
# count number of enhancers within 25kb (regions are defined as 1kb regions already centered on tss)
bedtools slop -b 24500 -i $LPROM -g $CHRSIZES | bedtools intersect -wao -a stdin -b enh/distance/regions/both_L_prebw.bed | bedtools groupby -i stdin -g 1,2,3,4,5,6,7,8 -c 13 -o sum | bedtools intersect -wao -a stdin -b enh/distance/regions/Lon_L_prebw.bed | bedtools groupby -i stdin -g 1,2,3,4,5,6,7,8,9 -c 14 -o sum > enh/distance/count_enh/L_counts.bed

bedtools slop -b 24500 -i $SPROM -g $CHRSIZES | bedtools intersect -wao -a stdin -b enh/distance/regions/both_S_prebw.bed | bedtools groupby -i stdin -g 1,2,3,4,5,6,7,8 -c 13 -o sum | bedtools intersect -wao -a stdin -b enh/distance/regions/Son_S_prebw.bed | bedtools groupby -i stdin -g 1,2,3,4,5,6,7,8,9 -c 14 -o sum > enh/distance/count_enh/S_counts.bed

join -j5 -t $'\t' -o 1.5,1.1,1.2,1.3,1.4,1.6,2.1,2.2,2.3,2.4,2.6,1.7,1.8,1.9,1.10,2.9,2.10 <(sort -k5,5 -t "${TAB}" enh/distance/count_enh/L_counts.bed) <(sort -k5,5 -t "${TAB}" enh/distance/count_enh/S_counts.bed) | sort -t "${TAB}" -k13,13nr > enh/distance/count_enh/total_enh_counts.txt

# make meta plots
Rscript distance_meta.R

# P values: cor.test for differential enh and conserved enh
#[1] 1.308308e-16
#[1] 0.2015453


DISTANCE

:<<"MOTIF"

# motif search both on L v S for control
# /sandbox/projects/xl_zga/analysis/enh
mkdir -p enh/homer/both_LvS
mkdir -p enh/homer/both_SvL
mkdir -p enh/homer/Lon
mkdir -p enh/homer/Son
mkdir -p enh/homer/prox/Lon_prox
mkdir -p enh/homer/prox/Son_prox
mkdir -p enh/homer/tss/Lon_tss
mkdir -p enh/homer/tss/Son_tss
# make midpt function
function BEDMIDPT () {
	awk -v OFS="\t" '{midpt=int(($3+$2)/2); print $1,midpt,midpt+1,$4}' $1	
}
# Single on regions
# proximals
# L
awk -v OFS="\t" '{print $2,$3,$4,$1}' lift/combined/L_on.txt | BEDMIDPT | grep proximal | bedtools intersect -v -a stdin -b <(bedtools slop -b 100 -i $TSS -g $CHRSIZES) > lift/combined/prox_L_on_L_1bp.bed
awk -v OFS="\t" '{print $6,$7,$8,$1}' lift/combined/L_on.txt | BEDMIDPT | grep proximal | bedtools intersect -v -a stdin -b <(bedtools slop -b 100 -i $TSS -g $CHRSIZES) > lift/combined/prox_L_on_S_1bp.bed
# S
awk -v OFS="\t" '{print $2,$3,$4,$1}' lift/combined/S_on.txt | BEDMIDPT | grep proximal | bedtools intersect -v -a stdin -b <(bedtools slop -b 100 -i $TSS -g $CHRSIZES) > lift/combined/prox_S_on_S_1bp.bed
awk -v OFS="\t" '{print $6,$7,$8,$1}' lift/combined/S_on.txt | BEDMIDPT | grep proximal | bedtools intersect -v -a stdin -b <(bedtools slop -b 100 -i $TSS -g $CHRSIZES) > lift/combined/prox_S_on_L_1bp.bed
# tss overlap
awk -v OFS="\t" '{print $2,$3,$4,$1}' lift/combined/L_on.txt | BEDMIDPT | grep proximal | bedtools intersect -u -a stdin -b <(bedtools slop -b 100 -i $TSS -g $CHRSIZES) > lift/combined/tss_L_on_L_1bp.bed
awk -v OFS="\t" '{print $6,$7,$8,$1}' lift/combined/L_on.txt | BEDMIDPT | grep proximal | bedtools intersect -u -a stdin -b <(bedtools slop -b 100 -i $TSS -g $CHRSIZES) > lift/combined/tss_L_on_S_1bp.bed
# S
awk -v OFS="\t" '{print $2,$3,$4,$1}' lift/combined/S_on.txt | BEDMIDPT | grep proximal | bedtools intersect -u -a stdin -b <(bedtools slop -b 100 -i $TSS -g $CHRSIZES) > lift/combined/tss_S_on_S_1bp.bed
awk -v OFS="\t" '{print $6,$7,$8,$1}' lift/combined/S_on.txt | BEDMIDPT | grep proximal | bedtools intersect -u -a stdin -b <(bedtools slop -b 100 -i $TSS -g $CHRSIZES) > lift/combined/tss_S_on_L_1bp.bed
#
# +/- 100bp for Homer motif finding
# make prox regions
bedtools slop -b 100 -g $CHRSIZES -i lift/combined/prox_L_on_L_1bp.bed | bedtools getfasta -nameOnly -fi $GENOME -bed - -fo enh/homer/prox/prox_L_on_L_200bp.fa
bedtools slop -b 100 -g $CHRSIZES -i lift/combined/prox_L_on_S_1bp.bed | bedtools getfasta -nameOnly -fi $GENOME -bed - -fo enh/homer/prox/prox_L_on_S_200bp.fa
bedtools slop -b 100 -g $CHRSIZES -i lift/combined/prox_S_on_L_1bp.bed | bedtools getfasta -nameOnly -fi $GENOME -bed - -fo enh/homer/prox/prox_S_on_L_200bp.fa
bedtools slop -b 100 -g $CHRSIZES -i lift/combined/prox_S_on_S_1bp.bed | bedtools getfasta -nameOnly -fi $GENOME -bed - -fo enh/homer/prox/prox_S_on_S_200bp.fa
# make tss regions
bedtools slop -b 100 -g $CHRSIZES -i lift/combined/tss_L_on_L_1bp.bed | bedtools getfasta -nameOnly -fi $GENOME -bed - -fo enh/homer/tss/tss_L_on_L_200bp.fa
bedtools slop -b 100 -g $CHRSIZES -i lift/combined/tss_L_on_S_1bp.bed | bedtools getfasta -nameOnly -fi $GENOME -bed - -fo enh/homer/tss/tss_L_on_S_200bp.fa
bedtools slop -b 100 -g $CHRSIZES -i lift/combined/tss_S_on_L_1bp.bed | bedtools getfasta -nameOnly -fi $GENOME -bed - -fo enh/homer/tss/tss_S_on_L_200bp.fa
bedtools slop -b 100 -g $CHRSIZES -i lift/combined/tss_S_on_S_1bp.bed | bedtools getfasta -nameOnly -fi $GENOME -bed - -fo enh/homer/tss/tss_S_on_S_200bp.fa
# distal motif search
findMotifs.pl lift/combined/distal_L_on_L_200bp.fa fasta homer/Lon/ -fasta lift/combined/distal_L_on_S_200bp.fa -mknown files/data/vert_known.motifs -p 6
findMotifs.pl lift/combined/distal_S_on_S_200bp.fa fasta homer/Son/ -fasta lift/combined/distal_S_on_L_200bp.fa -mknown files/data/vert_known.motifs -p 6
findMotifs.pl lift/combined/distal_both_on_L_200bp.fa fasta homer/both_LvS/ -fasta lift/combined/distal_both_on_S_200bp.fa -mknown files/data/vert_known.motifs -p 6
findMotifs.pl lift/combined/distal_both_on_S_200bp.fa fasta homer/both_SvL/ -fasta lift/combined/distal_both_on_L_200bp.fa -mknown files/data/vert_known.motifs -p 6
# proximal motif search
findMotifs.pl enh/homer/prox/prox_L_on_L_200bp.fa fasta enh/homer/prox/Lon_prox/ -fasta enh/homer/prox/prox_L_on_S_200bp.fa -mknown files/data/vert_known.motifs -p 6
findMotifs.pl enh/homer/prox/prox_S_on_S_200bp.fa fasta enh/homer/prox/Son_prox/ -fasta enh/homer/prox/prox_S_on_L_200bp.fa -mknown files/data/vert_known.motifs -p 6
# tss motif search
findMotifs.pl enh/homer/tss/tss_L_on_L_200bp.fa fasta enh/homer/tss/Lon_tss/ -fasta enh/homer/tss/tss_L_on_S_200bp.fa -mknown files/data/vert_known.motifs -p 6
findMotifs.pl enh/homer/tss/tss_S_on_S_200bp.fa fasta enh/homer/tss/Son_tss/ -fasta enh/homer/tss/tss_S_on_L_200bp.fa -mknown files/data/vert_known.motifs -p 6

MOTIF

:<<"MAT_TF"

# pull out genes annotated as transcription factors and plot maternal expression
# GO:0003700 = DNA-binding transcription factor activity
# GO:0000981 = DNA-binding transcription factor activity, RNA polymerase II-specific
# download go annotations from xenbase
mkdir -p enh/tfs
grep -w GO:0003700 $GO > enh/tfs/GO0003700.txt
grep -w GO:0000981 $GO > enh/tfs/GO0000981.txt
cat enh/tfs/GO0003700.txt enh/tfs/GO0000981.txt | sort -k2,2 | uniq > enh/tfs/go_tfs.txt
# pull out L and S genes based on homeo_expanded doc
cut -f5-7 $HOMEO | sort -k1,1 > enh/tfs/Lhomeo.txt
cut -f8-10 $HOMEO | sort -k1,1 > enh/tfs/Shomeo.txt
cut -f1,39 rna/exon_tpm.txt > enh/tfs/mat_exp.txt
join -1 2 -2 1 -t $'\t' -o 2.2,2.3,1.4 enh/tfs/go_tfs.txt enh/tfs/Lhomeo.txt | sort -k1,1 | uniq | join -j1 -t $'\t' -o 2.1,2.2,1.2,2.3 <(sort -k1,1 enh/tfs/mat_exp.txt) - > enh/tfs/L_tf_exp.txt
join -1 2 -2 1 -t $'\t' -o 2.2,2.3,1.4 enh/tfs/go_tfs.txt enh/tfs/Shomeo.txt | sort -k1,1 | uniq | join -j1 -t $'\t' -o 2.1,2.2,1.2,2.3 <(sort -k1,1 enh/tfs/mat_exp.txt) - > enh/tfs/S_tf_exp.txt
# beeswarm in R
Rscript tf_bees.R

MAT_TF

:<<LANDMARK

# liftOver check
# make both, L, and S bed files
mkdir -p enh/landmark/closest
mkdir -p enh/landmark/agg
mkdir -p enh/landmark/match
LFC=rna/homeolog_lfc_all.txt
SELECT=tss/TSS_select.bed
# both - 1312
grep 'LS-' lift/combined/all_enh_homeologs.txt | grep distal | cut -f1-4 | awk -v OFS="\t" -v FS="\t" '{print $2,$3,$4,$1}' > enh/landmark/both_on_L.bed
grep 'LS-' lift/combined/all_enh_homeologs.txt | grep distal | cut -f1,6-8 | awk -v OFS="\t" -v FS="\t" '{print $2,$3,$4,$1}' > enh/landmark/both_on_S.bed
# L - 2561
grep 'L-' lift/combined/all_enh_homeologs.txt | grep distal | cut -f1-4 | awk -v OFS="\t" -v FS="\t" '{print $2,$3,$4,$1}' > enh/landmark/L_on_L.bed
grep 'L-' lift/combined/all_enh_homeologs.txt | grep distal | cut -f1,6-8 | awk -v OFS="\t" -v FS="\t" '{print $2,$3,$4,$1}' > enh/landmark/L_on_S.bed
# S - 1776
grep 'S-NA' lift/combined/all_enh_homeologs.txt | grep distal | cut -f1-4 | awk -v OFS="\t" -v FS="\t" '{print $2,$3,$4,$1}' > enh/landmark/S_on_L.bed
grep 'S-NA' lift/combined/all_enh_homeologs.txt | grep distal | cut -f1,6-8 | awk -v OFS="\t" -v FS="\t" '{print $2,$3,$4,$1}' > enh/landmark/S_on_S.bed
# make homeo bed file
# L
awk -v OFS="\t" -v FS="\t" '{if($2 == "LS") print $0}' $LFC | cut -f3,4 | join -1 1 -2 4 -t $'\t' -o 2.1,2.2,2.3,2.4,1.2,2.6 <(sort -k1,1 -) <(sort -t "${TAB}" -k4,4 $SELECT) | cut -f1-6 > enh/landmark/homeo_L.bed
# S
awk -v OFS="\t" -v FS="\t" '{if($2 == "LS") print $0}' $LFC | cut -f26,27 | join -1 1 -2 4 -t $'\t' -o 2.1,2.2,2.3,2.4,1.2,2.6 <(sort -k1,1 -) <(sort -t "${TAB}" -k4,4 $SELECT) | cut -f1-6 > enh/landmark/homeo_S.bed
# put together
cat enh/landmark/homeo_L.bed enh/landmark/homeo_S.bed | sort -k4,4 > enh/landmark/homeos.bed
rm enh/landmark/homeo_L.bed enh/landmark/homeo_S.bed
#
# identify 5 closest upstream and downstream genes with respect to enhancer. Take all ties
# both - upstream
bedtools closest -k 5 -D a -t all -id -a <(sort -k1,1 -k2,2n enh/landmark/both_on_L.bed) -b <(sort -k1,1 -k2,2n enh/landmark/homeos.bed) | sed s/'\.L'//g | sed s/'\.S'//g > enh/landmark/closest/both_on_L_us.bed
bedtools closest -k 5 -D a -t all -id -a <(sort -k1,1 -k2,2n enh/landmark/both_on_S.bed) -b <(sort -k1,1 -k2,2n enh/landmark/homeos.bed) | sed s/'\.L'//g | sed s/'\.S'//g > enh/landmark/closest/both_on_S_us.bed
# L - upstream
bedtools closest -k 5 -D a -t all -id -a <(sort -k1,1 -k2,2n enh/landmark/L_on_L.bed) -b <(sort -k1,1 -k2,2n enh/landmark/homeos.bed) | sed s/'\.L'//g | sed s/'\.S'//g > enh/landmark/closest/L_on_L_us.bed
bedtools closest -k 5 -D a -t all -id -a <(sort -k1,1 -k2,2n enh/landmark/L_on_S.bed) -b <(sort -k1,1 -k2,2n enh/landmark/homeos.bed) | sed s/'\.L'//g | sed s/'\.S'//g > enh/landmark/closest/L_on_S_us.bed
# S - upstream
bedtools closest -k 5 -D a -t all -id -a <(sort -k1,1 -k2,2n enh/landmark/S_on_L.bed) -b <(sort -k1,1 -k2,2n enh/landmark/homeos.bed) | sed s/'\.L'//g | sed s/'\.S'//g > enh/landmark/closest/S_on_L_us.bed
bedtools closest -k 5 -D a -t all -id -a <(sort -k1,1 -k2,2n enh/landmark/S_on_S.bed) -b <(sort -k1,1 -k2,2n enh/landmark/homeos.bed) | sed s/'\.L'//g | sed s/'\.S'//g > enh/landmark/closest/S_on_S_us.bed
# both - downstream
bedtools closest -k 5 -D a -t all -iu -a <(sort -k1,1 -k2,2n enh/landmark/both_on_L.bed) -b <(sort -k1,1 -k2,2n enh/landmark/homeos.bed) | sed s/'\.L'//g | sed s/'\.S'//g > enh/landmark/closest/both_on_L_ds.bed
bedtools closest -k 5 -D a -t all -iu -a <(sort -k1,1 -k2,2n enh/landmark/both_on_S.bed) -b <(sort -k1,1 -k2,2n enh/landmark/homeos.bed) | sed s/'\.L'//g | sed s/'\.S'//g > enh/landmark/closest/both_on_S_ds.bed
# L - downstream
bedtools closest -k 5 -D a -t all -iu -a <(sort -k1,1 -k2,2n enh/landmark/L_on_L.bed) -b <(sort -k1,1 -k2,2n enh/landmark/homeos.bed) | sed s/'\.L'//g | sed s/'\.S'//g > enh/landmark/closest/L_on_L_ds.bed
bedtools closest -k 5 -D a -t all -iu -a <(sort -k1,1 -k2,2n enh/landmark/L_on_S.bed) -b <(sort -k1,1 -k2,2n enh/landmark/homeos.bed) | sed s/'\.L'//g | sed s/'\.S'//g > enh/landmark/closest/L_on_S_ds.bed
# S - downstream
bedtools closest -k 5 -D a -t all -iu -a <(sort -k1,1 -k2,2n enh/landmark/S_on_L.bed) -b <(sort -k1,1 -k2,2n enh/landmark/homeos.bed) | sed s/'\.L'//g | sed s/'\.S'//g > enh/landmark/closest/S_on_L_ds.bed
bedtools closest -k 5 -D a -t all -iu -a <(sort -k1,1 -k2,2n enh/landmark/S_on_S.bed) -b <(sort -k1,1 -k2,2n enh/landmark/homeos.bed) | sed s/'\.L'//g | sed s/'\.S'//g > enh/landmark/closest/S_on_S_ds.bed
# aggregate up and downstream
Rscript landmark/aggregate_homeo.R
# reconstruct both, L, and S files for input into python script
# both
join -1 4 -2 1 -t $'\t' -o 1.1,1.2,1.3,1.4,2.2,2.3 <(sort -k4,4 enh/landmark/both_on_L.bed) <(sort -k1,1 enh/landmark/agg/bothL_agg.txt) | join -j 4 -t $'\t' -o 1.4,1.1,1.2,1.3,1.5,1.6,2.1,2.2,2.3 <(sort -t "${TAB}" -k4,4 -) <(sort -k4,4 enh/landmark/both_on_S.bed) | join -j 1 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.2,2.3 <(sort -t "${TAB}" -k1,1 -) <(sort -k1,1 enh/landmark/agg/bothS_agg.txt) > enh/landmark/agg/both_on_total.txt
# Lon
join -1 4 -2 1 -t $'\t' -o 1.1,1.2,1.3,1.4,2.2,2.3 <(sort -k4,4 enh/landmark/L_on_L.bed) <(sort -k1,1 enh/landmark/agg/LonL_agg.txt) | join -j 4 -t $'\t' -o 1.4,1.1,1.2,1.3,1.5,1.6,2.1,2.2,2.3 <(sort -t "${TAB}" -k4,4 -) <(sort -k4,4 enh/landmark/L_on_S.bed) | join -j 1 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.2,2.3 <(sort -t "${TAB}" -k1,1 -) <(sort -k1,1 enh/landmark/agg/LonS_agg.txt) > enh/landmark/agg/L_on_total.txt
# Son
join -1 4 -2 1 -t $'\t' -o 1.1,1.2,1.3,1.4,2.2,2.3 <(sort -k4,4 enh/landmark/S_on_L.bed) <(sort -k1,1 enh/landmark/agg/SonL_agg.txt) | join -j 4 -t $'\t' -o 1.4,1.1,1.2,1.3,1.5,1.6,2.1,2.2,2.3 <(sort -t "${TAB}" -k4,4 -) <(sort -k4,4 enh/landmark/S_on_S.bed) | join -j 1 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.2,2.3 <(sort -t "${TAB}" -k1,1 -) <(sort -k1,1 enh/landmark/agg/SonS_agg.txt) > enh/landmark/agg/S_on_total.txt
#
# input into python and convert HQ pairs to fasta format
# !!!! manually alter script for each set !!!!
python3 enh/landmark/match_homeo.py > enh/landmark/match/both_enh_match.txt
python3 enh/landmark/match_homeo.py > enh/landmark/match/L_enh_match.txt
python3 enh/landmark/match_homeo.py > enh/landmark/match/S_enh_match.txt
# calculate percentage correctly mapped (at least 1 up and downstream)
# barplot results
Rscript enh/landmark/landmark_bars.R

LANDMARK

exit

