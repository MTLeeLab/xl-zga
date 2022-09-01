# commands for generating pou/sox heatmaps and plots - fig4
# look at coverage over tf peaks

TAB=$'\t'

TPM=rna/exon_tpm.txt
HOMEO=rna/homeolog_lfc_all.txt
TSS_SELECT=tss/TSS_select.bed
POU_PEAK=files/data/pou_peaks.bed
SOX_PEAK=files/data/sox_peaks.bed
POU_COV=files/beds/s8_pou_hq.bed
SOX_COV=files/beds/s8_sox_hq.bed
NOAB_COV=files/beds/s8_tf_noab_hq.bed
POU_BW=files/bigwigs/s8_pou_hq.bw
SOX_BW=files/bigwigs/s8_sox_hq.bw
NOAB_POU=files/bigwigs/s8_noab_pou_hq.bw
NOAB_SOX=files/bigwigs/s8_noab_sox_hq.bw
SIZE=files/anno/chr_sizes.txt
POU_ENR=files/bigwigs/s8_pou_enrich.bw
SOX_ENR=files/bigwigs/s8_sox_enrich.bw
ENH=lift/combined/all_enh_homeologs.txt
GENOME=files/anno/genome.fa
K27_BG=files/data/distal_k27low.bed
LPROM=prom/heatmaps/regions/ls_fw_conf_chx_Lprom_noscaff.bed
SPROM=prom/heatmaps/regions/ls_fw_conf_chx_Sprom_noscaff.bed
CATS=rna/ls_fw_conf_chx_ps.txt
MOTIF=files/data/vert_known.motifs
PS_ENRICH=files/data/all_enh_homeologs_ps_enrich.txt
GENOME=files/anno/genome.fa

:<<TFPEAKS

# make output directories
mkdir -p tf/heatmaps/regions
mkdir -p tf/heatmaps/matrix
mkdir -p tf/heatmaps/pdfs

# intersect and take the lowest number of unique overlapping regions
# 1427 pousox regions
bedtools intersect -u -a <(cut -f1-3,5,6 $POU_PEAK) -b <(cut -f1-3,5,6 $SOX_PEAK) > tf/heatmaps/regions/pou_sox_regions.bed
# 3979 pou regions
bedtools intersect -v -a <(cut -f1-3,5,6 $POU_PEAK) -b <(cut -f1-3,5,6 $SOX_PEAK) > tf/heatmaps/regions/pou_regions.bed
# 24610 sox regions
bedtools intersect -v -b <(cut -f1-3,5,6 $POU_PEAK) -a <(cut -f1-3,5,6 $SOX_PEAK) > tf/heatmaps/regions/sox_regions.bed

# coverage over central regions for sort order
cut -f4,5 tf/heatmaps/regions/pou_regions.bed | sed s/'-'/"${TAB}"/g | sed s/':'/"${TAB}"/g | awk -v OFS="\t" -v FS="\t" '{midpt=int(($4+$3)/2); print $2,midpt,midpt+1,$1}' | bedtools slop -b 100 -i stdin -g $SIZE | sort -k4,4n > tf/heatmaps/regions/pou_regions.pouAscSort.bed
cut -f4,5 tf/heatmaps/regions/pou_sox_regions.bed | sed s/'-'/"${TAB}"/g | sed s/':'/"${TAB}"/g | awk -v OFS="\t" -v FS="\t" '{midpt=int(($4+$3)/2); print $2,midpt,midpt+1,$1}' | bedtools slop -b 100 -i stdin -g $SIZE | sort -k4,4nr > tf/heatmaps/regions/pou_sox_regions.pouDecSort.bed
cut -f4,5 tf/heatmaps/regions/sox_regions.bed | sed s/'-'/"${TAB}"/g | sed s/':'/"${TAB}"/g | awk -v OFS="\t" -v FS="\t" '{midpt=int(($4+$3)/2); print $2,midpt,midpt+1,$1}' | bedtools slop -b 100 -i stdin -g $SIZE | sort -k4,4nr > tf/heatmaps/regions/sox_regions.soxDecSort.bed
cat tf/heatmaps/regions/pou_regions.pouAscSort.bed tf/heatmaps/regions/pou_sox_regions.pouDecSort.bed tf/heatmaps/regions/sox_regions.soxDecSort.bed > tf/heatmaps/regions/tf_regions.psSort.bed

#
# computeMatrix on sorted regions
TF_PEAKS=tf/heatmaps/regions/tf_regions.psSort.bed
computeMatrix reference-point --referencePoint center -a 5000 -b 5000 -R $TF_PEAKS -S $POU_BW $NOAB_POU --binSize 125 --sortRegions keep --missingDataAsZero -o heatmaps/matrix/pou.pousoxsort.matrix -p8
computeMatrix reference-point --referencePoint center -a 5000 -b 5000 -R $TF_PEAKS -S $SOX_BW $NOAB_SOX --binSize 125 --sortRegions keep --missingDataAsZero -o heatmaps/matrix/sox.pousoxsort.matrix -p8

# plot heatmaps
plotHeatmap -m tf/heatmaps/matrix/pou.pousoxsort.matrix -o tf/heatmaps/pdfs/pou.pousoxsort.pdf --sortRegions keep --dpi 1000 --heatmapWidth 2 --heatmapHeight 7 --xAxisLabel '' --regionsLabel 'tf peaks' --samplesLabel 'pou' 'noab' --colorList white,maroon
plotHeatmap -m tf/heatmaps/matrix/sox.pousoxsort.matrix -o tf/heatmaps/pdfs/sox.pousoxsort.pdf --sortRegions keep --dpi 1000 --heatmapWidth 2 --heatmapHeight 7 --xAxisLabel '' --regionsLabel 'tf peaks' --samplesLabel 'sox' 'noab' --colorList white,darkgreen

TFPEAKS

:<<TFMOTIF

# motif call top 500 tf peaks
mkdir -p tf/homer/fasta
mkdir -p tf/homer/pou
mkdir -p tf/homer/sox

# make fasta files
sort -k5,5nr $POU_PEAK | head -500 | bedtools getfasta -nameOnly -fi $GENOME -fo tf/homer/fasta/pou_top500.fa -bed stdin
sort -k5,5nr $SOX_PEAK | head -500 | bedtools getfasta -nameOnly -fi $GENOME -fo tf/homer/fasta/sox_top500.fa -bed stdin
sort -k5,5nr $K27_BG | bedtools getfasta -nameOnly -fi $GENOME -fo tf/homer/fasta/distal_k27low.fa -bed stdin

# run homer
findMotifs.pl tf/homer/fasta/pou_top500.fa fasta tf/homer/pou/ -fasta tf/homer/fasta/distal_k27low.fa -mknown $MOTIF -p 6
findMotifs.pl tf/homer/fasta/sox_top500.fa fasta tf/homer/sox/ -fasta tf/homer/fasta/distal_k27low.fa -mknown $MOTIF -p 6

TFMOTIF

:<<ENH

# Pou/Sox enrichment over LvS enhancers
# make directories
mkdir -p LS_heatmaps/regions
mkdir -p LS_heatmaps/matrix
mkdir -p LS_heatmaps/pdfs
mkdir -p LS_heatmaps/bigwig

# sort regions by tf coverage
LONL=enh/heatmaps/regions/L_on_L.k27Sort.bed
LONS=enh/heatmaps/regions/L_on_S.k27Sort.bed
SONS=enh/heatmaps/regions/S_on_S.k27Sort.bed
SONL=enh/heatmaps/regions/S_on_L.k27Sort.bed
BOTHL=enh/heatmaps/regions/both_on_L.k27Sort.bed
BOTHS=enh/heatmaps/regions/both_on_S.k27Sort.bed

# calculate tf coverage over enh regions
join -j5 -t $'\t' -o 1.1,1.2,1.3,1.4,2.1,2.2,2.3,2.4,1.5 <(cut -f1-5 $LONL) <(cut -f1-5 $LONS) | sort -k1,1 -k2,2n | bedtools coverage -sorted -counts -a stdin -b $POU_COV | bedtools coverage -sorted -counts -a stdin -b $SOX_COV | awk -v OFS="\t" -v FS="\t" '{print $0,$10+$11}' | sort -k12,12nr | cut -f1-4,9,12 > tf/LS_heatmaps/regions/L_on_L.tfSort.bed
join -j5 -t $'\t' -o 1.1,1.2,1.3,1.4,2.1,2.2,2.3,2.4,1.5 <(cut -f1-5 $LONL) <(cut -f1-5 $LONS) | sort -k1,1 -k2,2n | bedtools coverage -sorted -counts -a stdin -b $POU_COV | bedtools coverage -sorted -counts -a stdin -b $SOX_COV | awk -v OFS="\t" -v FS="\t" '{print $0,$10+$11}' | sort -k12,12nr | cut -f5-9,12 > tf/LS_heatmaps/regions/L_on_S.tfSort.bed
join -j5 -t $'\t' -o 1.1,1.2,1.3,1.4,2.1,2.2,2.3,2.4,1.5 <(cut -f1-5 $SONS) <(cut -f1-5 $SONL) | sort -k1,1 -k2,2n | bedtools coverage -sorted -counts -a stdin -b $POU_COV | bedtools coverage -sorted -counts -a stdin -b $SOX_COV | awk -v OFS="\t" -v FS="\t" '{print $0,$10+$11}' | sort -k12,12nr | cut -f1-4,9,12 > tf/LS_heatmaps/regions/S_on_S.tfSort.bed
join -j5 -t $'\t' -o 1.1,1.2,1.3,1.4,2.1,2.2,2.3,2.4,1.5 <(cut -f1-5 $SONS) <(cut -f1-5 $SONL) | sort -k1,1 -k2,2n | bedtools coverage -sorted -counts -a stdin -b $POU_COV | bedtools coverage -sorted -counts -a stdin -b $SOX_COV | awk -v OFS="\t" -v FS="\t" '{print $0,$10+$11}' | sort -k12,12nr | cut -f5-9,12 > tf/LS_heatmaps/regions/S_on_L.tfSort.bed

# calculate coverage for both on enh separately then pair them
sort -k1,1 -k2,2n $BOTHL | cut -f1-5 | bedtools coverage -sorted -counts -a stdin -b $POU_COV | bedtools coverage -sorted -counts -a stdin -b $SOX_COV | awk -v OFS="\t" -v FS="\t" '{print $0,$6+$7}' | cut -f1-5,8 > tf/both_on_L.tmp
sort -k1,1 -k2,2n $BOTHS | cut -f1-5 | tf/bedtools coverage -sorted -counts -a stdin -b $POU_COV | bedtools coverage -sorted -counts -a stdin -b $SOX_COV | awk -v OFS="\t" -v FS="\t" '{print $0,$6+$7}' | cut -f1-5,8 > tf/both_on_S.tmp
join -j5 -t $'\t' <(sort -k5,5 tf/both_on_L.tmp) <(sort -k5,5 tf/both_on_S.tmp) | awk -v OFS="\t" -v FS="\t" '{print $0,$6+$11}' | awk -v OFS="\t" -v FS="\t" '{print $2,$3,$4,$5,$1,$12}' > tf/LS_heatmaps/regions/both_on_L.tfSort.bed
join -j5 -t $'\t' <(sort -k5,5 tf/both_on_L.tmp) <(sort -k5,5 tf/both_on_S.tmp) | awk -v OFS="\t" -v FS="\t" '{print $0,$6+$11}' | awk -v OFS="\t" -v FS="\t" '{print $7,$8,$9,$10,$1,$12}' > tf/LS_heatmaps/regions/both_on_S.tfSort.bed

# remove intermediate files
rm tf/both_on_L.tmp tf/both_on_S.tmp

#
# sort regions by 2-fold category and then tf coverage
LONL=tf/LS_heatmaps/regions/L_on_L.tfSort.bed
LONS=tf/LS_heatmaps/regions/L_on_S.tfSort.bed
SONS=tf/LS_heatmaps/regions/S_on_S.tfSort.bed
SONL=tf/LS_heatmaps/regions/S_on_L.tfSort.bed
BOTHL=tf/LS_heatmaps/regions/both_on_L.tfSort.bed
BOTHS=tf/LS_heatmaps/regions/both_on_S.tfSort.bed

# join in tf enrichment over noab
# sort by enrichment
join -j5 -t $'\t' <(sort -k5,5 $BOTHL) <(sort -k5,5 $BOTHS) | join -j1 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,2.12,2.15,2.18,2.21 - <(sort -k1,1 $PS_ENRICH) | awk -v OFS="\t" -v FS="\t" '{if($12 != "." || $13 != "." || $14 != "." || $15 != ".") {print $2,$3,$4,$5,$1,$6,$7,$8,$9,$10,$1,$11,"PS"} else \
{print $2,$3,$4,$5,$1,$6,$7,$8,$9,$10,$1,$11,"."}}' | cut -f1-6,13 | sort -k7,7r -k6,6nr > tf/LS_heatmaps/regions/both_on_L.tfcatSort.bed
join -j5 -t $'\t' <(sort -k5,5 $BOTHL) <(sort -k5,5 $BOTHS) | join -j1 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,2.12,2.15,2.18,2.21 - <(sort -k1,1 $PS_ENRICH) | awk -v OFS="\t" -v FS="\t" '{if($12 != "." || $13 != "." || $14 != "." || $15 != ".") {print $2,$3,$4,$5,$1,$6,$7,$8,$9,$10,$1,$11,"PS"} else \
{print $2,$3,$4,$5,$1,$6,$7,$8,$9,$10,$1,$11,"."}}' | cut -f7-13 | sort -k7,7r -k6,6nr > tf/LS_heatmaps/regions/both_on_S.tfcatSort.bed
join -j5 -t $'\t' <(sort -k5,5 $LONL) <(sort -k5,5 $LONS) | join -j1 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,2.12,2.15,2.18,2.21 - <(sort -k1,1 $PS_ENRICH) | awk -v OFS="\t" -v FS="\t" '{if($12 != "." || $13 != "." || $14 != "." || $15 != ".") {print $2,$3,$4,$5,$1,$6,$7,$8,$9,$10,$1,$11,"PS"} else \
{print $2,$3,$4,$5,$1,$6,$7,$8,$9,$10,$1,$11,"."}}' | cut -f1-6,13 | sort -k7,7r -k6,6nr > tf/LS_heatmaps/regions/L_on_L.tfcatSort.bed
join -j5 -t $'\t' <(sort -k5,5 $LONL) <(sort -k5,5 $LONS) | join -j1 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,2.12,2.15,2.18,2.21 - <(sort -k1,1 $PS_ENRICH) | awk -v OFS="\t" -v FS="\t" '{if($12 != "." || $13 != "." || $14 != "." || $15 != ".") {print $2,$3,$4,$5,$1,$6,$7,$8,$9,$10,$1,$11,"PS"} else \
{print $2,$3,$4,$5,$1,$6,$7,$8,$9,$10,$1,$11,"."}}' | cut -f7-13 | sort -k7,7r -k6,6nr > tf/LS_heatmaps/regions/L_on_S.tfcatSort.bed
join -j5 -t $'\t' <(sort -k5,5 $SONL) <(sort -k5,5 $SONS) | join -j1 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,2.12,2.15,2.18,2.21 - <(sort -k1,1 $PS_ENRICH) | awk -v OFS="\t" -v FS="\t" '{if($12 != "." || $13 != "." || $14 != "." || $15 != ".") {print $2,$3,$4,$5,$1,$6,$7,$8,$9,$10,$1,$11,"PS"} else \
{print $2,$3,$4,$5,$1,$6,$7,$8,$9,$10,$1,$11,"."}}' | cut -f1-6,13 | sort -k7,7r -k6,6nr > tf/LS_heatmaps/regions/S_on_L.tfcatSort.bed
join -j5 -t $'\t' <(sort -k5,5 $SONL) <(sort -k5,5 $SONS) | join -j1 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,2.12,2.15,2.18,2.21 - <(sort -k1,1 $PS_ENRICH) | awk -v OFS="\t" -v FS="\t" '{if($12 != "." || $13 != "." || $14 != "." || $15 != ".") {print $2,$3,$4,$5,$1,$6,$7,$8,$9,$10,$1,$11,"PS"} else \
{print $2,$3,$4,$5,$1,$6,$7,$8,$9,$10,$1,$11,"."}}' | cut -f7-13 | sort -k7,7r -k6,6nr > tf/LS_heatmaps/regions/S_on_S.tfcatSort.bed

# compute matrices with tf enrichment
LONL=tf/LS_heatmaps/regions/L_on_L.tfcatSort.bed
LONS=tf/LS_heatmaps/regions/L_on_S.tfcatSort.bed
SONS=tf/LS_heatmaps/regions/S_on_S.tfcatSort.bed
SONL=tf/LS_heatmaps/regions/S_on_L.tfcatSort.bed
BOTHL=tf/LS_heatmaps/regions/both_on_L.tfcatSort.bed
BOTHS=tf/LS_heatmaps/regions/both_on_S.tfcatSort.bed

computeMatrix reference-point --referencePoint center -a 5000 -b 5000 -R $LONL $LONS $BOTHL $BOTHS $SONL $SONS -S $POU_ENR --binSize 50 --sortRegions keep --missingDataAsZero -o tf/LS_heatmaps/matrix/LvS_pouenrich.tfcatSort.matrix -p8
computeMatrix reference-point --referencePoint center -a 5000 -b 5000 -R $LONL $LONS $BOTHL $BOTHS $SONL $SONS -S $SOX_ENR --binSize 50 --sortRegions keep --missingDataAsZero -o tf/LS_heatmaps/matrix/LvS_soxenrich.tfcatSort.matrix -p8

# plot heatmap in viridis
plotHeatmap -m tf/LS_heatmaps/matrix/LvS_pouenrich.tfcatSort.matrix -o tf/LS_heatmaps/pdfs/LvS_pouenrich.tfcatSort.pdf --sortRegions keep --dpi 1000 --heatmapWidth 2 --heatmapHeight 7 --xAxisLabel '' --regionsLabel 'lonL' 'lonS' 'bothL' 'bothS' 'sonL' 'sonS' --samplesLabel 'miler pou enrichment' --colorMap viridis --zMin 0.7 --zMax 2.2
plotHeatmap -m tf/LS_heatmaps/matrix/LvS_soxenrich.tfcatSort.matrix -o tf/LS_heatmaps/pdfs/LvS_soxenrich.tfcatSort.pdf --sortRegions keep --dpi 1000 --heatmapWidth 2 --heatmapHeight 7 --xAxisLabel '' --regionsLabel 'lonL' 'lonS' 'bothL' 'bothS' 'sonL' 'sonS' --samplesLabel 'miler sox enrichment' --colorMap viridis --zMin 0.7 --zMax 2.2

# quantify differences
mkdir -p tf/quant/regions

# calculate tf coverage
bedtools coverage -counts -a tf/LS_heatmaps/regions/both_on_L.k27Sort.bed -b $POU_COV | bedtools coverage -counts -a stdin -b $SOX_COV | cut -f1-5,7,8 | sort -k5,5 > tf/quant/regions/both_on_L.ps.bed
bedtools coverage -counts -a tf/LS_heatmaps/regions/both_on_S.k27Sort.bed -b $POU_COV | bedtools coverage -counts -a stdin -b $SOX_COV | cut -f1-5,7,8 | sort -k5,5 > tf/quant/regions/both_on_S.ps.bed
bedtools coverage -counts -a tf/LS_heatmaps/regions/L_on_L.k27Sort.bed -b $POU_COV | bedtools coverage -counts -a stdin -b $SOX_COV | cut -f1-5,7,8 | sort -k5,5 > tf/quant/regions/L_on_L.ps.bed
awk -v OFS="\t" -v FS="\t" 'midpt=int(($3+$2)/2) {print $1,midpt,midpt+1,$4,$5}' tf/LS_heatmaps/regions/L_on_S.k27Sort.bed | bedtools slop -b 250 -i stdin -g $CHRSIZES | bedtools coverage -counts -a stdin -b $POU_COV | bedtools coverage -counts -a stdin -b $SOX_COV | sort -k5,5 > tf/quant/regions/L_on_S.ps.bed
bedtools coverage -counts -a tf/LS_heatmaps/regions/S_on_S.k27Sort.bed -b $POU_COV | bedtools coverage -counts -a stdin -b $SOX_COV | cut -f1-5,7,8 | sort -k5,5 > tf/quant/regions/S_on_S.ps.bed
awk -v OFS="\t" -v FS="\t" 'midpt=int(($3+$2)/2) {print $1,midpt,midpt+1,$4,$5}' tf/LS_heatmaps/regions/S_on_L.k27Sort.bed | bedtools slop -b 250 -i stdin -g $CHRSIZES | bedtools coverage -counts -a stdin -b $POU_COV | bedtools coverage -counts -a stdin -b $SOX_COV | sort -k5,5 > tf/quant/regions/S_on_L.ps.bed

# normalize to RPM
awk -v OFS="\t" -v FS="\t" 'Ppm=1000000/11972449, Spm=1000000/18510017 {print $0,$6*Ppm,$7*Spm}' tf/quant/regions/both_on_L.ps.bed > tf/quant/regions/both_on_L.rpm.bed
awk -v OFS="\t" -v FS="\t" 'Ppm=1000000/11972449, Spm=1000000/18510017 {print $0,$6*Ppm,$7*Spm}' tf/quant/regions/both_on_S.ps.bed > tf/quant/regions/both_on_S.rpm.bed
awk -v OFS="\t" -v FS="\t" 'Ppm=1000000/11972449, Spm=1000000/18510017 {print $0,$6*Ppm,$7*Spm}' tf/quant/regions/L_on_L.ps.bed > tf/quant/regions/L_on_L.rpm.bed
awk -v OFS="\t" -v FS="\t" 'Ppm=1000000/11972449, Spm=1000000/18510017 {print $0,$6*Ppm,$7*Spm}' tf/quant/regions/L_on_S.ps.bed > tf/quant/regions/L_on_S.rpm.bed
awk -v OFS="\t" -v FS="\t" 'Ppm=1000000/11972449, Spm=1000000/18510017 {print $0,$6*Ppm,$7*Spm}' tf/quant/regions/S_on_L.ps.bed > tf/quant/regions/S_on_L.rpm.bed
awk -v OFS="\t" -v FS="\t" 'Ppm=1000000/11972449, Spm=1000000/18510017 {print $0,$6*Ppm,$7*Spm}' tf/quant/regions/S_on_S.ps.bed > tf/quant/regions/S_on_S.rpm.bed

# plot distributions
# stat LvS - pou;sox - L;LS;S
Rscript scripts_R/ps_quant.R

ENH

:<<PS_MOTIF

# call tf motifs at enhancers and test significance
# make directories
mkdir tf/homer_psenh/regions
mkdir tf/homer_psenh/fasta
mkdir tf/homer_psenh/results

# pull out regions
grep PS $LONL > tf/homer_psenh/regions/L_on_L_PS.bed
grep PS $LONS | awk -v OFS="\t" -v FS="\t" 'midpt=int(($3+$2)/2) {print $1,midpt,midpt+1,$4,$5,$6,$7}' | bedtools slop -b 250 -i stdin -g $SIZE > tf/homer_psenh/regions/L_on_S_PS.bed
grep -v PS $LONL > tf/homer_psenh/regions/L_on_L_non.bed
grep -v PS $LONS | awk -v OFS="\t" -v FS="\t" 'midpt=int(($3+$2)/2) {print $1,midpt,midpt+1,$4,$5,$6,$7}' | bedtools slop -b 250 -i stdin -g $SIZE > tf/homer_psenh/regions/L_on_S_non.bed

#
grep PS $SONS > tf/homer_psenh/regions/S_on_S_PS.bed
grep PS $SONL | awk -v OFS="\t" -v FS="\t" 'midpt=int(($3+$2)/2) {print $1,midpt,midpt+1,$4,$5,$6,$7}' | bedtools slop -b 250 -i stdin -g $SIZE > tf/homer_psenh/regions/S_on_L_PS.bed
grep -v PS $SONS > tf/homer_psenh/regions/S_on_S_non.bed
grep -v PS $SONL | awk -v OFS="\t" -v FS="\t" 'midpt=int(($3+$2)/2) {print $1,midpt,midpt+1,$4,$5,$6,$7}' | bedtools slop -b 250 -i stdin -g $SIZE > tf/homer_psenh/regions/S_on_L_non.bed

#
grep PS $BOTHL > tf/homer_psenh/regions/both_on_L_PS.bed
grep PS $BOTHS > tf/homer_psenh/regions/both_on_S_PS.bed
grep -v PS $BOTHL > tf/homer_psenh/regions/both_on_L_non.bed
grep -v PS $BOTHS > tf/homer_psenh/regions/both_on_S_non.bed

# turn beds into fasta
for i in tf/homer_psenh/regions/*.bed; do bedtools getfasta -nameOnly -fi $GENOME -fo tf/homer_psenh/fasta/`basename $i .bed`.fa -bed $i; done

# run homer
POU_MOTIF=files/data/pou_motif.motif
SOX_MOTIF=files/data/sox_motif.motif
cat $POU_MOTIF $SOX_MOTIF > tf/homer_psenh/ps_motifs.motif
for i in tf/homer_psenh/fasta/*.fa; do findMotifs.pl $i fasta tf/homer_psenh/results -find tf/homer_psenh/ps_motifs.motif -nomotif -p 6 > tf/homer_psenh/results/`basename $i .fa`.txt; done

# count instances of motifs
wc tf/homer_psenh/regions/*.bed
for i in tf/homer_psenh/results/*.txt; do grep Oct6 $i | cut -f1 | sort | uniq | wc; done
for i in tf/homer_psenh/results/*.txt; do grep Sox3 $i | cut -f1 | sort | uniq | wc; done
#
# stat tests
# pou motif
# Lon - PS = p-value = 7.823e-11
#chisq.test(matrix(c(479,1250,317,1412),2))
# Lon - non = p-value = 0.5011
#chisq.test(matrix(c(107,725,97,735),2))
# Son - PS = p-value = 8.353e-13
#chisq.test(matrix(c(383,827,229,981),2))
# Son - non = p-value = 0.004598
#chisq.test(matrix(c(84,482,52,514),2))
# both - PS = p-value = 0.6535
#chisq.test(matrix(c(206,744,197,753),2))
# both - non = p-value = 0.3926
#chisq.test(matrix(c(21,341,15,347),2))
#
# sox motif
# Lon - PS = p-value = 5.419e-10
#chisq.test(matrix(c(893,836,710,1019),2))
# Lon - non = p-value = 0.05855
#chisq.test(matrix(c(259,573,223,609),2))
# Son - PS = p-value = 1.53e-10
#chisq.test(matrix(c(479,731,637,573),2))
# Son - non = p-value = 0.948
#chisq.test(matrix(c(166,400,168,398),2))
# both - PS = p-value = 1
#chisq.test(matrix(c(406,544,407,543),2))
# both - non = p-value = 0.9269
#chisq.test(matrix(c(76,286,74,288),2))

PS_MOTIF

:<<FWSW

# heatmaps over TSS for tf enrichment
# make output directories
mkdir -p tf/fsm_heatmaps/regions
mkdir -p tf/fsm_heatmaps/matrix
mkdir -p tf/fsm_heatmaps/pdf

# join in tss coordinates and remove unk and mito
join -t $'\t' -1 4 -2 1 <(cut -f1-6 tss/TSS_select.bed | sort -k4,4) <(cut -f1,6,7,14,15 rna/gene_ps_cats.txt | sort -k1,1) -o 1.1,1.2,1.3,1.4,1.5,1.6,2.2,2.3,2.4,2.5 | grep -v Unk | grep -v chrUn | grep -v chrM > tf/fsm_heatmaps/regions/gene_ps_cats.bed

# pull out firstwave, secondwave, and maternal regions
# 1867 fw genes - use same number for inact
awk -v OFS="\t" -v FS="\t" '{if($10 == "PSstrong" || $10 == "PSweaker" || $10 == "FWnon") print $0}' tf/fsm_heatmaps/regions/gene_ps_cats.bed | sort -k9,9nr > tf/fsm_heatmaps/regions/fw_ps_cats.bed
awk -v OFS="\t" -v FS="\t" '{if($10 == "SWnon") print $0}' tf/fsm_heatmaps/regions/gene_ps_cats.bed | sort -k9,9nr > tf/fsm_heatmaps/regions/sw_ps_cats.bed
awk -v OFS="\t" -v FS="\t" '{if($10 == "NoTxMaternal" || $10 == "NoTx") print $0}' tf/fsm_heatmaps/regions/gene_ps_cats.bed | sort -k9,9nr | grep -v "e-0" | head -1867 > tf/fsm_heatmaps/regions/notx_top_ps_cats.bed

# plot tf enrichment over fw sw and notx sample
for TF in s8_pou s8_sox
do
    BW=files/bigwigs/${TF}_enrich.bw
    WD=tf/fsm_heatmaps/regions
    WD2=tf/fsm_heatmaps/matrix
    WD3=tf/fsm_heatmaps/pdf
    PREFIX=$WD2/${TF}_enrich_2500
    computeMatrix reference-point -p8 --referencePoint center -a 50000 -b 50000 -R $WD/fw_ps_cats.bed $WD/sw_ps_cats.bed $WD/notx_top_ps_cats.bed -S $BW --binSize 2500 --averageTypeBins mean --sortRegions keep --missingDataAsZero -o ${PREFIX}.matrix
    PREFIX2=$WD3/${TF}_enrich_2500
    plotProfile -m ${PREFIX}.matrix -o ${PREFIX2}.profile.pdf --plotHeight 5 --plotWidth 8
    plotHeatmap -m ${PREFIX}.matrix -o ${PREFIX2}.heatmap.pdf --dpi 1000 --heatmapWidth 4 --heatmapHeight 7 --xAxisLabel '' --regionsLabel 'fw' 'sw' 'notx' --colorMap viridis --dpi 150 --sortRegions keep --zMin 0.7 --zMax 1.2
done

# make rna-seq hm to accompany
Rscript scripts_R/fwsw_rna.R

FWSW


:<<COUNT

# plot binding site density between L and S homeologs
# make directories
mkdir -p tf/LvS/psAff_sections
mkdir -p tf/LvS/matrix
mkdir -p tf/LvS/pdf

# identify ps affected homeologs
# pull out ps affected genes from fw conf chx ps
awk -v OFS="\t" -v FS="\t" '{if(($18 == "PSstrong" || $18 == "PSweaker") || ($20 == "PSstrong" || $20 == "PSweaker")) print $0}' $CATS > tf/LvS/ps_homeo.txt

# pull out homeologs based on activation categories
awk -v OFS="\t" -v FS="\t" '{if(($18 == "PSstrong" || $18 == "PSweaker") && ($20 == "FWnon" || $20 == "SWnon" || $20 == "NoTx" || $20 == "NoTxMaternal")) print $0}' tf/LvS/ps_homeo.txt > tf/LvS/lon_ps_homeo.txt
awk -v OFS="\t" -v FS="\t" '{if(($18 == "PSstrong" || $18 == "PSweaker") && ($20 == "PSstrong" || $20 == "PSweaker")) print $0}' tf/LvS/ps_homeo.txt > tf/LvS/both_ps_homeo.txt
awk -v OFS="\t" -v FS="\t" '{if(($20 == "PSstrong" || $20 == "PSweaker") && ($18 == "FWnon" || $18 == "SWnon" || $18 == "NoTx" || $18 == "NoTxMaternal")) print $0}' tf/LvS/ps_homeo.txt > tf/LvS/son_ps_homeo.txt

# add bed intervals
# L regions
join -1 4 -2 2 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.17,2.18 <(sort -t "${TAB}" -k4,4 $LPROM) <(sort -t "${TAB}" -k2,2 tf/LvS/lon_ps_homeo.txt) | sort -t "${TAB}" -k8,8nr > tf/LvS/psAff_sections/LpromL_psAff.txt
join -1 4 -2 2 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.17,2.18 <(sort -t "${TAB}" -k4,4 $LPROM) <(sort -t "${TAB}" -k2,2 tf/LvS/both_ps_homeo.txt) | sort -t "${TAB}" -k8,8nr > tf/LvS/psAff_sections/LpromLS_psAff.txt
join -1 4 -2 2 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.17,2.18 <(sort -t "${TAB}" -k4,4 $LPROM) <(sort -t "${TAB}" -k2,2 tf/LvS/son_ps_homeo.txt) | sort -t "${TAB}" -k8,8nr > tf/LvS/psAff_sections/LpromS_psAff.txt

# S regions
join -1 4 -2 1 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.19,2.20 <(sort -t "${TAB}" -k4,4 $SPROM) <(sort -t "${TAB}" -k1,1 tf/LvS/lon_ps_homeo.txt) | sort -t "${TAB}" -k8,8nr > tf/LvS/psAff_sections/SpromL_psAff.txt
join -1 4 -2 1 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.19,2.20 <(sort -t "${TAB}" -k4,4 $SPROM) <(sort -t "${TAB}" -k1,1 tf/LvS/both_ps_homeo.txt) | sort -t "${TAB}" -k8,8nr > tf/LvS/psAff_sections/SpromLS_psAff.txt
join -1 4 -2 1 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.19,2.20 <(sort -t "${TAB}" -k4,4 $SPROM) <(sort -t "${TAB}" -k1,1 tf/LvS/son_ps_homeo.txt) | sort -t "${TAB}" -k8,8nr > tf/LvS/psAff_sections/SpromS_psAff.txt

# make matrices
LONL=tf/LvS/psAff_sections/LpromL_psAff.txt
LONS=tf/LvS/psAff_sections/LpromS_psAff.txt
SONS=tf/LvS/psAff_sections/SpromS_psAff.txt
SONL=tf/LvS/psAff_sections/SpromL_psAff.txt
BOTHL=tf/LvS/psAff_sections/LpromLS_psAff.txt
BOTHS=tf/LvS/psAff_sections/SpromLS_psAff.txt

computeMatrix reference-point --referencePoint center -a 25000 -b 25000 -R $LONL $BOTHL $LONS $SONL $BOTHS $SONS -S files/bigwigs/pou_diff.bw files/bigwigs/sox_diff.bw files/bigwigs/pou_sox_diff_union.bw --averageTypeBins max --binSize 500 --sortRegions keep --missingDataAsZero -o tf/LvS/matrix/ps_diff.matrix -p2
computeMatrix reference-point --referencePoint center -a 25000 -b 25000 -R $LONL $BOTHL $LONS $SONL $BOTHS $SONS -S files/bigwigs/pou_cons.bw files/bigwigs/sox_cons.bw files/bigwigs/pou_sox_cons_union.bw --averageTypeBins max --binSize 500 --sortRegions keep --missingDataAsZero -o tf/LvS/matrix/ps_cons.matrix -p2

# make plotheatmaps
plotHeatmap -m  tf/LvS/matrix/ps_diff.matrix -o tf/LvS/pdf/ps_diff.pdf --sortRegions keep --dpi 1000 --heatmapWidth 21 --heatmapHeight 21 --xAxisLabel '' --regionsLabel 'LpromL' 'LpromLS' 'LpromS' 'SpromL' 'SpromLS' 'SpromS' --samplesLabel 'pou diff' 'sox diff' 'tf diff' --colorList white,black --zMax 1 --zMin 0 --whatToShow 'heatmap and colorbar'
plotHeatmap -m  tf/LvS/matrix/ps_cons.matrix -o tf/LvS/pdf/ps_cons.pdf --sortRegions keep --dpi 1000 --heatmapWidth 21 --heatmapHeight 21 --xAxisLabel '' --regionsLabel 'LpromL' 'LpromLS' 'LpromS' 'SpromL' 'SpromLS' 'SpromS' --samplesLabel 'pou cons' 'sox cons' 'tf cons' --colorList white,black --zMax 1 --zMin 0 --whatToShow 'heatmap and colorbar'

COUNT

exit
