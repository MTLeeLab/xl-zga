# figure 5 enhancer cmds
# lift lae regions to trop and then trop regions to zf

S9PEAKS=hom/xt_chip/GSM1350515_Tropicalis_H3K27Ac_stage9_broad_peaks.homer.txt
S8PEAKS=hom/xt_chip/GSM1350514_Tropicalis_H3K27Ac_stage8_broad_peaks.homer.txt
TRO27=files/chains/xenTro2ToXenTro7.over.chain
TRO79=files/chains/xenTro7ToXenTro9.over.chain
LAE2TRO9=files/chains/xenLae2ToXenTro9.over.chain
GENOME=files/anno/genome.fa
CHRSIZES=files/anno/chr_sizes.txt
TRO9TRO7=files/chains/xenTro9ToXenTro7.over.chain.gz
TRO7RER10=files/chains/xenTro7ToDanRer10.over.chain.gz
RER10RER11=files/chains/danRer10ToDanRer11.over.chain.gz
GENOME2=files/anno/genome_trop.fa

:<<"TROP"

# lift trop v2 k27ac peaks to v9
# v2 to v7
mkdir -p hom/xt_chip/lift_9
mkdir -p hom/xt_chip/lift_8
# !!! add trop peak files from Gupta et al into xt_chip!!!

cut -f1-6 $S9PEAKS | awk -v OFS="\t" -v FS="\t" '{print $2,$3,$4,$1,$6,$5}' | tail -11314 | liftOver -minMatch=0.9 stdin $TRO27 hom/xt_chip/lift_9/trop_H3K27Ac_stage9_lift_v7.bed hom/xt_chip/lift_9/unmapped_st9_v7
cut -f1-6 $S8PEAKS | awk -v OFS="\t" -v FS="\t" '{print $2,$3,$4,$1,$6,$5}' | tail -8078 | liftOver -minMatch=0.9 stdin $TRO27 hom/xt_chip/lift_8/trop_H3K27Ac_stage8_lift_v7.bed hom/xt_chip/lift_8/unmapped_st8_v7

# v7 to v9
liftOver -minMatch=0.9 hom/xt_chip/lift_9/trop_H3K27Ac_stage9_lift_v7.bed $TRO79 hom/xt_chip/lift_9/trop_H3K27Ac_stage9_lift_v9.bed hom/xt_chip/lift_9/unmapped_st9_v9
liftOver -minMatch=0.9 hom/xt_chip/lift_8/trop_H3K27Ac_stage8_lift_v7.bed $TRO79 hom/xt_chip/lift_8/trop_H3K27Ac_stage8_lift_v9.bed hom/xt_chip/lift_8/unmapped_st8_v9

TROP

:<<LAETRO

# pull out L (N=18137) and S (N=16259) regions
mkdir -p hom/anno
ln -s lift/combined/all_enh_homeologs.txt hom/anno/enh_tf_anno.txt
grep -e 'LS-' -e 'L-' -e 'S-NA' -e 'l-' hom/anno/enh_tf_anno.txt | awk -v OFS="\t" -v FS="\t" '{print $2,$3,$4,$1}' > hom/anno/Lenh_anno.bed
grep -e 'LS-' -e 'L-' -e 'S-NA' -e 's-' hom/anno/enh_tf_anno.txt | awk -v OFS="\t" -v FS="\t" '{print $6,$7,$8,$1}' > hom/anno/Senh_anno.bed

# lift laevis enhancers to trop (N=19233 initially)
# make directories
mkdir -p hom/lift
mkdir -p hom/lift/500bp
mkdir -p hom/lift/rescue

# liftover 500bp for L (N=15305) and S (N=13685)
liftOver -minMatch=0.1 hom/anno/Lenh_anno.bed $LAE2TRO9 hom/lift/500bp/Lenh.bed hom/lift/500bp/unmapped_Lenh
liftOver -minMatch=0.1 hom/anno/Senh_anno.bed $LAE2TRO9 hom/lift/500bp/Senh.bed hom/lift/500bp/unmapped_Senh

#
### PHASE 2 LIFTOVER
# Rescue the remaining enhancers (N=1800) that didn't liftover with 500bp liftover using upstream and downstream 2500bp liftover
# remove any enhancers already accounted for above due to asymmetric liftover among the both on (L to S but not S to L, e.g.)
mkdir -p hom/rescue

# pull out unlifted regions for L (N = 2832) and S (N=2574)
grep -v "#" hom/lift/500bp/unmapped_Lenh > hom/rescue/unmappedL_filtered.bed
grep -v "#" hom/lift/500bp/unmapped_Senh > hom/rescue/unmappedS_filtered.bed

## Liftover full 5kb region around L enh (N=1754 lift) and S enh (N=1645 lift)
bedtools slop -b 2250 -g $CHRSIZES -i hom/rescue/unmappedL_filtered.bed > hom/rescue/unmappedL_filtered_slop5k.bed
bedtools slop -b 2250 -g $CHRSIZES -i hom/rescue/unmappedS_filtered.bed > hom/rescue/unmappedS_filtered_slop5k.bed
liftOver -minMatch=0.1 hom/rescue/unmappedL_filtered_slop5k.bed $LAE2TRO9 hom/rescue/enhL_ext_lift.bed hom/rescue/unmappedL.txt
liftOver -minMatch=0.1 hom/rescue/unmappedS_filtered_slop5k.bed $LAE2TRO9 hom/rescue/enhS_ext_lift.bed hom/rescue/unmappedS.txt

## Next, liftover left and right flanks separately, then join together
## for L; left, N =1624; right, N = 1613 lifted
## for S; left, N =1502; right, N = 1551 lifted
# make arms
bedtools flank -l 2250 -r 0 -g $CHRSIZES -i hom/rescue/unmappedL_filtered.bed > hom/rescue/unmappedL_leftarm.bed
bedtools flank -l 2250 -r 0 -g $CHRSIZES -i hom/rescue/unmappedS_filtered.bed > hom/rescue/unmappedS_leftarm.bed
bedtools flank -l 0 -r 2250 -g $CHRSIZES -i hom/rescue/unmappedL_filtered.bed > hom/rescue/unmappedL_rightarm.bed
bedtools flank -l 0 -r 2250 -g $CHRSIZES -i hom/rescue/unmappedS_filtered.bed > hom/rescue/unmappedS_rightarm.bed

# lift over arms
liftOver -minMatch=0.1 hom/rescue/unmappedL_leftarm.bed $LAE2TRO9 hom/rescue/enhL_leftarm_lift.bed hom/rescue/unmappedL_leftarm.txt
liftOver -minMatch=0.1 hom/rescue/unmappedS_leftarm.bed $LAE2TRO9 hom/rescue/enhS_leftarm_lift.bed hom/rescue/unmappedS_leftarm.txt
liftOver -minMatch=0.1 hom/rescue/unmappedL_rightarm.bed $LAE2TRO9 hom/rescue/enhL_rightarm_lift.bed hom/rescue/unmappedL_rightarm.txt
liftOver -minMatch=0.1 hom/rescue/unmappedS_rightarm.bed $LAE2TRO9 hom/rescue/enhS_rightarm_lift.bed hom/rescue/unmappedS_rightarm.txt

# Pair the arms together and find the midpoint for L (N = 992 successfully joined) and S (N = 973 successfully joined)
sort -k4,4 hom/rescue/enhL_leftarm_lift.bed | join -j4 -t $'\t' - <(sort -k4,4 hom/rescue/enhL_rightarm_lift.bed) | awk 'OFS="\t", FS="\t" {
if ($3>=$6 && $3<=$7 && $4>=$6 && $4<=$7) {print $2,$3,$4,$1} \
else if ($6>=$3 && $6<=$4 && $7>=$3 && $7<=$4) {print $2,$6,$7,$1} \
else if ($4<=$6 && $4<=$7) {print $2,$4,$6,$1} \
else if ($7<=$3 && $7<=$4) {print $2,$7,$3,$1} \
else if ($3>$6 && $3 <$7) {print $2,$3,$7,$1} \
else if ($4>$6 && $4<$7) {print $2,$6,$4,$1} \
}' | sort -k4,4 > hom/rescue/enhL_left_right_merge.bed
sort -k4,4 hom/rescue/enhS_leftarm_lift.bed | join -j4 -t $'\t' - <(sort -k4,4 hom/rescue/enhS_rightarm_lift.bed) | awk 'OFS="\t", FS="\t" {
if ($3>=$6 && $3<=$7 && $4>=$6 && $4<=$7) {print $2,$3,$4,$1} \
else if ($6>=$3 && $6<=$4 && $7>=$3 && $7<=$4) {print $2,$6,$7,$1} \
else if ($4<=$6 && $4<=$7) {print $2,$4,$6,$1} \
else if ($7<=$3 && $7<=$4) {print $2,$7,$3,$1} \
else if ($3>$6 && $3 <$7) {print $2,$3,$7,$1} \
else if ($4>$6 && $4<$7) {print $2,$6,$4,$1} \
}' | sort -k4,4 > hom/rescue/enhS_left_right_merge.bed

# Join the full 5k liftover with merged arms -- any successful left/right liftover will need to have a gap in the middle less than the size of the full 5k liftover (678 L regions and 682 S regions)
sort -k4,4 hom/rescue/enhL_ext_lift.bed | join -j4 -t $'\t' hom/rescue/enhL_left_right_merge.bed - | awk 'OFS="\t" {if ($7-$6 >= $4-$3) print $2,$3,$4,$1}' > hom/rescue/enhL_combined_lift.bed
sort -k4,4 hom/rescue/enhS_ext_lift.bed | join -j4 -t $'\t' hom/rescue/enhS_left_right_merge.bed - | awk 'OFS="\t" {if ($7-$6 >= $4-$3) print $2,$3,$4,$1}' > hom/rescue/enhS_combined_lift.bed

# Join everything together with the original 500bp liftovers
# Each pair is given an ID derived from the individual enhancer IDs
# 15983 unique L regions and 14367 unique S regions
cat hom/lift/500bp/Lenh.bed hom/rescue/enhL_combined_lift.bed | sort -k4,4 | uniq > hom/rescue/enhL_combined_full.bed
cat hom/lift/500bp/Senh.bed hom/rescue/enhS_combined_lift.bed | sort -k4,4 | uniq > hom/rescue/enhS_combined_full.bed

# make paired table with trop acetylation status
# intersect with xt acetylation
mkdir -p hom/enh
bedtools intersect -wao -a hom/rescue/enhL_combined_full.bed -b hom/xt_chip/lift_9/trop_H3K27Ac_stage9_lift_v9.bed | cut -f1-5 | awk -v OFS="\t" -v FS="\t" '{if($5 != ".") {print$1,$2,$3,$4,"X"} \
else {print $0}}' | sort -k4,4 | uniq > hom/enh/Lenh_xts9_anno.bed
bedtools intersect -wao -a hom/rescue/enhS_combined_full.bed -b hom/xt_chip/lift_9/trop_H3K27Ac_stage9_lift_v9.bed | cut -f1-5 | awk -v OFS="\t" -v FS="\t" '{if($5 != ".") {print$1,$2,$3,$4,"X"} \
else {print $0}}' | sort -k4,4 | uniq > hom/enh/Senh_xts9_anno.bed
# merge trop regions into all enhancer table
sort -k1,1 hom/anno/enh_tf_anno.txt | join -1 1 -2 4 -a 1 -t $'\t' -o auto -e NA - hom/enh/Lenh_xts9_anno.bed | join -1 1 -2 4 -a 1 -t $'\t' -o auto -e NA - hom/enh/Senh_xts9_anno.bed > hom/xl_xt_enh_anno.txt

LAETRO

:<<TRORER

# liftover full trop regions from L (N=15983) and S (N=14367) to zebrafish
# make working directories
mkdir -p hom/lift_zf
mkdir -p hom/lift_zf/xt7
mkdir -p hom/lift_zf/dr10
mkdir -p hom/lift_zf/dr11
mkdir -p hom/rescue_zf

# lift xentro9 regions from before to xentro7 for L (N=15910 lifted) and S (N=14288 lifted)
liftOver -minMatch=0.9 hom/rescue/enhL_combined_full.bed $TRO9TRO7 hom/lift_zf/xt7/xt7_lae_enhL.bed hom/lift_zf/xt7/unmappedL_xt7
liftOver -minMatch=0.9 hom/rescue/enhS_combined_full.bed $TRO9TRO7 hom/lift_zf/xt7/xt7_lae_enhS.bed hom/lift_zf/xt7/unmappedS_xt7

# lift xentro7 to danrer10 for L (N=3318) and S (N=3055)
liftOver -minMatch=0.1 hom/lift_zf/xt7/xt7_lae_enhL.bed $TRO7RER10 hom/lift_zf/dr10/dr10_lae_enhL.bed hom/lift_zf/dr10/unmappedL_dr10
liftOver -minMatch=0.1 hom/lift_zf/xt7/xt7_lae_enhS.bed $TRO7RER10 hom/lift_zf/dr10/dr10_lae_enhS.bed hom/lift_zf/dr10/unmappedS_dr10

#
### PHASE 2 LIFTOVER
# Rescue the remaining enhancers (N=1800) that didn't liftover with 500bp liftover using upstream and downstream 2500bp liftover
# remove any enhancers already accounted for above due to asymmetric liftover among the both on (L to S but not S to L, e.g.)
# make xt7 size file
samtools faidx $GENOME2
cut -f1,2 files/anno/genome_trop.fa.fai | sort -k1,1 > files/anno/xt7_chr_size.txt
CHRSIZES2=files/anno/xt7_chr_size.txt

# pull out unlifted regions for L (N = 12592) and S (N=11233)
grep -v "#" hom/lift_zf/dr10/unmappedL_dr10 > hom/rescue_zf/unmappedL_filtered.bed
grep -v "#" hom/lift_zf/dr10/unmappedS_dr10 > hom/rescue_zf/unmappedS_filtered.bed

## Liftover full 5kb region around L enh (N=2236 lift) and S enh (N=1989 lift)
bedtools slop -b 2250 -g $CHRSIZES2 -i hom/rescue_zf/unmappedL_filtered.bed > hom/rescue_zf/unmappedL_filtered_slop5k.bed
bedtools slop -b 2250 -g $CHRSIZES2 -i hom/rescue_zf/unmappedS_filtered.bed > hom/rescue_zf/unmappedS_filtered_slop5k.bed
liftOver -minMatch=0.1 hom/rescue_zf/unmappedL_filtered_slop5k.bed $TRO7RER10 hom/rescue_zf/enhL_ext_lift.bed hom/rescue_zf/unmappedL.txt
liftOver -minMatch=0.1 hom/rescue_zf/unmappedS_filtered_slop5k.bed $TRO7RER10 hom/rescue_zf/enhS_ext_lift.bed hom/rescue_zf/unmappedS.txt

## Next, liftover left and right flanks separately, then join together
## for L; left, N =2123; right, N = 2116 lifted
## for S; left, N =1864; right, N = 1911 lifted
bedtools flank -l 2250 -r 0 -g $CHRSIZES2 -i hom/rescue_zf/unmappedL_filtered.bed > hom/rescue_zf/unmappedL_leftarm.bed
bedtools flank -l 2250 -r 0 -g $CHRSIZES2 -i hom/rescue_zf/unmappedS_filtered.bed > hom/rescue_zf/unmappedS_leftarm.bed
bedtools flank -l 0 -r 2250 -g $CHRSIZES2 -i hom/rescue_zf/unmappedL_filtered.bed > hom/rescue_zf/unmappedL_rightarm.bed
bedtools flank -l 0 -r 2250 -g $CHRSIZES2 -i hom/rescue_zf/unmappedS_filtered.bed > hom/rescue_zf/unmappedS_rightarm.bed

# lift over arms from xt to zf
liftOver -minMatch=0.1 hom/rescue_zf/unmappedL_leftarm.bed $TRO7RER10 hom/rescue_zf/enhL_leftarm_lift.bed hom/rescue_zf/unmappedL_leftarm.txt
liftOver -minMatch=0.1 hom/rescue_zf/unmappedS_leftarm.bed $TRO7RER10 hom/rescue_zf/enhS_leftarm_lift.bed hom/rescue_zf/unmappedS_leftarm.txt
liftOver -minMatch=0.1 hom/rescue_zf/unmappedL_rightarm.bed $TRO7RER10 hom/rescue_zf/enhL_rightarm_lift.bed hom/rescue_zf/unmappedL_rightarm.txt
liftOver -minMatch=0.1 hom/rescue_zf/unmappedS_rightarm.bed $TRO7RER10 hom/rescue_zf/enhS_rightarm_lift.bed hom/rescue_zf/unmappedS_rightarm.txt

# Pair the arms together and find the midpoint for L (N = 534 successfully joined) and S (N = 434 successfully joined)
sort -k4,4 hom/rescue_zf/enhL_leftarm_lift.bed | join -j4 -t $'\t' - <(sort -k4,4 hom/rescue_zf/enhL_rightarm_lift.bed) | awk 'OFS="\t", FS="\t" {
if ($3>=$6 && $3<=$7 && $4>=$6 && $4<=$7) {print $2,$3,$4,$1} \
else if ($6>=$3 && $6<=$4 && $7>=$3 && $7<=$4) {print $2,$6,$7,$1} \
else if ($4<=$6 && $4<=$7) {print $2,$4,$6,$1} \
else if ($7<=$3 && $7<=$4) {print $2,$7,$3,$1} \
else if ($3>$6 && $3 <$7) {print $2,$3,$7,$1} \
else if ($4>$6 && $4<$7) {print $2,$6,$4,$1} \
}' | sort -k4,4 > hom/rescue_zf/enhL_left_right_merge.bed
sort -k4,4 hom/rescue_zf/enhS_leftarm_lift.bed | join -j4 -t $'\t' - <(sort -k4,4 hom/rescue_zf/enhS_rightarm_lift.bed) | awk 'OFS="\t", FS="\t" {
if ($3>=$6 && $3<=$7 && $4>=$6 && $4<=$7) {print $2,$3,$4,$1} \
else if ($6>=$3 && $6<=$4 && $7>=$3 && $7<=$4) {print $2,$6,$7,$1} \
else if ($4<=$6 && $4<=$7) {print $2,$4,$6,$1} \
else if ($7<=$3 && $7<=$4) {print $2,$7,$3,$1} \
else if ($3>$6 && $3 <$7) {print $2,$3,$7,$1} \
else if ($4>$6 && $4<$7) {print $2,$6,$4,$1} \
}' | sort -k4,4 > hom/rescue_zf/enhS_left_right_merge.bed

# Join the full 5k liftover with merged arms -- any successful left/right liftover will need to have a gap in the middle less than the size of the full 5k liftover (183 L regions and 148 S regions)
sort -k4,4 hom/rescue_zf/enhL_ext_lift.bed | join -j4 -t $'\t' hom/rescue_zf/enhL_left_right_merge.bed - | awk 'OFS="\t" {if ($7-$6 >= $4-$3) print $2,$3,$4,$1}' > hom/rescue_zf/enhL_combined_lift.bed
sort -k4,4 hom/rescue_zf/enhS_ext_lift.bed | join -j4 -t $'\t' hom/rescue_zf/enhS_left_right_merge.bed - | awk 'OFS="\t" {if ($7-$6 >= $4-$3) print $2,$3,$4,$1}' > hom/rescue_zf/enhS_combined_lift.bed

# Join everything together with the original 500bp liftovers
# Each pair is given an ID derived from the individual enhancer IDs
# 3501 unique L regions and 3203 unique S regions
cat hom/lift_zf/dr10/dr10_lae_enhL.bed hom/rescue_zf/enhL_combined_lift.bed | sort -k4,4 | uniq > hom/rescue_zf/enhL_combined_full.bed
cat hom/lift_zf/dr10/dr10_lae_enhS.bed hom/rescue_zf/enhS_combined_lift.bed | sort -k4,4 | uniq > hom/rescue_zf/enhS_combined_full.bed

# lift danRer10 regions to danRer11
# lifted L (N=3488) and S (N=3196) regions
liftOver -minMatch=0.9 hom/rescue_zf/enhL_combined_full.bed $RER10RER11 hom/lift_zf/dr11/dr11_lae_enhL.bed hom/lift_zf/dr11/unmappedL_dr11
liftOver -minMatch=0.9 hom/rescue_zf/enhS_combined_full.bed $RER10RER11 hom/lift_zf/dr11/dr11_lae_enhS.bed hom/lift_zf/dr11/unmappedS_dr11

TRORER

:<<TABLE

# call peaks on bogdanovic k27ac dome data for zebrafish
mkdir -p hom/zf_chip
bamToBed -i files/bams/zf_dome_k27ac.bam > hom/zf_chip/zf_chip_reads.bed
macs2 callpeak -t hom/zf_chip/zf_chip_reads.bed -f BED -g 4.59e8 --name dome_k27ac --outdir hom/zf_chip --verbose 3

# intersect with zf acetylation
# acetylated L (N=1291) and S (N=1180) regions
bedtools intersect -wao -a hom/lift_zf/dr11/dr11_lae_enhL.bed -b hom/zf_chip/dome_k27ac_peaks.narrowPeak | cut -f1-5 | awk -v OFS="\t" -v FS="\t" '{if($5 != ".") {print$1,$2,$3,$4,"X"} \
else {print $0}}' | sort -k4,4 | uniq > hom/enh/Lenh_zfdome_anno.bed
bedtools intersect -wao -a hom/lift_zf/dr11/dr11_lae_enhS.bed -b hom/zf_chip/dome_k27ac_peaks.narrowPeak | cut -f1-5 | awk -v OFS="\t" -v FS="\t" '{if($5 != ".") {print$1,$2,$3,$4,"X"} \
else {print $0}}' | sort -k4,4 | uniq > hom/enh/Senh_zfdome_anno.bed

# merge trop regions into all enhancer table
sort -k1,1 hom/xl_xt_enh_anno.txt | join -1 1 -2 4 -a 1 -t $'\t' -o auto -e NA - hom/enh/Lenh_zfdome_anno.bed | join -1 1 -2 4 -a 1 -t $'\t' -o auto -e NA - hom/enh/Senh_zfdome_anno.bed > hom/xl_xt_zf_enh_anno.txt

# reconcile trop regions into one
awk -v OFS="\t" -v FS="\t" '{if($2 ~ $10 && $6 == "NA") {print $0,$10":"$11"-"$12","$13} else \
if($2 !~ $10 && $6 == "NA") {print $0,"NA,NA"} else \
if($6 ~ $14 && $2 == "NA") {print $0,$14":"$15"-"$16","$17} else \
if($6 !~ $14 && $2 == "NA") {print $0,"NA,NA"} else \
if($2 ~ $10 && $2 !~ $14) {print $0,$10":"$11"-"$12","$13} else \
if($2 ~ $14 && $2 !~ $10) {print $0,$14":"$15"-"$16","$17} else \
if(($2 == "chr9_10L" || $6 == "chr9_10S") && $10 == "chr10" && $14 != "chr10") {print $0,$10":"$11"-"$12","$13} else \
if(($2 == "chr9_10L" || $6 == "chr9_10S") && $14 == "chr10" && $10 != "chr10") {print $0,$14":"$15"-"$16","$17} else \
if(($2 == "chr9_10L" || $6 == "chr9_10S") && $10 == "chr10" && $14 == "chr10" && $11 <= $15 && $16 >= $12) {print $0,$10":"$11"-"$16","$13$17} else \
if(($2 == "chr9_10L" || $6 == "chr9_10S") && $10 == "chr10" && $14 == "chr10" && $15 <= $11 && $12 >= $16) {print $0,$10":"$15"-"$12","$13$17} else \
if(($2 == "chr9_10L" || $6 == "chr9_10S") && $10 == "chr10" && $14 == "chr10" && $12 < $15 && $15-$12 < 500) {print $0,$10":"$11"-"$16","$13$17} else \
if(($2 == "chr9_10L" || $6 == "chr9_10S") && $10 == "chr10" && $14 == "chr10" && $16 < $11 && $11-$16 < 500) {print $0,$10":"$15"-"$12","$13$17} else \
if(($2 == "chr9_10L" || $6 == "chr9_10S") && $10 == "chr10" && $14 == "chr10" && $11 < $15 && $12 > $16) {print $0,$10":"$11"-"$12","$13$17} else \
if(($2 == "chr9_10L" || $6 == "chr9_10S") && $10 == "chr10" && $14 == "chr10" && $15 < $11 && $16 > $12) {print $0,$10":"$15"-"$16","$13$17} else \
if($2 !~ $10 && $2 !~ $14) {print $0,"NA,NA"} else \
if($11 <= $15 && $16 >= $12) {print $0,$10":"$11"-"$16","$13$17} else \
if($15 <= $11 && $12 >= $16) {print $0,$10":"$15"-"$12","$13$17} else \
if($12 < $15 && $15-$12 <= 500) {print $0,$10":"$11"-"$16","$13$17} else \
if($16 < $11 && $11-$16 <= 500) {print $0,$10":"$15"-"$12","$13$17} else \
if($12 < $15 && $15-$12 > 500) {print $0,"unk,unk"} else \
if($16 < $11 && $11-$16 > 500) {print $0,"unk,unk"} else \
if($11 < $15 && $12 > $16) {print $0,$10":"$11"-"$12","$13$17} else \
if($15 < $11 && $16 > $12) {print $0,$10":"$15"-"$16","$13$17} else \
{print $0,"."}}' hom/xl_xt_zf_enh_anno.txt > hom/xl_xt_zf_enh_full.txt_tmp

# replace various combined LS annotation into a single trop annotation
cut -f1-9,18- hom/xl_xt_zf_enh_full.txt_tmp | awk -v OFS="\t" -v FS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$18,$10,$11,$12,$13,$14,$15,$16,$17}' | sed s/','/'\t'/g | sed s/'\.\.'/'\.'/g | sed s/'NANA'/'NA'/g | sed -e s/'\.X'/'X'/g -e s/'XX'/'X'/g -e s/'X\.'/'X'/g  > hom/xl_xt_zf_enh_full.txt_tmp2

# reconcile zebrafish regions into one
awk -v OFS="\t" -v FS="\t" '{if($12 == "NA" && $16 == "NA") {print $0,"NA,NA"} else \
if($12 != "NA" && $16 == "NA") {print $0,$12":"$13"-"$14","$15} else \
if($16 != "NA" && $12 == "NA") {print $0,$16":"$17"-"$18","$19} else \
if($12 != $16) {print $0,"NA,NA"} else \
if($13 <= $17 && $18 >= $14) {print $0,$12":"$13"-"$18","$15$19} else \
if($17 <= $13 && $14 >= $18) {print $0,$12":"$17"-"$14","$15$19} else \
if($14 < $17 && $17-$14 <= 500) {print $0,$12":"$13"-"$18","$15$19} else \
if($18 < $13 && $13-$18 <= 500) {print $0,$12":"$17"-"$14","$15$19} else \
if($14 < $17 && $17-$14 > 500) {print $0,"unk,unk"} else \
if($18 < $13 && $13-$18 > 500) {print $0,"unk,unk"} else \
if($13 < $17 && $14 > $18) {print $0,$12":"$13"-"$14","$15$19} else \
if($17 < $13 && $18 > $14) {print $0,$12":"$17"-"$18","$15$19} else \
{print $0,"."}}' hom/xl_xt_zf_enh_full.txt_tmp2 > hom/xl_xt_zf_enh_full.txt_tmp3

# replace various combined LS annotation into a single zf annotation
cut -f1-11,20 hom/xl_xt_zf_enh_full.txt_tmp3 | sed s/','/'\t'/g | sed s/'\.\.'/'\.'/g | sed -e s/'\.X'/'X'/g -e s/'XX'/'X'/g -e s/'X\.'/'X'/g > hom/xl_xt_zf_enh_full.txt

# remove tmps
rm *.txt_tmp*

#
# make barcharts
Rscript cons_acetyl.R

TABLE

exit
