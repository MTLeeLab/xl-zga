GENOME=files/anno/genome.fa
CHRSIZES=files/anno/chr_sizes.txt
CHAIN=files/chains/xlL_xlS_med.liftOver.chain


### PHASE 1 LIFTOVER
# L/S Liftover command 1: +/- 250 bp on enhancer atac peak center

:<<"LIFT500"

#Combine proximal and distal elements
mkdir -p lift/500bp
mkdir -p lift/5kb
mkdir -p lift/combined

cut -f1-4 ../reg_thresh/distal_k27high.bed ../reg_thresh/proximal_k27high.bed | grep -v chrUn | sort -k4,4 > lift/enh.bed
liftOver -minMatch=0.1 enh.bed $CHAIN lift/500bp/enh_lift.bed lift/500bp/unmapped.txt

# For the S->L liftover, match with original L enhancers
# Redundant because matching enhancers can be many-many
# Output L BED4 columns followed by S BED4 columns

grep L lift/500bp/enh_lift.bed | bedtools intersect -wo -b - -a lift/enh.bed | cut -f1-4,8 | sort -k5,5 | join -1 5 -2 4 -o 1.1,1.2,1.3,1.4,2.1,2.2,2.3,2.4 -t $'\t' - lift/enh.bed > lift/500bp/both_on_redundant.bed_temp

# Repeat for L->S liftover, concat, unique
# Also output L BED4 columns followed by S BED4 columns (so, join is in opposite order)

grep S lift/500bp/enh_lift.bed | bedtools intersect -wo -b - -a lift/enh.bed | cut -f1-4,8 | sort -k5,5 | join -1 4 -2 5 -o 1.1,1.2,1.3,1.4,2.1,2.2,2.3,2.4 -t $'\t' lift/enh.bed - >> lift/500bp/both_on_redundant.bed_temp

# Union = 734 *redundant* regions -- the same L or S region may appear more than once
# due to many-to-many liftovers. Will need to reconcile in further analyses

sort lift/500bp/both_on_redundant.bed_temp | uniq > lift/500bp/both_on_redundant.bed
cat <(cut -f4 lift/500bp/both_on_redundant.bed) <(cut -f8 lift/500bp/both_on_redundant.bed) | sort | uniq > lift/500bp/both_on_ids.txt

# Also make a file that has the actual lifted over both-on regions
grep -f lift/500bp/both_on_ids.txt lift/500bp/enh_lift.bed > lift/500bp/both_on_lifted_regions.bed

# Single on -- filter out anything that appears in the both on list
# S on L off regions
grep L lift/500bp/enh_lift.bed | join -t $'\t' -j4 -o 1.1,1.2,1.3,1.4,2.1,2.2,2.3 lift/enh.bed - | grep -w -v -f lift/500bp/both_on_ids.txt > lift/500bp/S_on.bed

# L on S off regions
grep S lift/500bp/enh_lift.bed | join -t $'\t' -j4 -o 1.1,1.2,1.3,1.4,2.1,2.2,2.3 lift/enh.bed - | grep -w -v -f lift/500bp/both_on_ids.txt > lift/500bp/L_on.bed

LIFT500


### PHASE 2 LIFTOVER
# Rescue the remaining enhancers (N=2780) that didn't liftover with 500bp liftover using upstream and downstream 2500bp liftover; remove any enhancers already accounted for above due to asymmetric liftover among the both on (L to S but not S to L, e.g.)

:<<"LIFT5K"

grep -v "#" lift/500bp/unmapped.txt | grep -v -w -f lift/500bp/both_on_ids.txt > lift/500bp/unmapped_filtered.bed

## Liftover full 5kb region around enh
bedtools slop -b 2250 -g $CHRSIZES -i lift/500bp/unmapped_filtered.bed > lift/500bp/unmapped_filtered_slop5k.bed
liftOver -minMatch=0.1 lift/500bp/unmapped_filtered_slop5k.bed $CHAIN lift/5kb/enh_ext_lift.bed lift/5kb/unmapped.txt

## Next, liftover left and right flanks separately, then join together
bedtools flank -l 2250 -r 0 -g $CHRSIZES -i lift/500bp/unmapped_filtered.bed > lift/500bp/unmapped_leftarm.bed

bedtools flank -l 0 -r 2250 -g $CHRSIZES -i lift/500bp/unmapped_filtered.bed > lift/500bp/unmapped_rightarm.bed


liftOver -minMatch=0.1 lift/500bp/unmapped_leftarm.bed $CHAIN lift/5kb/enh_leftarm_lift.bed lift/5kb/unmapped_leftarm.txt

liftOver -minMatch=0.1 lift/500bp/unmapped_rightarm.bed $CHAIN lift/5kb/enh_rightarm_lift.bed lift/5kb/unmapped_rightarm.txt

# Pair the arms together and find the midpoint (N = 633 successfully joined)
sort -k4,4 lift/5kb/enh_leftarm_lift.bed | join -j4 -t $'\t' - <(sort -k4,4 lift/5kb/enh_rightarm_lift.bed) | awk 'OFS="\t" {
if ($3>=$6 && $3<=$7 && $4>=$6 && $4<=$7) {print $2,$3,$4,$1} \
else if ($6>=$3 && $6<=$4 && $7>=$3 && $7<=$4) {print $2,$6,$7,$1} \
else if($4<=$6 && $4<=$7) {print $2,$4,$6,$1} \
else if ($7<=$3 && $7<=$4) {print $2,$7,$3,$1} \
else if ($3>$6 && $3 <$7) {print $2,$3,$7,$1} \
else if ($4>$6 && $4<$7) {print $2,$6,$4,$1} \
}' | sort -k4,4 > lift/5kb/enh_left_right_merge.bed

# Join the full 5k liftover with merged arms -- any successful left/right liftover will need to have a gap in the middle less than the size of the full 5k liftover (535 regions)

sort -k4,4 lift/5kb/enh_ext_lift.bed | join -j4 -t $'\t' lift/5kb/enh_left_right_merge.bed - | awk 'OFS="\t" {if ($7-$6 >= $4-$3) print $2,$3,$4,$1}' > lift/5kb/enh_combined_lift.bed


# Identify Both on enhancers, one on enhancers, same as above
grep L lift/5kb/enh_combined_lift.bed | bedtools intersect -wo -b - -a lift/enh.bed | cut -f1-4,8 | sort -k5,5 | join -1 5 -2 4 -o 1.1,1.2,1.3,1.4,2.1,2.2,2.3,2.4 -t $'\t' - lift/enh.bed > lift/5kb/both_on_redundant.bed_temp

grep S lift/5kb/enh_combined_lift.bed | bedtools intersect -wo -b - -a lift/enh.bed | cut -f1-4,8 | sort -k5,5 | join -1 4 -2 5 -o 1.1,1.2,1.3,1.4,2.1,2.2,2.3,2.4 -t $'\t' lift/enh.bed - >> lift/5kb/both_on_redundant.bed_temp

cat <(cut -f4 lift/5kb/both_on_redundant.bed_temp) <(cut -f8 lift/5kb/both_on_redundant.bed_temp) | sort | uniq | grep -v -f lift/500bp/both_on_ids.txt > lift/5kb/both_on_ids.txt

# Make a file that has the actual lifted over both-on regions
grep -f lift/5kb/both_on_ids.txt lift/5kb/enh_combined_lift.bed > lift/5kb/both_on_lifted_regions.bed

# Identify L on (N = 299)
grep S lift/5kb/enh_combined_lift.bed | bedtools intersect -v -a - -b lift/enh.bed | sort -k4,4 | join -t $'\t' -j4 -o 1.1,1.2,1.3,1.4,2.1,2.2,2.3 lift/enh.bed - > lift/5kb/L_on_rescue.bed

# Identify S on (N = 220)
grep L lift/5kb/enh_combined_lift.bed | bedtools intersect -v -a - -b lift/enh.bed | sort -k4,4 | join -t $'\t' -j4 -o 1.1,1.2,1.3,1.4,2.1,2.2,2.3 lift/enh.bed - > lift/5kb/S_on_rescue.bed

LIFT5K


# Join 500bp and 5k liftovers
:<<"COMBINE"

# Combined regions, both on; but some regions listed multiple times (redundant)
# Each pair is given an ID derived from the individual enhancer IDs

cat lift/500bp/both_on_redundant.bed lift/5kb/both_on_redundant.bed_temp | awk -v OFS="\t" '{print "LS-"$4"-"$8,$0}' | sort | uniq > lift/combined/both_on_redundant.txt_temp

# Sort based on both L and S ids (ensures pairs involving distals are first)
cut -f5,9 lift/combined/both_on_redundant.txt_temp | awk '{if ($1 < $2) {print $1} else {print $2}}' | paste lift/combined/both_on_redundant.txt_temp - > lift/combined/both_on_redundant.txt_temp2
sort -k10,10 lift/combined/both_on_redundant.txt_temp2 | cut -f1-9 > lift/combined/both_on_redundant.txt
rm lift/combined/both_on_redundant.txt_temp*

# For visualization/statistics purposes, make a non-redundant file of LS paired enhancers
# where each enhancer appears only once (one to one)
# (this involves column rearrangement and uniqs)

# Uniq first on S ID (column 8), sort by L ID (column 5)
sort -k9,9 lift/combined/both_on_redundant.txt | uniq -f8 | sort -k5,5 > lift/combined/both_tempLS.txt

# Move the L id to the end, uniq on the L ID, restore column order
awk -v OFS="\t" '{print $1,$2,$3,$4,$6,$7,$8,$9,$5}' lift/combined/both_tempLS.txt | uniq -f8 | awk 'OFS="\t" {print $1,$2,$3,$4,$9,$5,$6,$7,$8}' > lift/combined/both_on_nonredundant.txt_temp

# Sort based on both L and S ids (ensures pairs involving distals are first)
cut -f5,9 lift/combined/both_on_nonredundant.txt_temp | awk '{if ($1 < $2) {print $1} else {print $2}}' | paste lift/combined/both_on_nonredundant.txt_temp - > lift/combined/both_on_nonredundant.txt_temp2
sort -k10,10 lift/combined/both_on_nonredundant.txt_temp2 | cut -f1-9 > lift/combined/both_on_nonredundant.txt

rm lift/combined/both_tempLS.txt lift/combined/both_on_nonredundant.txt_temp*

# Make a list of both on enhancer IDs
cat <(cut -f5 lift/combined/both_on_redundant.txt) <(cut -f9 lift/combined/both_on_redundant.txt) | sort | uniq > lift/combined/both_on_ids.txt

# Single on enhancers (create a 'NA' ID for the homeologous region)
# L on S off
cat lift/500bp/L_on.bed lift/5kb/L_on_rescue.bed | grep -w -v -f lift/combined/both_on_ids.txt | awk -v OFS="\t" '{print "L-"$4"-NA",$0,"NA"}' > lift/combined/L_on.txt

# S on L off
cat lift/500bp/S_on.bed lift/5kb/S_on_rescue.bed | grep -w -v -f lift/combined/both_on_ids.txt | awk -v OFS="\t" '{print "S-NA-"$4,$0,"NA"}'> lift/combined/S_on.txt

# Concatenate all the files together (maintaining L in cols 2-5 and S in cols 6-9)
cat lift/combined/both_on_redundant.txt lift/combined/L_on.txt <(awk -v OFS="\t" '{print $1,$6,$7,$8,$9,$2,$3,$4,$5}' lift/combined/S_on.txt) > lift/combined/all_enh_homeologs.txt_temp

# Identify the enhancers that failed to liftover, append to master list
cat <(cut -f5 lift/combined/all_enh_homeologs.txt_temp) <(cut -f9 lift/combined/all_enh_homeologs.txt_temp) | grep -v -f - lift/enh.bed > lift/combined/unlifted.bed
grep L lift/combined/unlifted.bed | awk -v OFS="\t" '{print "l-"$4,$1,$2,$3,$4,"NA","NA","NA","NA"}' > lift/combined/unlifted.bed_templ
grep S lift/combined/unlifted.bed | awk -v OFS="\t" '{print "s-"$4,"NA","NA","NA","NA",$1,$2,$3,$4}' > lift/combined/unlifted.bed_temps
cat lift/combined/all_enh_homeologs.txt_temp lift/combined/unlifted.bed_templ lift/combined/unlifted.bed_temps > lift/combined/all_enh_homeologs.txt

rm lift/combined/unlifted.bed_temp* lift/combined/all_enh_homeologs.txt_temp


# GRAND TOTAL (involving 7562 total distal enhancers)
# 865+909=1774 conserved distal enhancers (226 L and 272 S proximal elements match to distal elements)
# 2561 L on S off enhancers
# 1776 S on L off enhancers
# 1451 failed to liftover nicely (999 L, 452 S)
COMBINE


# UCSC bigbed
:<<"BIGBED"
egrep "^LS|^L-" lift/combined/all_enh_homeologs.txt | awk -v OFS="\t" '{print $2,$3,$4,$1,1000}' > lift/combined/enh_ucsc.bed
egrep "^LS|^S-" lift/combined/all_enh_homeologs.txt | awk -v OFS="\t" '{print $6,$7,$8,$1,1000}' >> lift/combined/enh_ucsc.bed

grep "^L-" lift/combined/all_enh_homeologs.txt | awk -v OFS="\t" '{print $6,$7,$8,$1,300}' >> lift/combined/enh_ucsc.bed
grep "^S-" lift/combined/all_enh_homeologs.txt | awk -v OFS="\t" '{print $2,$3,$4,$1,300}' >> lift/combined/enh_ucsc.bed

sort -k1,1 -k2,2n lift/combined/enh_ucsc.bed > lift/combined/enh_ucsc_sorted.bed
bedToBigBed -type=bed5 lift/combined/enh_ucsc_sorted.bed $CHRSIZES lift/combined/enh_ucsc.bigbed

BIGBED


# Make versions of the BED4 files centered on each region

#:<<"SEQFASTA"

function BEDMIDPT () {
	awk -v OFS="\t" '{midpt=int(($3+$2)/2); print $1,midpt,midpt+1,$4}' $1	
}

# both on distal regions (maintain sort order) (non redundant and redundant versions)
grep distal lift/combined/both_on_nonredundant.txt | awk -v OFS="\t" '{print $2,$3,$4,$1}' | BEDMIDPT > lift/combined/distal_both_on_L_1bp.bed
grep distal lift/combined/both_on_nonredundant.txt | awk -v OFS="\t" '{print $6,$7,$8,$1}' | BEDMIDPT > lift/combined/distal_both_on_S_1bp.bed

grep distal lift/combined/both_on_redundant.txt | awk -v OFS="\t" '{print $2,$3,$4,$1}' | BEDMIDPT | sort -k1,1 -k2,2n | rev | uniq -f1 | rev > lift/combined/distal_both_on_L_1bp_allunique.bed
grep distal lift/combined/both_on_redundant.txt | awk -v OFS="\t" '{print $6,$7,$8,$1}' | BEDMIDPT | sort -k1,1 -k2,2n | rev | uniq -f1 | rev > lift/combined/distal_both_on_S_1bp_allunique.bed

# Single on regions
grep distal lift/combined/L_on.txt | awk -v OFS="\t" '{print $2,$3,$4,$1}' | BEDMIDPT > lift/combined/distal_L_on_L_1bp.bed
grep distal lift/combined/L_on.txt | awk -v OFS="\t" '{print $6,$7,$8,$1}' | BEDMIDPT > lift/combined/distal_L_on_S_1bp.bed

grep distal lift/combined/S_on.txt | awk -v OFS="\t" '{print $2,$3,$4,$1}' | BEDMIDPT > lift/combined/distal_S_on_S_1bp.bed
grep distal lift/combined/S_on.txt | awk -v OFS="\t" '{print $6,$7,$8,$1}' | BEDMIDPT > lift/combined/distal_S_on_L_1bp.bed


# Generate sequence files for these regions of different sizes

# +/- 100bp from center for Homer motif finding
bedtools slop -b 100 -g $CHRSIZES -i lift/combined/distal_both_on_L_1bp_allunique.bed | bedtools getfasta -nameOnly -fi $GENOME -bed - -fo lift/combined/distal_both_on_L_200bp.fa
bedtools slop -b 100 -g $CHRSIZES -i lift/combined/distal_both_on_S_1bp_allunique.bed | bedtools getfasta -nameOnly -fi $GENOME -bed - -fo lift/combined/distal_both_on_S_200bp.fa

bedtools slop -b 100 -g $CHRSIZES -i lift/combined/distal_L_on_L_1bp.bed | grep distal | bedtools getfasta -nameOnly -fi $GENOME -bed - -fo lift/combined/distal_L_on_L_200bp.fa
bedtools slop -b 100 -g $CHRSIZES -i lift/combined/distal_L_on_S_1bp.bed | grep distal | bedtools getfasta -nameOnly -fi $GENOME -bed - -fo lift/combined/distal_L_on_S_200bp.fa

bedtools slop -b 100 -g $CHRSIZES -i lift/combined/distal_S_on_L_1bp.bed | grep distal | bedtools getfasta -nameOnly -fi $GENOME -bed - -fo lift/combined/distal_S_on_L_200bp.fa
bedtools slop -b 100 -g $CHRSIZES -i lift/combined/distal_S_on_S_1bp.bed | grep distal | bedtools getfasta -nameOnly -fi $GENOME -bed - -fo lift/combined/distal_S_on_S_200bp.fa


exit
