# commands for selecting TSS sites and homeolog coverage of k4me3
# generate a TSS file for each transcript
# add GTF from xenbase and genome file to files/anno

#### we converted gtf to ucsc nomenclature! Find that code! ####

TAB=$'\t'

GTF=files/anno/XENLA_9.2_Xenbase_ucsc.gtf
GENOME=files/anno/genome.fa

:<<TSSANNO

# make chromosome size file
samtools faidx $GENOME
cut -f1,2 files/anno/genome.fa.fai > files/anno/chr_sizes.txt
SIZE=files/anno/chr_sizes.txt

# make output directory
mkdir -p tss

# Full TSS bed file from GTF
grep "Xenbase${TAB}transcript" $GTF | sed 's/ gene_name.*$//' | awk 'OFS="\t" {print $1,$4-1,$5,$NF,".",$7}' | sed 's/[\";]//g' | sort -k4,4 | bedtools flank -l 1 -r 0 -s -i stdin -g $SIZE > tss/transcript_TSS.bed

# TSS subset choosing most likely isoform

# calculate zygotic RNA-seq coverage in a 50bp window into the gene body of each TSS
DMSO1=files/bams/s9.5_dmso_a1
DMSO2=files/bams/s9.5_dmso_a2
DMSO3=files/bams/s9.5_dmso_d1
DMSO4=files/bams/s9.5_dmso_d2

bedtools slop -l 0 -r 50 -s -i tss/transcript_TSS.bed -g $SIZE | sort -k1,1 -k2,2n | bedtools coverage -sorted -s -a stdin -b $DMSO1 $DMSO2 $DMSO3 $DMSO4 -g $SIZE | sort -k4,4 > tss/TSS_50bp_cov.bed

# annotate nonzero RNA-seq coverage TSSs
awk -v OFS="\t" -v FS="\t" '{if ($7 > 0) $5="X"; print $1,$2,$3,$5,$6,$7,$8,$9,$10,$4;}' tss/TSS_50bp_cov.bed > tss/TSS_nonzero_anno

# Pull out the most upstream, nonzero TSS for each gene (sorted by end coord for -, start coord for +)
grep -w '-' tss/TSS_nonzero_anno | sort -k10,10 -k4,4r -k3,3nr | uniq -f9 | awk -v OFS="\t" -v FS="\t" '{print($1,$2,$3,$10,$4,$5,$6,$7,$8,$9)}' > tss/TSS_nonzero_neg.txt_

grep -w '+' tss/TSS_nonzero_anno | sort -k10,10 -k4,4r -k2,2n | uniq -f9 | awk -v OFS="\t" -v FS="\t" '{print($1,$2,$3,$10,$4,$5,$6,$7,$8,$9)}' > tss/TSS_nonzero_pos.txt_

cat tss/TSS_nonzero_pos.txt_ tss/TSS_nonzero_neg.txt_ | bedtools flank -l 1 -r 0 -s -i stdin -g $SIZE | bedtools flank -l 0 -r 1 -s -i stdin -g $SIZE > tss/TSS_select.bed
rm tss/TSS_nonzero_neg.txt_ tss/TSS_nonzero_pos.txt_

TSSANNO


exit
