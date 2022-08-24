# commands to generate the rna lfc biplots and paired L and S heatmaps
# add homeolog lfc table to files/anno
#### uses xl_gene_chrs.txt. Find that code! ####

# make output directory
mkdir -p rna/pdfs

# run R script to generate plots and ls rna lfc data table (used in promoter plots)
Rscript scripts_R/ls_rna_seq_plots.R

exit
