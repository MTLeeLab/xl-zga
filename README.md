# Hybridization led to a rewired pluripotency network in the allotetraploid Xenopus laevis - Code

This repository will contain the unix and R scripts used to conduct the computational analyses within the paper "Hybridization led to a rewired pluripotency network in the allotetraploid Xenopus laevis" published at ...

## Description

Our project aims to understand how recent genome upheaval has affected the pluripotency gene regulatory network in the allotetraploid Xenopus laevis. Here we will provide the unix and R scripts used to compare L&S subgenome gene expression and regulation, as well as our homology analysis with Xenopus tropicalis and zebrafish. These scripts will also output the data visualizations we've presented in each figure.

## Getting Started

### Dependencies

* Programs: PAML v4.9f, pal2nal v14, EMBOSS needle v6.6.0.0, R v4.0.4, BEDtools v2.30.0, Samtools v1.12, deepTools v3.5.1, MACS2 v2.2.7.1, SEACR v1.3, Homer v4.11.1, Genome Browser utilities, python3.
* Annotation files: Xenbase v9.2 gtf, Xenbase GO annotations, X. laevis v9.2 genome, X. tropicalis v7, homer known vertebrate motif database, zebrafish gene annotations from Lee & Bonneau et al (GEO GSE47558).
* Chain files: danRer10ToDanRer11, xenTro2ToXenTro7, xenTro7ToXenTro9, xenLae2ToXenTro9, xenTro7ToDanRer10, and xenTro9ToXenTro7 (from UCSC genome browser). xlL_xlS_med (from this study at OFS ...).
* Data files: all provided bigwig, bed, bam, and txt files from this study (deposited at GEO ...), published stage 9 H3K27ac peaks from Gupta et al (GEO: GSM1350514 & GSM1350515), published dome H3K27ac data from Bogdanovic et al (GEO: GSM915197, align and use BAM file).

### Installing

* Download all unix scripts and the scripts_R directory into the working directory. All required data files can be downloaded from the provided GEO accession pages. All xenopus annotation files can be downloaded directly from xenbase.org and known homer motifs can be downloaded from homer.ucsd.edu/homer/.
* As of 9/1/22, the provided python script must be manually altered at the last line and run three times total to properly match the both on total, L on total, and S on total files that are generated in the enh cmds script.

### Executing program

*  Within the working directory, create a 'files' directory and populate it with the directories 'anno', 'bams', 'beds', 'bigwigs', and 'chains'. All above annotation files should be placed in 'anno' and all data files provided/required in this study should be placed in the appropriate files directory before running any scripts.
* Scripts should be executable from the command line at the working directory level and will call individual R scripts from the scripts_R directory (also in the working directory). Scripts will continually be updated until they are finalized at the date of manuscript publication. 
* If the unix scripts are not executable for whatever reason, the following command will make them executable:
```
chmod u+x $SCRIPT
```

## Help

Helper advice will be added here as issues arise
```
command to run if program contains helper info
```

## Authors

Wesley A. Phelps (wap22@pitt.edu)
Miler T. Lee (miler@pitt.edu)

## Version History

* 0.1
    * Initial Release

## License

This project is licensed under the MTLeeLab/xl-zga License - see the LICENSE.md file for details

## Acknowledgments

Inspiration, code snippets, etc.
* [DomPizzie/README-Template.md](https://gist.github.com/DomPizzie/7a5ff55ffa9081f2de27c315f5018afc)
