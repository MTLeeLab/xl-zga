# Hybridization led to a rewired pluripotency network in the allotetraploid Xenopus laevis - Analysis scripts

This repository contains unix and R scripts used to conduct the computational analyses for the first submission of Phelps et al 2023, Hybridization led to a rewired pluripotency network in the allotetraploid Xenopus laevis, eLife https://doi.org/10.7554/eLife.83952 .

### Dependencies

* Programs: PAML v4.9f, pal2nal v14, EMBOSS needle v6.6.0.0, R v4.0.4, BEDtools v2.30.0, Samtools v1.12, deepTools v3.5.1, MACS2 v2.2.7.1, SEACR v1.3, Homer v4.11.1, Genome Browser utilities, python3.
* Annotation files: Xenbase v9.2 gtf, Xenbase GO annotations, X. laevis v9.2 genome, X. tropicalis v7, homer known vertebrate motif database, zebrafish gene annotations from Lee & Bonneau et al (GEO GSE47558).
* Liftover chain files: danRer10ToDanRer11, xenTro2ToXenTro7, xenTro7ToXenTro9, xenLae2ToXenTro9, xenTro7ToDanRer10, and xenTro9ToXenTro7 (from UCSC genome browser). xlL_xlS_med (from this study at OFS ...).
* Data files: all provided bigwig, bed, bam, and txt files from this study (deposited at GEO ...), published stage 9 H3K27ac peaks from Gupta et al (GEO: GSM1350514 & GSM1350515), published dome H3K27ac data from Bogdanovic et al (GEO: GSM915197, align and use BAM file).

### Notes

* Within the working directory, create a 'files' directory and populate it with the directories 'anno', 'bams', 'beds', 'bigwigs', and 'chains'. All above annotation files should be placed in 'anno' and all data files provided/required in this study should be placed in the appropriate files directory before running any scripts.


## Authors

Wesley A. Phelps (wap22@pitt.edu)
Miler T. Lee (miler@pitt.edu)

## Version History

* 0.1
    * Analyses for first submission

## License

This project is licensed under the MTLeeLab/xl-zga License - see the LICENSE.md file for details

## Acknowledgments

Inspiration, code snippets, etc.
* [DomPizzie/README-Template.md](https://gist.github.com/DomPizzie/7a5ff55ffa9081f2de27c315f5018afc)
