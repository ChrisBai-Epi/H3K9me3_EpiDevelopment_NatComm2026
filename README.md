# H3K9me3_EpiDevelopment_NatComm2026

This repository contains analysis scripts used in the study investigating the role of H3K9me3 in epidermal development (Bai et al., 2026).

## Repository structure

scRNA-seq/<br>
R scripts for scRNA-seq data preprocessing, clustering, differential analysis, GO analysis and related figure generation.<br>
<br>
bulk-RNA-seq/<br>
R scripts for bulk RNA-seq data analysis with DESeq2, GO analysis and figure generation.<br>
<br>
PReCIS-seq/<br>
R scripts for PReCIS-seq data preprocessing, differential analysis, pausing anlaysis, GO analysis and related figure generation.<br>
<br>
CUT&RUN_H3K9me3_E15/<br>
R scripts for analyzing CUT&RUN data for E15.5 epidermal cells. Analyses include methylation pattern clustering, distance analysis, GO analysis and related figure generation.<br>
<br>
CUT&RUN_H3K9me3_Temporal/<br>
R scripts for analyzing temporal CUT&RUN data for E12.5, E14.5 and E16.5. Analyses include dynamic peak analysis, gene-peak association analysis, GO analysis and related figure generation.<br>


## Data availibility
Data for reproducing the analysis is available from GEO as described in the study.

## Notes

This code is intended to support reproducibility of the analyses presented in
the paper and is not designed as a standalone software package.

## License

This project is licensed under the MIT License.
