# Enhancing Chemotherapy Response Prediction via Matched Colorectal Tumor-Organoid Gene Expression Analysis and Network-Based Biomarker Selection 

Wei Zhang, Chao Wu, Hanchen Huang, Paulina Bleu, Wini Zambare, Janet Alvarez, Lily Wang, Philip B. Paty, Paul B. Romesser, J. Joshua Smith, X. Steven Chen



This github repository includes scripts used for the analyses in the above manuscript.

## Abstract

This study presents an innovative methodology for predicting chemotherapy responses in colorectal cancer patients by integrating gene expression data from matched colorectal tumor and organoid samples. Employing Consensus Weighted Gene Co-expression Network Analysis (WGCNA) across multiple datasets, we identified key gene modules and hub genes linked to patient response to chemotherapy, focusing on 5-fluorouracil (5-FU). This integrative approach marks a significant advancement in precision medicine, enhancing the specificity and accuracy of chemotherapy regimen selection based on individual tumor profiles. Our predictive model, validated by independent datasets demonstrated improved accuracy over traditional methods. This strategy shows promise in overcoming typical challenges in high-dimensional genomic data analysis for cancer biomarker research. 

## Folder Directory

| File                    | Description                                        |
| ----------------------- | -------------------------------------------------- |
| load_package.R          | Session information and all required packages      |
| code/data_preprocessing | Preprocessing of all the data used in the analysis |
| code/WGCNA              | Consensus WGCNA model                              |
| code/organoid_model     | Organoid drug-response model training              |
| code/others             | Other gene selection processes                     |
| code/utility            | Functions used for analysis                        |
| results                 | Genes in consensus WGCNA key modules               |

## For Reproducible Research

Install all the R packages from the load_package.R file.

The platform information is:

```R
# For reproducible research, please install the following R packages 
# and make sure the R and BiocManager versions are correct
# Session Info ----------------------------------------------------------------------------------------------
# setting  value 
# version  R version 4.3.1 (2023-06-16)
# os       macOS Sonoma 14.3
# system   x86_64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/New_York
# date     2024-01-27
# rstudio  2023.09.1+494 Desert Sunflower (desktop)
# pandoc   3.1.6.1 @ /usr/local/bin/ (via rmarkdown)
```

