# -----------------------------------------------------------------------------------------------------------
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
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# The following packages should be load before running the analysis
# -----------------------------------------------------------------------------------------------------------
list.of.packages <- c(
  "clusterProfiler",
  "devtools",
  "dplyr",
  "doParallel",
  "genefilter",
  "GEOquery",
  "glmnet",
  "ggplot2",
  "ggpubr",
  "janitor",
  "matrixStats",
  "msigdbr", 
  "plyr",
  "purrr",
  "pheatmap",
  "randomForestSRC",
  "readr",
  "rstatix",
  "survival",
  "survminer",
  "SummarizedExperiment",
  "TCGAbiolinks",
  "tibble",
  "WGCNA"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)){
  for(new in new.packages){
    if(new %in% available.packages()[,1]){
      install.packages(new)
    } else BiocManager::install(new)
  }
} 