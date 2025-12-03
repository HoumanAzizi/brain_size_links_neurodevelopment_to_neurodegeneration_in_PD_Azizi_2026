## NOTE: should gunzip and gzip before and after


rm(list = ls())
cat("\014")

library(dplyr)

setwd('/Users/HoumanAzizi/Library/Mobile\ Documents/com~apple~CloudDocs/Education/GitHub/Mendelian_Randomization_UKB/')


GWAS <- read.delim("Data/GP2_PD_GWAS/GP2_euro_ancestry_meta_analysis_2024_originalDownload/GP2_ALL_EUR_ALL_DATASET_HG38_12162024.rsid.txt", header = TRUE, sep = "\t", quote = "", fill = TRUE)


N_case <- 63555
N_proxy <- 17700
N_case_total <- N_case + N_proxy
N_control <- 1746386


GWAS$N_case <- N_case
GWAS$N_proxy <- N_proxy
GWAS$N_case_proxy <- N_case_total
GWAS$N_control <- N_control


write.table(GWAS, "Data/GP2_PD_GWAS/GP2_ALL_EUR_ALL_DATASET_HG38_12162024.rsid_withN.txt", sep = "\t", row.names = FALSE, quote = FALSE)