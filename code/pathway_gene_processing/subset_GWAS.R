# Load required libraries
library(dplyr)

rm(list = ls())
cat("\014")

setwd("/Users/houmanazizi/Library/Mobile\ Documents/com~apple~CloudDocs/Education/GitHub/Pathway_PD_PRS/")




#### Read PD GWAS
PD_GWAS <- "/Users/HoumanAzizi/Library/Mobile\ Documents/com~apple~CloudDocs/Education/GitHub/Mendelian_Randomization_UKB/Data/GP2_PD_GWAS/GP2_ALL_EUR_ALL_DATASET_HG38_12162024.rsid_withN.txt"
PD_GP2 <- read.table(PD_GWAS, sep = '\t', header = TRUE)


#### Read pathway SNPs
all_pathways_snps <- read.table("/Users/HoumanAzizi/Library/Mobile\ Documents/com~apple~CloudDocs/Education/GitHub/UKB_Analysis/Data/Genetic/PD_PRS/PD_GP2_PRScs_pathway_SNPs/Mitochondrial_Lysosomal_Autophagy_direct_genes_9_PD_GWAS.snplist", 
                                header = TRUE, stringsAsFactors = FALSE)




#### Subset the GWAS
PD_GWAS_nonPathway <- PD_GP2 %>% filter( !(rsID %in% all_pathways_snps$rsID) )
write.table(PD_GWAS_nonPathway, 
            file = "Outputs/GWAS_subset/GP2_ALL_EUR_ALL_DATASET_HG38_12162024_nonPathwaySubset.txt", 
            sep = '\t', 
            row.names = FALSE, 
            quote = FALSE)

PD_GWAS_onlyPathway <- PD_GP2 %>% filter(rsID %in% all_pathways_snps$rsID)
write.table(PD_GWAS_onlyPathway, 
            file = "Outputs/GWAS_subset/GP2_ALL_EUR_ALL_DATASET_HG38_12162024_onlyPathwaySubset.txt", 
            sep = '\t', 
            row.names = FALSE, 
            quote = FALSE)




#### Prepare GWAS for MAGMA by removing "" rsID values
PD_GP2_clean <- PD_GP2 %>% filter(rsID != "" & !is.na(rsID))
write.table(PD_GP2_clean, 
            file = "/Users/HoumanAzizi/Library/Mobile\ Documents/com~apple~CloudDocs/Education/GitHub/UKB_Analysis/Data/Genetic/GWAS/GP2_ALL_EUR_ALL_DATASET_HG38_12162024.rsid_withN_rsIDcleaned.txt", 
            sep = '\t', 
            row.names = FALSE, 
            quote = FALSE)

PD_GWAS_nonPathway_clean <- PD_GWAS_nonPathway %>% filter(rsID != "" & !is.na(rsID))
write.table(PD_GWAS_nonPathway_clean, 
            file = "Outputs/GWAS_subset/GP2_ALL_EUR_ALL_DATASET_HG38_12162024_nonPathwaySubset_rsIDcleaned.txt", 
            sep = '\t', 
            row.names = FALSE, 
            quote = FALSE)

PD_GWAS_onlyPathway_clean <- PD_GWAS_onlyPathway %>% filter(rsID != "" & !is.na(rsID))
write.table(PD_GWAS_onlyPathway_clean, 
            file = "Outputs/GWAS_subset/GP2_ALL_EUR_ALL_DATASET_HG38_12162024_onlyPathwaySubset_rsIDcleaned.txt", 
            sep = '\t', 
            row.names = FALSE, 
            quote = FALSE)