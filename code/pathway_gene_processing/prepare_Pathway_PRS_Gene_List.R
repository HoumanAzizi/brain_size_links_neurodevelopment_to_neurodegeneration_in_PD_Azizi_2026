library(tidyr)
library(dplyr)
library(readxl)

# if (!require("BiocManager")) install.packages("BiocManager")
# BiocManager::install("biomaRt")
# library(biomaRt)

rm(list = ls())
cat("\014")

setwd("/Users/houmanazizi/Library/Mobile\ Documents/com~apple~CloudDocs/Education/GitHub/Pathway_PD_PRS/")


########## Read the original data ##########
# Mitochondrial genes
mito_Carta <- read_xlsx('Data/Original_gene_lists/Human_MitoCarta3_genesOnly.xlsx', col_types = 'text')
mito_Billingsley <- read_xlsx('Data/Original_gene_lists/MitoAssociatedGenes_Billingsley.xlsx')
# Lysosomal genes
lyso_GO_direct <- read.csv('Data/Original_gene_lists/Lysosomal_GO.csv')
lyso_GO_any <- read.csv('Data/Original_gene_lists/Lysosomal_GO_withChildTerms.csv')
# Autophagy genes
autophagy_GO_direct <- read.csv('Data/Original_gene_lists/Autophagy_GO.csv')
autophagy_GO_any <- read.csv('Data/Original_gene_lists/Autophagy_GO_withChildTerms.csv')
autophagy_HADb <- read_xlsx('Data/Original_gene_lists/HADb_Autophagy_Genes.xlsx')
# Lysosomal AND Autophagy genes
# lyso_autophagy_GO_direct <- read_xlsx('Data/Original_gene_lists/Autophagy_Lysosome_GO_directAnnotation.xlsx') # now calculating it directly in the code





########## Function to map genes to BP position ##########
get_positions_mixed <- function(gene_list, entrez_ids = FALSE, ensembl_ids = FALSE, gene_names = FALSE) {
  library(biomaRt) # loading biomaRT and unloading at the end of function because it prevents dplyr to work
  
  # Get gene databases
  # previous method - server problem
  # ensembl38 <- useMart(biomart = "ensembl", dataset="hsapiens_gene_ensembl")
  # ensembl37 <- useMart(host="https://grch37.ensembl.org", biomart="ensembl", dataset="hsapiens_gene_ensembl")
  # new method
  ensembl38 <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "asia")
  ensembl37 <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37)
  
  
  # Check if exactly one method is selected
  selected_methods <- sum(c(entrez_ids, ensembl_ids, gene_names))
  
  if (selected_methods != 1) {
    stop("Error: Please select exactly ONE method (entrez_ids, ensembl_ids, or gene_names)")
  }
  
  # Determine which filter to use based on method
  if (entrez_ids) {
    filter_type <- "entrezgene_id"
    join_column <- "entrezgene_id"
  } else if (ensembl_ids) {
    filter_type <- "ensembl_gene_id"
    join_column <- "ensembl_gene_id"
  } else if (gene_names) {
    filter_type <- "hgnc_symbol"
    join_column <- "hgnc_symbol"
  }
  
  # Get positions from both versions
  results38 <- getBM(attributes = c("entrezgene_id", 
                                    "ensembl_gene_id",
                                    "hgnc_symbol", 
                                    "chromosome_name",
                                    "start_position", 
                                    "end_position"),
                     filters = c(filter_type, 
                                 "chromosome_name"),
                     values = list(gene_list, 
                                   c(1:22, "X")),
                     mart = ensembl38)
  results38 <- results38[!duplicated(results38$ensembl_gene_id), ]
  
  results37 <- getBM(attributes = c("entrezgene_id", 
                                    "ensembl_gene_id",
                                    "hgnc_symbol", 
                                    "chromosome_name",
                                    "start_position", 
                                    "end_position"),
                     filters = c(filter_type, 
                                 "chromosome_name"),
                     values = list(gene_list, 
                                   c(1:22, "X")),
                     mart = ensembl37)
  results37 <- results37[!duplicated(results37$ensembl_gene_id), ]
  results37 <- results37[, c("ensembl_gene_id", "chromosome_name", "start_position", "end_position")]
  
  # Add version labels to position columns
  pos_cols <- c("chromosome_name", "start_position", "end_position")
  colnames(results38)[colnames(results38) %in% pos_cols] <- 
    paste0(colnames(results38)[colnames(results38) %in% pos_cols], "_GRCh38")
  colnames(results37)[colnames(results37) %in% pos_cols] <- 
    paste0(colnames(results37)[colnames(results37) %in% pos_cols], "_GRCh37")
  
  # Full join the results
  merged_results <- merge(results38, results37, 
                          by = 'ensembl_gene_id', 
                          all = TRUE)
  
  # Order by chromosome and position (using GRCh38, if available)
  merged_results <- merged_results[order(merged_results$chromosome_name_GRCh38, 
                                         merged_results$start_position_GRCh38), ]
  
  rm(ensembl37, ensembl38)
  detach("package:biomaRt", unload = TRUE)
  return(merged_results)
}











########## Prepare and merge genes ##########
### Gene list descriptions:
## Single processes
# 1. Mitochondrial only -> mito_Carta + mito_Billingsley
# 2. Lysosomal direct GO annotation only -> lyso_GO_direct
# 3. Lysosomal any GO annotation only -> lyso_GO_any
# 4. Autophagy direct GO only -> autophagy_HADb + autophagy_GO_direct
# 5. Autophagy any GO only -> autophagy_HADb + autophagy_GO_any
# 6. Autophagy AND Lysosomal direct GO only -> lyso genes that are ALSO autophagy (GO:0006914) based on direct GO terms
## Grouped processes
# 7. Mitochondrial + lysosomal_direct -> 1 + 2
# 8. Mitochondrial + lysosomal_any -> 1 + 3
# 9. Mitochondrial + lysosomal_direct + autophagy_direct -> 1 + 2 + 4
# 10. Mitochondrial + lysosomal_any + autophagy_any -> 1 + 3 + 5
# 11. Mitochondrial + autophagy-lysosomal -> 1 + 6 (i.e. taking genes that belong to BOTH lysosomal and autophagy terms)
# 12. Mitochondrial + nonAutophagy-lysosomal -> 1 + 2 - 6 (i.e. taking genes that belong to lysosomal GO term, but not autophagy GO term)
# 13. nonAutophagy-lysosomal only -> 2 - 6
# 14. Autophagy OR Lysosomal direct GO -> 2 + 4
# 15. Autophagy OR Lysosomal any GO -> 3 + 5


## 1 ##
## 1. Mitochondrial only -> mito_Carta + mito_Billingsley
# fix the multiple Ensembl ID cases
mito_Carta$EnsemblGeneID <- sapply(mito_Carta$EnsemblGeneID_mapping_version_20200130, function(x) {
  if (grepl("\\|", x)) { strsplit(x, "\\|")[[1]][1] } else { x } })

# add BP locations
mito_Carta_withBP <- get_positions_mixed(gene_list = mito_Carta$EnsemblGeneID, ensembl_ids = TRUE)
mito_Carta <- mito_Carta %>% left_join(mito_Carta_withBP, by = c('EnsemblGeneID'='ensembl_gene_id'))
rm(mito_Carta_withBP)

mito_Billingsley_withBP <- get_positions_mixed(mito_Billingsley$gene, gene_names = TRUE)
# mito_Billingsley_withBP <- mito_Billingsley_withBP %>% filter(!is.na(hgnc_symbol)) %>% filter(ensembl_gene_id != 'ENSG00000293451')
mito_Billingsley <- mito_Billingsley %>% left_join(mito_Billingsley_withBP, by = c('gene'='hgnc_symbol'))
rm(mito_Billingsley_withBP)

# get missing GRCh37 position from the original file if missing
mito_Billingsley <- mito_Billingsley %>%
  mutate(chromosome_name_GRCh37 = coalesce(chromosome_name_GRCh37, chr),
         start_position_GRCh37 = coalesce(start_position_GRCh37, start),
         end_position_GRCh37 = coalesce(end_position_GRCh37, stop))

# Keep important columns + rename
mito_Carta <- mito_Carta %>% mutate(gene_names = Symbol,
                                    ensembl_ID = EnsemblGeneID,
                                    chr_GRCh37 = chromosome_name_GRCh37,
                                    start_GRCh37 = start_position_GRCh37,
                                    end_GRCh37 = end_position_GRCh37) %>% select(gene_names, ensembl_ID, chr_GRCh37, start_GRCh37, end_GRCh37) %>% filter(!is.na(start_GRCh37)) %>% unique()
mito_Billingsley <- mito_Billingsley %>% mutate(gene_names = gene,
                                                ensembl_ID = ensembl_gene_id,
                                                chr_GRCh37 = chromosome_name_GRCh37,
                                                start_GRCh37 = start_position_GRCh37,
                                                end_GRCh37 = end_position_GRCh37) %>% select(gene_names, ensembl_ID, chr_GRCh37, start_GRCh37, end_GRCh37) %>% filter(!is.na(start_GRCh37)) %>% unique()
# merge and save
Mitochondrial_genes <- rbind(mito_Carta, mito_Billingsley) %>% unique()
write.csv(Mitochondrial_genes, 'Outputs/Mitochondrial_genes_1.csv', row.names = FALSE)


## 2 ##
## 2. Lysosomal direct GO annotation only -> lyso_GO_direct
lyso_GO_direct_withBP <- get_positions_mixed(lyso_GO_direct$Gene.stable.ID, ensembl_ids = TRUE)
lyso_GO_direct <- lyso_GO_direct %>% left_join(lyso_GO_direct_withBP, by = c('Gene.stable.ID'='ensembl_gene_id'))
rm(lyso_GO_direct_withBP)
Lysosomal_genes_direct <- lyso_GO_direct %>% mutate(gene_names = Gene.name,
                                                    ensembl_ID = Gene.stable.ID,
                                                    chr_GRCh37 = chromosome_name_GRCh37,
                                                    start_GRCh37 = start_position_GRCh37,
                                                    end_GRCh37 = end_position_GRCh37) %>% select(gene_names, ensembl_ID, chr_GRCh37, start_GRCh37, end_GRCh37) %>% filter(!is.na(start_GRCh37)) %>% unique()
# save
write.csv(Lysosomal_genes_direct, 'Outputs/Lysosomal_direct_genes_2.csv', row.names = FALSE)


## 3 ##
# 3. Lysosomal any GO annotation only -> lyso_GO_any
lyso_GO_any_withBP <- get_positions_mixed(lyso_GO_any$Gene.stable.ID, ensembl_ids = TRUE)
lyso_GO_any <- lyso_GO_any %>% left_join(lyso_GO_any_withBP, by = c('Gene.stable.ID'='ensembl_gene_id'))
rm(lyso_GO_any_withBP)
Lysosomal_genes_any <- lyso_GO_any %>% mutate(gene_names = Gene.name,
                                                    ensembl_ID = Gene.stable.ID,
                                                    chr_GRCh37 = chromosome_name_GRCh37,
                                                    start_GRCh37 = start_position_GRCh37,
                                                    end_GRCh37 = end_position_GRCh37) %>% select(gene_names, ensembl_ID, chr_GRCh37, start_GRCh37, end_GRCh37) %>% filter(!is.na(start_GRCh37)) %>% unique()
# save
write.csv(Lysosomal_genes_any, 'Outputs/Lysosomal_any_genes_3.csv', row.names = FALSE)


## 4 ##
# 4. Autophagy direct GO only -> autophagy_HADb + autophagy_GO_direct

autophagy_HADb_withBP <- get_positions_mixed(autophagy_HADb$Gene, gene_names = TRUE)
autophagy_HADb <- autophagy_HADb %>% left_join(autophagy_HADb_withBP, by = c('Gene'='hgnc_symbol'))
rm(autophagy_HADb_withBP)
autophagy_HADb <- autophagy_HADb %>% mutate(gene_names = Gene,
                                            ensembl_ID = ensembl_gene_id,
                                            chr_GRCh37 = chromosome_name_GRCh37,
                                            start_GRCh37 = start_position_GRCh37,
                                            end_GRCh37 = end_position_GRCh37) %>% select(gene_names, ensembl_ID, chr_GRCh37, start_GRCh37, end_GRCh37) %>% filter(!is.na(start_GRCh37)) %>% unique()

autophagy_GO_direct_withBP <- get_positions_mixed(autophagy_GO_direct$Gene.stable.ID, ensembl_ids = TRUE)
autophagy_GO_direct <- autophagy_GO_direct %>% left_join(autophagy_GO_direct_withBP, by = c('Gene.stable.ID'='ensembl_gene_id'))
rm(autophagy_GO_direct_withBP)
Autophagy_genes_direct <- autophagy_GO_direct %>% mutate(gene_names = Gene.name,
                                                         ensembl_ID = Gene.stable.ID,
                                                         chr_GRCh37 = chromosome_name_GRCh37,
                                                         start_GRCh37 = start_position_GRCh37,
                                                         end_GRCh37 = end_position_GRCh37) %>% select(gene_names, ensembl_ID, chr_GRCh37, start_GRCh37, end_GRCh37) %>% filter(!is.na(start_GRCh37)) %>% unique()
# merge and save
Autophagy_genes_direct <- rbind(Autophagy_genes_direct, autophagy_HADb) %>% unique()
write.csv(Autophagy_genes_direct, 'Outputs/Autophagy_direct_genes_4.csv', row.names = FALSE)


## 5 ##
# 5. Autophagy any GO only -> autophagy_HADb + autophagy_GO_any
autophagy_GO_any_withBP <- get_positions_mixed(autophagy_GO_any$Gene.stable.ID, ensembl_ids = TRUE)
autophagy_GO_any <- autophagy_GO_any %>% left_join(autophagy_GO_any_withBP, by = c('Gene.stable.ID'='ensembl_gene_id'))
rm(autophagy_GO_any_withBP)
Autophagy_genes_any <- autophagy_GO_any %>% mutate(gene_names = Gene.name,
                                                         ensembl_ID = Gene.stable.ID,
                                                         chr_GRCh37 = chromosome_name_GRCh37,
                                                         start_GRCh37 = start_position_GRCh37,
                                                         end_GRCh37 = end_position_GRCh37) %>% select(gene_names, ensembl_ID, chr_GRCh37, start_GRCh37, end_GRCh37) %>% filter(!is.na(start_GRCh37)) %>% unique()
# merge and save
Autophagy_genes_any <- rbind(Autophagy_genes_any, autophagy_HADb) %>% unique()
write.csv(Autophagy_genes_any, 'Outputs/Autophagy_any_genes_5.csv', row.names = FALSE)


## 6 ##
# 6. Autophagy AND Lysosomal direct GO only -> lyso genes that are ALSO autophagy (GO:0006914) based on direct GO terms
L_A_common_genes <- intersect(Lysosomal_genes_direct$ensembl_ID, Autophagy_genes_direct$ensembl_ID)
Lysosomal_AND_Autophagy_genes <- Lysosomal_genes_direct %>% filter(ensembl_ID %in% L_A_common_genes) %>% unique()
write.csv(Lysosomal_AND_Autophagy_genes, 'Outputs/Lysosomal_AND_Autophagy_genes_6.csv', row.names = FALSE)


## 7 ##
# 7. Mitochondrial + lysosomal_direct -> 1 + 2
Mitochondrial_Lysosomal_direct_genes <- rbind(Mitochondrial_genes, Lysosomal_genes_direct) %>% unique()
write.csv(Mitochondrial_Lysosomal_direct_genes, 'Outputs/Mitochondrial_Lysosomal_direct_genes_7.csv', row.names = FALSE)


## 8 ##
# 8. Mitochondrial + lysosomal_any -> 1 + 3
Mitochondrial_Lysosomal_any_genes <- rbind(Mitochondrial_genes, Lysosomal_genes_any) %>% unique()
write.csv(Mitochondrial_Lysosomal_any_genes, 'Outputs/Mitochondrial_Lysosomal_any_genes_8.csv', row.names = FALSE)

## 9 ##
# 9. Mitochondrial + lysosomal_direct + autophagy_direct -> 1 + 2 + 4
Mitochondrial_Lysosomal_Autophagy_direct_genes <- rbind(Mitochondrial_genes, Lysosomal_genes_direct, Autophagy_genes_direct) %>% unique()
write.csv(Mitochondrial_Lysosomal_Autophagy_direct_genes, 'Outputs/Mitochondrial_Lysosomal_Autophagy_direct_genes_9.csv', row.names = FALSE)

## 10 ##
# 10. Mitochondrial + lysosomal_any + autophagy_any -> 1 + 3 + 5
Mitochondrial_Lysosomal_Autophagy_any_genes <- rbind(Mitochondrial_genes, Lysosomal_genes_any, Autophagy_genes_any) %>% unique()
write.csv(Mitochondrial_Lysosomal_Autophagy_any_genes, 'Outputs/Mitochondrial_Lysosomal_Autophagy_any_genes_10.csv', row.names = FALSE)

## 11 ##
# 11. Mitochondrial + autophagy-lysosomal -> 1 + 6 (i.e. taking genes that belong to BOTH lysosomal and autophagy terms)
Mitochondrial_LysoAuto_genes <- rbind(Mitochondrial_genes, Lysosomal_AND_Autophagy_genes) %>% unique()
write.csv(Mitochondrial_LysoAuto_genes, 'Outputs/Mitochondrial_LysoAuto_genes_11.csv', row.names = FALSE)

## 12 ##
# 12. Mitochondrial + nonAutophagy-lysosomal -> 1 + 2 - 6 (i.e. taking genes that belong to lysosomal GO term, but not autophagy GO term)
Mitochondrial_Lysosomal_noAutophagy_genes <- rbind(Mitochondrial_genes, Lysosomal_genes_direct)
Mitochondrial_Lysosomal_noAutophagy_genes <- Mitochondrial_Lysosomal_noAutophagy_genes %>% filter(!(ensembl_ID %in% Lysosomal_AND_Autophagy_genes$ensembl_ID)) %>% unique()
write.csv(Mitochondrial_Lysosomal_noAutophagy_genes, 'Outputs/Mitochondrial_Lysosomal_noAutophagy_genes_12.csv', row.names = FALSE)

## 13 ## 
# 13. nonAutophagy-lysosomal only -> 2 - 6
Lysosomal_noAutophagy_genes <- Lysosomal_genes_direct %>% filter(!(gene_names %in% Lysosomal_AND_Autophagy_genes$gene_names)) %>% unique()
write.csv(Lysosomal_noAutophagy_genes, 'Outputs/Lysosomal_noAutophagy_genes_13.csv', row.names = FALSE)

## 14 ##
# 14. Autophagy OR Lysosomal direct GO -> 2 + 4
Lysosomal_Autophagy_direct_genes <- rbind(Lysosomal_genes_direct, Autophagy_genes_direct) %>% unique()
write.csv(Lysosomal_Autophagy_direct_genes, 'Outputs/Lysosomal_Autophagy_direct_genes_14.csv', row.names = FALSE)

## 15 ##
# 15. Autophagy OR Lysosomal any GO -> 3 + 5
Lysosomal_Autophagy_any_genes <- rbind(Lysosomal_genes_any, Autophagy_genes_any) %>% unique()
write.csv(Lysosomal_Autophagy_any_genes, 'Outputs/Lysosomal_Autophagy_any_genes_15.csv', row.names = FALSE)

