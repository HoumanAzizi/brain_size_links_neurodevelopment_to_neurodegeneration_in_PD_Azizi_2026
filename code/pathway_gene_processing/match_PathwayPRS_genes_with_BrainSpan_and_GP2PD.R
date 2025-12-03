library(tidyr)
library(dplyr)

# if (!require("BiocManager")) install.packages("BiocManager")
# BiocManager::install("biomaRt")
# library(biomaRt)

rm(list = ls())
cat("\014")

setwd("/Users/houmanazizi/Library/Mobile\ Documents/com~apple~CloudDocs/Education/GitHub/Pathway_PD_PRS/")


################### Function to get gene synonyms ###################
# get_gene_synonyms <- function(gene_list) {
#   library(biomaRt) # loading biomaRT and unloading at the end of function because it prevents dplyr to work
#   
#   # Get Ensembl database
#   ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#   
#   # Retrieve synonyms using gene name as filter
#   synonyms <- getBM(attributes = c("hgnc_symbol", "external_synonym"),
#                     filters = "hgnc_symbol",
#                     values = gene_list,
#                     mart = ensembl)
#   
#   # Add original gene name to this dictionary as well
#   original_gene_names <- unique(synonyms$hgnc_symbol)
#   original_gene_names <- data.frame(hgnc_symbol = original_gene_names,
#                                     external_synonym = original_gene_names)
#   synonym_list <- rbind(synonyms, original_gene_names)
#   colnames(synonym_list) <- c('gene_name', 'gene_synonym')
#   
#   # Remove biomaRT and return value
#   rm(ensembl)
#   detach("package:biomaRt", unload = TRUE)
#   return(synonym_list)
# }

get_gene_synonyms <- function(gene_list) {
  library(biomaRt)
  
  # Remove '' values for genes
  gene_list <- gene_list[gene_list != ""]
  
  # Get Ensembl database
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Initialize vector to store ensembl IDs
  ensembl_ids <- c()
  
  # First try to find genes by their HGNC symbols
  direct_matches <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                          filters = "hgnc_symbol",
                          values = gene_list,
                          mart = ensembl)
  
  # Save ensembl IDs for direct matches
  ensembl_ids <- c(ensembl_ids, direct_matches$ensembl_gene_id)
  
  # For genes not found in direct matches, look in synonyms
  unmatched_genes <- setdiff(gene_list, direct_matches$hgnc_symbol)
  
  if(length(unmatched_genes) > 0) {
    # Get genes and their synonyms along with ensembl IDs
    synonym_matches <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "external_synonym"),
                             filters = "external_synonym",
                             values = unmatched_genes,
                             mart = ensembl)
    
    # Add these ensembl IDs
    ensembl_ids <- c(ensembl_ids, synonym_matches$ensembl_gene_id)
  }
  
  # Remove duplicates from ensembl_ids
  ensembl_ids <- unique(ensembl_ids)
  
  # Now get all synonyms for these specific ensembl IDs
  if(length(ensembl_ids) > 0) {
    final_synonyms <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "external_synonym"),
                            filters = "ensembl_gene_id",
                            values = ensembl_ids,
                            mart = ensembl)
    
    # Add the HGNC symbols as synonyms too
    original_names <- data.frame(
      ensembl_gene_id = final_synonyms$ensembl_gene_id[!duplicated(final_synonyms$hgnc_symbol)],
      hgnc_symbol = unique(final_synonyms$hgnc_symbol),
      external_synonym = unique(final_synonyms$hgnc_symbol)
    )
    
    final_synonyms <- rbind(final_synonyms, original_names)
    
    # Create the final dataframe
    result <- data.frame(
      gene_name = final_synonyms$hgnc_symbol,
      gene_synonym = final_synonyms$external_synonym
    )
    
    # Remove any NA values, empty strings, and duplicates
    result <- result[complete.cases(result) & result$gene_synonym != '' & result$gene_name != '', ]
    result <- unique(result)
  } else {
    result <- data.frame(gene_name = character(0), 
                         gene_synonym = character(0))
  }
  
  # Remove biomaRt and return value
  rm(ensembl)
  detach("package:biomaRt", unload = TRUE)
  
  return(result)
}





################### Get file list, get synonyms, save ###################
# Get gene list files
all_files <- list.files('./Outputs/gene_Lists/in_GRCh37', pattern = '*.csv', full.names = TRUE)

# Loop through each file and get the gene synonyms - NEEDS to be done only once
# for (file in all_files) {
#   genes <- read.csv(file)
#   gene_synonyms <- get_gene_synonyms(genes$gene_names)
#   
#   # Create the output file path and name
#   output_file <- file.path("./Outputs/gene_lists_with_synonym", basename(file))
#   output_file <- sub("\\.csv$", "_withSynonyms.csv", output_file)
#   
#   # Write the gene synonyms to the output file
#   write.csv(gene_synonyms, output_file, row.names = FALSE)
#   rm(genes, gene_synonyms)
# }





################### Find which PD genes are in each pathway ###################
all_files <- list.files('./Outputs/gene_Lists/in_GRCh38', pattern = '*.csv', full.names = TRUE)

## FINAL VERSION: Eric Yu's Nalls Nomination + remaining genes directly from GP2
Yu_loci <- read.csv('./Data/PD_GWAS_NallsPrioritizationbyYu_and_GP2genes/EricYu2024_Nalls_PD_Prioritization.csv', header = TRUE) %>% select(-c("Nearest.gene.GWAS", "Nearest.gene.based.on.distance", "Nearest.gene.based.on.TSS"))
GP2_loci_original <- read.csv('./Data/PD_GWAS_NallsPrioritizationbyYu_and_GP2genes/PD_GP2_sig_genes_OriginalTable.csv', header = TRUE) %>% select(Nearest.Gene, Locus.Number, Known.or.Novel)
colnames(Yu_loci) <- c("Nearest.Gene", "Locus_Nalls", "Rank", "Prediction_probability") # NOTE: not really the nearest gene, but named here to work with the rest of the code
colnames(GP2_loci_original) <- c("Nearest.Gene", "Locus_GP2", "Known_or_Novel_inGP2")

 

## Choose which ranks we want - has ranks 1-8
# rank_limit <- 3
# GP2_loci <- GP2_loci %>% filter(Rank <= rank_limit)

## Merge genes and get a final list
rank_limit <- 0.5
Yu_loci <- Yu_loci %>% filter(Prediction_probability > rank_limit)
Yu_loci <- unique(Yu_loci)
GP2_loci_original <- unique(GP2_loci_original)
GP2_loci <- Yu_loci %>% full_join(GP2_loci_original, by = 'Nearest.Gene')
rm(GP2_loci_original, Yu_loci)





## V2: using gene names - for Eric's nomination of GP2 loci
for (file in all_files) {
  # Get the basename of the current file (without the full path and extension)
  file_name <- tools::file_path_sans_ext(basename(file))
  
  # Read the current gene list file with synonyms
  current_gene_list <- read.csv(paste0('./Outputs/gene_Lists_with_Synonym/', file_name, '_withSynonyms.csv'))
  
  # Add a new column to the GP2_loci dataframe with the file name
  GP2_loci[[file_name]] <- NA
  
  # Loop through each row in the current_current_gene_list dataframe
  for (i in 1:nrow(current_gene_list)) {
    # Get the current gene's name
    gene_name <- current_gene_list$gene_synonym[i]
    
    # Check if the current gene's name match any rows in GP2_loci
    matching_rows <- GP2_loci$Nearest.Gene == gene_name
    
    # If there's a match, set the value in the new column to 1
    if (!is.na(all(matching_rows)) & any(matching_rows)) {
      GP2_loci[[file_name]][matching_rows] <- 1
    }
  }
}








### Create a list of genes in brainspan to be included in each pathway
# Get synonyms for all genes in Nalls
GP2_loci_syn <- get_gene_synonyms(GP2_loci$Nearest.Gene)

#GP2_loci_syn <- rbind(GP2_loci_syn, missing_genes_in_GP2)
# Save this separately as well
write.csv(GP2_loci_syn, paste0('./Outputs/gene_lists_with_synonym/GP2_genes_withSynonyms_rank', rank_limit,'.csv'), row.names = FALSE)
# List of genes that we cannot locate:
GP2_loci$Nearest.Gene[(GP2_loci$Nearest.Gene %in% GP2_loci_syn$gene_name == FALSE) & (GP2_loci$Nearest.Gene %in% GP2_loci_syn$gene_synonym == FALSE)]

## Add columns for each gene name and initialize them to 0
# Get all unique gene names from GP2_loci_syn
all_GP2_genes <- unique(c(GP2_loci_syn$gene_name, GP2_loci_syn$gene_synonym))
all_GP2_genes <- all_GP2_genes[!is.na(all_GP2_genes) & all_GP2_genes != ""]
# Create and initialize the result dataframe with gene_list_name column and all gene columns set to 0
pathway_genes_in_GP2_results <- data.frame(gene_list_name = character(), stringsAsFactors = FALSE)
pathway_genes_in_GP2_results <- data.frame(
  gene_list_name = tools::file_path_sans_ext(basename(all_files)),
  matrix(0, 
         nrow = length(all_files), 
         ncol = length(all_GP2_genes) #, dimnames = list(NULL, all_GP2_genes))
))
colnames(pathway_genes_in_GP2_results)[-1] <- all_GP2_genes

# Loop through each file to process pathways -> make values 1 if they exist in the pathway
for (i in 1:length(all_files)) {
  file <- all_files[i]
  file_name <- tools::file_path_sans_ext(basename(file))
  
  # Get genes where the pathway column has value 1
  pathway_genes_initial <- GP2_loci$Nearest.Gene[!is.na(GP2_loci[[file_name]]) & GP2_loci[[file_name]] == 1]
  
  # Initialize vector to store all genes including synonyms
  all_genes <- c()
  
  # For each gene in the pathway
  for (gene in pathway_genes_initial) {
    # Add the original gene
    all_genes <- c(all_genes, gene)
    
    # Find synonyms where gene is in gene_name column
    synonyms1 <- GP2_loci_syn$gene_synonym[GP2_loci_syn$gene_name == gene]
    all_genes <- c(all_genes, synonyms1)
    
    # Find original genes where gene is in synonym column
    original_genes <- GP2_loci_syn$gene_name[GP2_loci_syn$gene_synonym == gene]
    all_genes <- c(all_genes, original_genes)
  }
  
  # Get unique genes
  unique_genes <- unique(all_genes[!is.na(all_genes)])
  
  # Set value to 1 for all genes in this pathway
  pathway_genes_in_GP2_results[i, unique_genes] <- 1
}

# add a separate row with all 1 values which would be the complete Nalls list
row_GP2 <- as.data.frame(t(rep(1, ncol(pathway_genes_in_GP2_results))))
colnames(row_GP2) <- colnames(pathway_genes_in_GP2_results)
row_GP2$gene_list_name[1] <- 'GP2_PD_GWAS_Significant_Genes'
pathway_genes_in_GP2_results <- rbind(pathway_genes_in_GP2_results, row_GP2)

# Save
write.csv(pathway_genes_in_GP2_results, paste0('./Outputs/pathway_genes_in_GP2_PD_GWAS_rank', rank_limit, '.csv'), row.names = FALSE)







################### Find % matching with brainspan ###################
# Read gene lists
target_genes <- read.csv("./Data/Brainspan/brainspan_gene_list.csv")
all_files <- list.files('./Outputs/gene_Lists/in_GRCh38', pattern = '*.csv', full.names = TRUE)

# Calculate total unique genes in GP2_loci
n_genes_GP2_total <- length(unique(na.omit(GP2_loci$Nearest.Gene)))

# Calculate total GP2 genes available in brainspan
GP2_total_genes <- unique(na.omit(GP2_loci$Nearest.Gene))
GP2_total_genes_syn <- read.csv(paste0('./Outputs/gene_lists_with_synonym/GP2_genes_withSynonyms_rank', rank_limit, '.csv'))
GP2_total_in_brainspan <- length(unique(GP2_total_genes_syn$gene_name[GP2_total_genes_syn$gene_synonym %in% target_genes$gene_name]))

# Create an empty data frame to store the summary
summary_table <- data.frame(
  Gene_list = character(),
  N_genes = integer(),
  N_pathway_genes_available_in_brainspan = integer(),
  N_pathway_genes_in_GP2 = integer(),
  N_genes_GP2_total = integer(),
  N_pathway_genes_in_GP2_available_in_brainspan = integer(),
  N_genes_GP2_total_available_in_brainspan = integer()
)

# Loop through each file and process
for (file in all_files) {
  # Get the gene list name
  current_gene_list_name <- substr(basename(file), 1, nchar(basename(file)) - 4)
  
  # Read the gene list files
  current_gene_list <- read.csv(file)
  current_gene_list_syn <- read.csv(paste0('./Outputs/gene_lists_with_synonym/', current_gene_list_name, '_withSynonyms.csv'))
  
  # Get the number of genes
  n_genes <- nrow(current_gene_list)
  
  # Count the number of genes available in the target list -> match based on their synonym to capture all
  available_genes <- unique(current_gene_list_syn$gene_name[current_gene_list_syn$gene_synonym %in% target_genes$gene_name])
  n_available <- length(available_genes)
  
  # Count genes in GP2
  n_genes_in_GP2 <- 0  # default value
  n_genes_in_GP2_in_brainspan <- 0  # default value
  if(current_gene_list_name %in% colnames(GP2_loci)) {
    # Get GP2 genes for this pathway
    GP2_pathway_genes <- unique(na.omit(GP2_loci$Nearest.Gene[GP2_loci[[current_gene_list_name]] == 1]))
    n_genes_in_GP2 <- length(GP2_pathway_genes)
    
    # Get synonyms for these GP2 pathway genes
    GP2_pathway_genes <- get_gene_synonyms(GP2_pathway_genes)
    
    # Find which of these GP2 pathway genes are in brainspan using current_gene_list_syn #### EDIT THISSSSS??
    GP2_genes_in_brainspan <- unique(GP2_pathway_genes$gene_name[GP2_pathway_genes$gene_synonym %in% target_genes$gene_name])
    n_genes_in_GP2_in_brainspan <- length(GP2_genes_in_brainspan)
  }
  
  # Add the row to the summary table
  summary_table <- rbind(summary_table, data.frame(
    Gene_list = current_gene_list_name,
    N_pathway_genes = n_genes,
    N_pathway_genes_available_in_brainspan = n_available,
    Percentage_pathway_genes_available_in_brainspan = (n_available*100)/n_genes,
    N_pathway_genes_in_GP2 = n_genes_in_GP2,
    N_genes_GP2_total = n_genes_GP2_total,
    Percentage_of_pathway_genes_in_GP2 = (n_genes_in_GP2*100)/n_genes_GP2_total,
    N_pathway_genes_in_GP2_available_in_brainspan = n_genes_in_GP2_in_brainspan,
    N_genes_GP2_total_available_in_brainspan = GP2_total_in_brainspan,
    Percentage_pathway_gene_in_GP2_available_in_brainspan = (n_genes_in_GP2_in_brainspan*100)/GP2_total_in_brainspan
  ))
}
write.csv(summary_table, paste0('./Outputs/pathway_Genes_GP2_Percentage_present_in_BrainSpan_rank', rank_limit, '.csv'), row.names = FALSE)








################### Create supplementary table to match GP2 genes to pathways ###################
colnames(GP2_loci)

GP2_loci <- GP2_loci %>% select(Nearest.Gene, Locus_Nalls, Locus_GP2, 
                                Known_or_Novel_inGP2, Prediction_probability, Rank,
                                Mitochondrial_genes_1, Autophagy_direct_genes_4, Lysosomal_direct_genes_2)
colnames(GP2_loci) <- c('Gene name', 
                        'Locus number in Nalls', 'Locus number in GP2', 
                        'Known or Novel in GP2', 
                        'Prediction probability', 'Prioritization Rank Yu2024',
                        'Mitochondrial pathway', 'Autophagy pathway', 'Lysosomal pathway')

GP2_loci <- GP2_loci %>%
  mutate(across(c(`Mitochondrial pathway`, `Autophagy pathway`, `Lysosomal pathway`), 
                ~ ifelse(is.na(.), "No", ifelse(. == 1, "Yes", .))))
# GP2_loci <- GP2_loci %>% mutate(Chromosome = gsub("chr", "", Chromosome))
# GP2_loci$Chromosome <- as.numeric(GP2_loci$Chromosome)

write.csv(GP2_loci, paste0('./Outputs/PD_GP2_genes_pathway_annotated_rank', rank_limit, '_noPosition.csv'), row.names = FALSE)







