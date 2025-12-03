# Script to create supplementary tables
library(readr)       # For reading CSVs efficiently
library(openxlsx)    # For Excel file manipulation
library(dplyr)       # For data manipulation
library(here)        # For path management
library(readxl)

rm(list = ls())
cat("\014")

# Set working directory
setwd("/Users/houmanazizi/Library/Mobile\ Documents/com~apple~CloudDocs/Education/GitHub/")






###################### ORIGINAL DATA ######################
## Read sample data to get the region orders
sample_tables <- read_excel(
  'UKB_Analysis/Outputs/Paper1_Tables/sample_Supplementary_Tables_old.xlsx', 
  sheet = 'Table S1',
  skip = 1,           # Skip the title row
  col_names = TRUE    # Use the first row after skipping as column names
) %>% filter(!is.na(beta))
DKT_region_order <- sample_tables$Region

## Create mapping between final region names and the region names in CSV file
region_mapping <- c(
  "caudalmiddlefrontal_L" = "L Caudal middle frontal",
  "entorhinal_L" = "L Entorhinal",
  "postcentral_L" = "L Postcentral",
  "parstriangularis_L" = "L Pars triangularis",
  "supramarginal_L" = "L Supramarginal",
  "insula_L" = "L Insula",
  "lateralorbitofrontal_L" = "L Lateral orbitofrontal",
  "parsorbitalis_L" = "L Pars orbitalis",
  "middletemporal_L" = "L Middle temporal",
  "pericalcarine_L" = "L Pericalcarine",
  "parahippocampal_L" = "L Parahippocampal",
  "paracentral_L" = "L Paracentral",
  "medialorbitofrontal_L" = "L Medial orbitofrontal",
  "cuneus_L" = "L Cuneus",
  "inferiortemporal_L" = "L Inferior temporal",
  "rostralmiddlefrontal_L" = "L Rostral middle frontal",
  "rostralanteriorcingulate_L" = "L Rostral anterior cingulate",
  "isthmuscingulate_L" = "L Isthmus cingulate",
  "lateraloccipital_L" = "L Lateral occipital",
  "lingual_L" = "L Lingual",
  "superiorparietal_L" = "L Superior parietal",
  "parsopercularis_L" = "L Pars opercularis",
  "fusiform_L" = "L Fusiform",
  "caudalanteriorcingulate_L" = "L Caudal anterior cingulate",
  "superiorfrontal_L" = "L Superior frontal",
  "precuneus_L" = "L Precuneus",
  "transversetemporal_L" = "L Transverse temporal",
  "precentral_L" = "L Precentral",
  "inferiorparietal_L" = "L Inferior parietal",
  "posteriorcingulate_L" = "L Posterior cingulate",
  "superiortemporal_L" = "L Superior temporal",
  "caudalmiddlefrontal_R" = "R Caudal middle frontal",
  "entorhinal_R" = "R Entorhinal",
  "postcentral_R" = "R Postcentral",
  "parstriangularis_R" = "R Pars triangularis",
  "supramarginal_R" = "R Supramarginal",
  "insula_R" = "R Insula",
  "lateralorbitofrontal_R" = "R Lateral orbitofrontal",
  "parsorbitalis_R" = "R Pars orbitalis",
  "middletemporal_R" = "R Middle temporal",
  "pericalcarine_R" = "R Pericalcarine",
  "parahippocampal_R" = "R Parahippocampal",
  "paracentral_R" = "R Paracentral",
  "medialorbitofrontal_R" = "R Medial orbitofrontal",
  "cuneus_R" = "R Cuneus",
  "inferiortemporal_R" = "R Inferior temporal",
  "rostralmiddlefrontal_R" = "R Rostral middle frontal",
  "rostralanteriorcingulate_R" = "R Rostral anterior cingulate",
  "isthmuscingulate_R" = "R Isthmus cingulate",
  "lateraloccipital_R" = "R Lateral occipital",
  "lingual_R" = "R Lingual",
  "superiorparietal_R" = "R Superior parietal",
  "parsopercularis_R" = "R Pars opercularis",
  "fusiform_R" = "R Fusiform",
  "caudalanteriorcingulate_R" = "R Caudal anterior cingulate",
  "superiorfrontal_R" = "R Superior frontal",
  "precuneus_R" = "R Precuneus",
  "transversetemporal_R" = "R Transverse temporal",
  "precentral_R" = "R Precentral",
  "inferiorparietal_R" = "R Inferior parietal",
  "posteriorcingulate_R" = "R Posterior cingulate",
  "superiortemporal_R" = "R Superior temporal"
)

## Read tract categories for WM
WM_categories <- read.csv('UKB_PLS/Data/White_Matter/ORG_Tract_Names_and_Categories_forPLS.csv', header = FALSE) %>% select(-V2)
colnames(WM_categories) <- c('Tract_file_name', 'Category')
WM_categories$Tract_file_name <- gsub("\\.", "-", WM_categories$Tract_file_name)









###################### FUNCTIONS ######################
### Function to format the decimals and scientific notations
## Function for SA and CT
format_prs_results_SA_CT <- function(df, decimals = list(t_value = 3, r2 = 3, r2_adj = 3), 
                               scientific = list(effect_size = 2, p_value = 2, p_FDR = 2, SE = 2)) {
  
  ## How to use?
  # # Basic
  # PRS_SA_formatted <- format_prs_results(PRS_SA)
  # # Customized
  # PRS_SA_formatted <- format_prs_results(PRS_SA, decimals = list(t_value = 3, r2 = 3, r2_adj = 3), scientific = list(effect_size = 2, p_value = 2, p_FDR = 2, SE = 2))
  
  
  ### Fix the SA/CT data
  # Add final region names and reorder
  df$Region_names_final <- region_mapping[df$Region]
  df <- df %>% 
    mutate(Region_names_final = factor(Region_names_final, levels = DKT_region_order)) %>% 
    arrange(Region_names_final)
  
  # Select required columns & rename them for final use
  df <- df %>% select(Region_names_final, effect_size, t_value, p_value, p_FDR, r2, r2_adj, SE)
  colnames(df)[1] <- 'Region'
  
  
  ### Reformat
  formatted_df <- df
  # Format decimal values
  for (col_name in names(decimals)) {
    if (col_name %in% colnames(formatted_df)) {
      formatted_df[[col_name]] <- formatC(df[[col_name]], format = "f", digits = decimals[[col_name]])
    }
  }
  # Format scientific values
  for (col_name in names(scientific)) {
    if (col_name %in% colnames(formatted_df)) {
      formatted_df[[col_name]] <- formatC(df[[col_name]], format = "e", digits = scientific[[col_name]])
    }
  }
  
  return(formatted_df)
}


## Function for WM
format_prs_results_WM <- function(df, decimals = list(t_value = 3, r2 = 3, r2_adj = 3), 
                                     scientific = list(effect_size = 2, p_value = 2, p_FDR = 2, SE = 2)) {
  
  ## How to use?
  # # Basic
  # PRS_SA_formatted <- format_prs_results(PRS_SA)
  # # Customized
  # PRS_SA_formatted <- format_prs_results(PRS_SA, decimals = list(t_value = 3, r2 = 3, r2_adj = 3), scientific = list(effect_size = 2, p_value = 2, p_FDR = 2, SE = 2))
  
  
  ### Fix the WM data
  # Add categories
  df$Tract_file_name <- gsub("\\.", "-", df$Tract_file_name)
  df$Tract_file_name <- gsub("\\&", "-", df$Tract_file_name)
  df <- df %>% left_join(WM_categories, by = 'Tract_file_name')
  # Rename the tracts
  laterality <- ifelse(grepl("_R$", df$Tract), "R", ifelse(grepl("_L$", df$Tract), "L", "")) # extract the R or L from the end if present
  df$Tract <- gsub("_(R|L)$", "", df$Tract) # remove the _R or _L from the end where present
  df$Tract <- gsub("_", " ", df$Tract) # replace all remaining underscores with spaces
  df$Tract <- ifelse(laterality != "", paste(laterality, df$Tract), df$Tract) # add the laterality to the beginning where applicable
  # Reorder
  df <- df %>% arrange(Category, Tract_file_name)
  # Select required columns & rename them for final use
  df <- df %>% select(Tract, Category, effect_size, t_value, p_value, p_FDR, r2, r2_adj, SE)
  
  
  ### Reformat
  formatted_df <- df
  # Format decimal values
  for (col_name in names(decimals)) {
    if (col_name %in% colnames(formatted_df)) {
      formatted_df[[col_name]] <- formatC(df[[col_name]], format = "f", digits = decimals[[col_name]])
    }
  }
  # Format scientific values
  for (col_name in names(scientific)) {
    if (col_name %in% colnames(formatted_df)) {
      formatted_df[[col_name]] <- formatC(df[[col_name]], format = "e", digits = scientific[[col_name]])
    }
  }
  
  return(formatted_df)
}


## Function for Subcortical Volume
format_prs_results_SV <- function(df, decimals = list(t_value = 3, r2 = 3, r2_adj = 3), 
                                  scientific = list(effect_size = 2, p_value = 2, p_FDR = 2, SE = 2)) {
  
  ## How to use?
  # # Basic
  # PRS_SA_formatted <- format_prs_results(PRS_SA)
  # # Customized
  # PRS_SA_formatted <- format_prs_results(PRS_SA, decimals = list(t_value = 3, r2 = 3, r2_adj = 3), scientific = list(effect_size = 2, p_value = 2, p_FDR = 2, SE = 2))
  
  
  ### Fix the SV data
  # Reorder
  df <- df %>% arrange(Region)
  # Rename the tracts
  laterality <- ifelse(grepl("_R$", df$Region), "R", ifelse(grepl("_L$", df$Region), "L", "")) # extract the R or L from the end if present
  df$Region <- gsub("_(R|L)$", "", df$Region) # remove the _R or _L from the end where present
  df$Region <- ifelse(laterality != "", paste(laterality, df$Region), df$Region) # add the laterality to the beginning where applicable
  # Select required columns & rename them for final use
  df <- df %>% select(Region, effect_size, t_value, p_value, p_FDR, r2, r2_adj, SE)
  
  
  ### Reformat
  formatted_df <- df
  # Format decimal values
  for (col_name in names(decimals)) {
    if (col_name %in% colnames(formatted_df)) {
      formatted_df[[col_name]] <- formatC(df[[col_name]], format = "f", digits = decimals[[col_name]])
    }
  }
  # Format scientific values
  for (col_name in names(scientific)) {
    if (col_name %in% colnames(formatted_df)) {
      formatted_df[[col_name]] <- formatC(df[[col_name]], format = "e", digits = scientific[[col_name]])
    }
  }
  
  return(formatted_df)
}


## Function for Subcortical Volume
format_prs_results_SWI <- function(df, decimals = list(t_value = 3, r2 = 3, r2_adj = 3), 
                                   scientific = list(effect_size = 2, p_value = 2, p_FDR = 2, SE = 2)) {
  
  ## How to use?
  # # Basic
  # PRS_SA_formatted <- format_prs_results(PRS_SA)
  # # Customized
  # PRS_SA_formatted <- format_prs_results(PRS_SA, decimals = list(t_value = 3, r2 = 3, r2_adj = 3), scientific = list(effect_size = 2, p_value = 2, p_FDR = 2, SE = 2))
  
  
  ### Fix the SWI data
  # Reorder
  df <- df %>% arrange(Region)
  # Rename the tracts
  df$Region <- gsub("Median_MagneticSusceptibility_", "", df$Region)
  laterality <- ifelse(grepl("_R$", df$Region), "R", ifelse(grepl("_L$", df$Region), "L", "")) # extract the R or L from the end if present
  df$Region <- gsub("_(R|L)$", "", df$Region) # remove the _R or _L from the end where present
  df$Region <- ifelse(laterality != "", paste(laterality, df$Region), df$Region) # add the laterality to the beginning where applicable
  # Select required columns & rename them for final use
  df <- df %>% select(Region, effect_size, t_value, p_value, p_FDR, r2, r2_adj, SE)
  
  
  ### Reformat
  formatted_df <- df
  # Format decimal values
  for (col_name in names(decimals)) {
    if (col_name %in% colnames(formatted_df)) {
      formatted_df[[col_name]] <- formatC(df[[col_name]], format = "f", digits = decimals[[col_name]])
    }
  }
  # Format scientific values
  for (col_name in names(scientific)) {
    if (col_name %in% colnames(formatted_df)) {
      formatted_df[[col_name]] <- formatC(df[[col_name]], format = "e", digits = scientific[[col_name]])
    }
  }
  
  return(formatted_df)
}






###################### DEFINE GROPS TO LOOP THROUGH ######################
sex_groups <- c(
  "with_all_pathwayPRS_withMotion_all",
  "with_all_pathwayPRS_withMotion_male",
  "with_all_pathwayPRS_withMotion_female"
)

prs_groups <- c(
  "PD_PRScs_GP2",
  "PD_PRScs_GP2_Mitochondrial_Lysosomal_Autophagy_direct_genes_9",
  "PD_PRScs_GP2_Mitochondrial_Lysosomal_Autophagy_direct_genes_9_exclude",
  "PD_PRScs_GP2_Autophagy_direct_genes_4",
  "PD_PRScs_GP2_Lysosomal_direct_genes_2",
  "PD_PRScs_GP2_Mitochondrial_genes_1"
)

measure_types <- c("SA_DKT_", "CT_DKT_", "FA_Tracts_", "MD_Tracts_", "SubcorticaVolume_", "SWI_")

# Map measure types to abbreviations for sheet names
measure_abbreviations <- c(
  "SA_DKT_" = "SA",
  "CT_DKT_" = "CT", 
  "FA_Tracts_" = "FA",
  "MD_Tracts_" = "MD",
  "SubcorticaVolume_" = "SV",
  "SWI_" = "SWI"
)

# Map sex group folder names to simple labels for titles
sex_labels <- c(
  "with_all_pathwayPRS_withMotion_all" = "All",
  "with_all_pathwayPRS_withMotion_male" = "Male",
  "with_all_pathwayPRS_withMotion_female" = "Female"
)

# Map PRS group names to more readable names for titles
prs_labels <- c(
  "PD_PRScs_GP2" = "all genes PD-PRS",
  "PD_PRScs_GP2_Mitochondrial_Lysosomal_Autophagy_direct_genes_9" = "combined mito-ALP PD-PRS",
  "PD_PRScs_GP2_Mitochondrial_Lysosomal_Autophagy_direct_genes_9_exclude" = "non-mito-ALP PD-PRS",
  "PD_PRScs_GP2_Autophagy_direct_genes_4" = "autophagy pathway PD-PRS",
  "PD_PRScs_GP2_Lysosomal_direct_genes_2" = "lysosomal pathway PD-PRS",
  "PD_PRScs_GP2_Mitochondrial_genes_1" = "mitochondrial pathway PD-PRS"
)

# Map measure types to readable names
measure_labels <- c(
  "SA_DKT_" = "cortical surface area",
  "CT_DKT_" = "cortical thickness",
  "FA_Tracts_" = "white matter fractional anisotropy (FA)",
  "MD_Tracts_" = "white matter mean diffusivity (MD)",
  "SubcorticaVolume_" = "subcortical volume",
  "SWI_" = "SWI"
)










###################### CREATE TABLES ######################
# Create a new workbook for all results
wb <- createWorkbook()
sheet_counter <- 1


# Loop through each sex group
for (sex_group in sex_groups) {
  sex_label <- sex_labels[sex_group]
  
  # Loop through each PRS group
  for (prs_group in prs_groups) {
    prs_label <- prs_labels[prs_group]
    
    # Loop through each measure type
    for (measure_type in measure_types) {
      measure_label <- measure_labels[measure_type]
      
      # Construct the file path
      file_path <- file.path("UKB_Analysis/Outputs/PD_PRS_Regression_Results/with_all_pathwayPRS_Brain", 
                             sex_group, prs_group, 
                             paste0(measure_type, "vs_", prs_group, "_Regression_Results.csv"))
      
      
      # Read the CSV file
      PRS_data <- read.csv(file_path)
      
      
      # Format the values
      if (measure_type == 'SA_DKT_' || measure_type == 'CT_DKT_') {
        PRS_formatted <- format_prs_results_SA_CT(PRS_data)
      } else if (measure_type == 'FA_Tracts_' || measure_type == 'MD_Tracts_') {
        PRS_formatted <- format_prs_results_WM(PRS_data)
      } else if (measure_type == 'SubcorticaVolume_') {
        PRS_formatted <- format_prs_results_SV(PRS_data)
      } else if (measure_type == 'SWI_') {
        PRS_formatted <- format_prs_results_SWI(PRS_data)
      }
      
      # Set titles and sheet name
      current_table_title <- paste0(prs_label, " vs ", measure_label, " (", sex_label, ")")
      current_sheet_name <- paste0(sheet_counter, "_", substr(prs_label, 1, nchar(prs_label) - 7), "_", measure_abbreviations[measure_type], "_", sex_label)
      
      # Ensure sheet name is valid (within Excel's 31 character limit)
      if (nchar(current_sheet_name) > 31) {
        current_sheet_name <- substr(current_sheet_name, 1, 31)
      }
      
      # Add a worksheet
      addWorksheet(wb, current_sheet_name)
      
      # Write title and values
      writeData(wb, current_sheet_name, current_table_title, startRow = 1, startCol = 1)
      writeData(wb, current_sheet_name, PRS_formatted, startRow = 2, startCol = 1, 
                headerStyle = createStyle(textDecoration = "bold"))
      
      # Set column width
      setColWidths(wb, current_sheet_name, cols = 1, widths = "auto")
      
      # Make significant p_FDR values bold
      p_FDR_col_index <- which(names(PRS_formatted) == "p_FDR")
      sig_rows <- which(as.numeric(PRS_formatted$p_FDR) <= 0.05)
      
      # Apply bold formatting to significant p-values
      sig_style <- createStyle(textDecoration = "bold")
      if (length(sig_rows) > 0) {
        addStyle(wb, current_sheet_name, style = sig_style, 
                 rows = sig_rows + 2, # Adding 2 to account for header row and title
                 cols = p_FDR_col_index, 
                 gridExpand = TRUE)
      }
      
      # Increment sheet counter
      sheet_counter <- sheet_counter + 1
      
      cat("Processed:", file_path, "\n")
    }
  }
}


# Save workbook
saveWorkbook(wb, "UKB_Analysis/Outputs/Paper1_Tables/all_Regression_results_table_supplementary_newNames.xlsx", overwrite = TRUE)




