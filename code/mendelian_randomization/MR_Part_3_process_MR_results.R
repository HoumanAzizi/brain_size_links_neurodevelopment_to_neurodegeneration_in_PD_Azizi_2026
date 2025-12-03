#################### Basics ####################
rm(list = ls())
gc()
cat("\014")

library(ggplot2)
library(dplyr)
#library(TwoSampleMR)

setwd('/Users/HoumanAzizi/Library/Mobile\ Documents/com~apple~CloudDocs/Education/GitHub/Mendelian_Randomization_UKB/')




#################### Functions ####################
calculate_mr_power <- function(r2) {
  # Fixed parameters
  OR <- 1.2
  n_cases <- 63555 + 17700
  n_controls <- 1746386
  alpha <- 0.05
  
  # Calculate total sample size
  N <- n_cases + n_controls
  
  # Calculate the proportion of cases (K in the mRnd code)
  K <- n_cases / N
  
  # Calculate the regression coefficient b_MR according to mRnd formula
  b_MR <- K * (OR / (1 + K * (OR - 1)) - 1)
  
  # Calculate the variance of the MR estimator
  v_MR <- (K * (1 - K) - b_MR^2) / (N * r2)
  
  # Calculate the non-centrality parameter
  ncp <- b_MR^2 / v_MR
  
  # Calculate power using the non-central chi-squared distribution
  crit_value <- qchisq(1 - alpha, df = 1)
  power <- 1 - pchisq(crit_value, df = 1, ncp = ncp)
  
  # Calculate F-statistic (instrument strength)
  f_value <- 1 + N * r2 / (1 - r2)
  
  # Return the power
  return(power)
}












#################### Read MR results ####################
DKT_SA <- read.csv('Outputs/result_summary_allMeasures/DKTatlas_SA_MR_Results_complete.csv', header = TRUE)[, -1]
DKT_CT <- read.csv('Outputs/result_summary_allMeasures/DKTatlas_CT_MR_Results_complete.csv', header = TRUE)[, -1]
QSM <- read.csv('Outputs/result_summary_allMeasures/QSM_MR_Results_complete.csv', header = TRUE)[, -1]

Warrier_GWAS <- read.csv('Outputs/result_summary_allMeasures/Warrier_GWAS_MR_Results_complete.csv', header = TRUE)[, -1]
Warrier_GWAS[1,1] <- 'SA_meta_v1'
Warrier_GWAS_2 <- read.csv('Outputs/result_summary_allMeasures/Warrier_GWAS_MR_Results_SA_fixed_v3.csv', header = TRUE)[, -1]
Warrier_GWAS_2 <- Warrier_GWAS_2[1, ]
Warrier_GWAS_2[1,1] <- 'SA_meta_v3'
Warrier_GWAS <- rbind(Warrier_GWAS, Warrier_GWAS_2)
rm(Warrier_GWAS_2)

Subcortical_aseg <- read.csv('Outputs/result_summary_allMeasures/subcortical_volume_aseg_MR_Results_complete.csv', header = TRUE)[, -1]
Subcortical_FIRST <- read.csv('Outputs/result_summary_allMeasures/subcortical_volume_MR_Results_complete.csv', header = TRUE)[1, -1]
WM_FA <- read.csv('Outputs/result_summary_allMeasures/WM_FA_MR_Results_complete.csv', header = TRUE)[, -1]
WM_MD <- read.csv('Outputs/result_summary_allMeasures/WM_MD_MR_Results_complete.csv', header = TRUE)[, -1]
global_WM <- read.csv('Outputs/result_summary_allMeasures/Zhao_GWAS_MR_Results_complete.csv', header = TRUE)[, -1]

reverse_MR <- read.csv('Outputs/result_summary_allMeasures/reverseMR_PD_vsBrain_MR_Results_complete.csv', header = TRUE)[, -1]
reverse_MR2 <- read.csv('Outputs/result_summary_allMeasures/reverseMR_PD_vsBrain_MR_Results_complete_QSM_Zhao.csv', header = TRUE)[, -1]







#################### Calculate FDR ####################
# List of dataframes to process
df_names <- c("DKT_SA", "DKT_CT", "QSM", "Subcortical_aseg", "WM_FA", "WM_MD")

for (df_name in df_names) {
  # Get the dataframe from the environment
  df <- get(df_name)
  
  # Find columns ending with '_pval'
  pval_cols <- grep("_pval$", names(df), value = TRUE)
  
  # Calculate N (number of rows where MR_error is NA)
  N <- sum(is.na(df$MR_error))
  
  # For each p-value column, calculate FDR-corrected p-values
  for (col in pval_cols) {
    # Create new column name by replacing '_pval' with '_pval_fdr'
    new_col <- gsub("_pval$", "_pval_fdr", col)
    
    # Calculate FDR-corrected p-values
    df[[new_col]] <- p.adjust(df[[col]], method = "fdr", n = N)
  }
  
  # Assign the modified dataframe back to the environment
  assign(df_name, df, envir = .GlobalEnv)
  
  # Print confirmation
  cat(sprintf("Processed %s: Added FDR-corrected p-values\n", df_name))
  rm(df, col, df_name, N, new_col, pval_cols)
}





#################### Merge MR results ####################
# Create MR_regionwise by merging regional measures
MR_regionwise <- rbind(DKT_SA, DKT_CT, Subcortical_aseg, QSM, WM_FA, WM_MD)

# Create MR_global by merging global measures
MR_global <- rbind(Warrier_GWAS, Subcortical_FIRST, global_WM)

# Create MR_reverse by merging reverse MR results
MR_reverse <- rbind(reverse_MR, reverse_MR2)

# Remove original dataframes from environment
rm(DKT_SA, DKT_CT, QSM, WM_FA, WM_MD, Subcortical_aseg, 
   Warrier_GWAS, Subcortical_FIRST, global_WM,
   reverse_MR, reverse_MR2,
   df_names)







#################### Sort Results ####################
# Bring fdr corrected values to front
MR_regionwise <- MR_regionwise %>%
  relocate(
    "Inverse_variance_weighted_pval_fdr",
    "Weighted_median_pval_fdr",
    "MR_Egger_pval_fdr",
    "Simple_mode_pval_fdr",
    "Weighted_mode_pval_fdr",
    .after = "MR_error"
  )

# Create significance column
MR_regionwise <- MR_regionwise %>% mutate(IVW_sig = case_when(Inverse_variance_weighted_pval_fdr <= 0.05 ~ 1,
                                                              TRUE ~ NA)) %>% relocate('IVW_sig', .after = 'MR_error')
MR_global <- MR_global %>% mutate(IVW_sig = case_when(Inverse_variance_weighted_pval <= 0.05 ~ 1,
                                                      TRUE ~ NA)) %>% relocate('IVW_sig', .after = 'MR_error')
MR_reverse <- MR_reverse %>% mutate(IVW_sig = case_when(Inverse_variance_weighted_pval <= 0.05 ~ 1,
                                                        TRUE ~ NA)) %>% relocate('IVW_sig', .after = 'MR_error')

# Reorder from low IVW p to high
MR_regionwise <- MR_regionwise %>% arrange(Inverse_variance_weighted_pval_fdr)
MR_global <- MR_global %>% arrange(Inverse_variance_weighted_pval)
MR_reverse <- MR_reverse %>% arrange(Inverse_variance_weighted_pval)





#################### Add Power ####################
# Add power column to each dataframe
MR_regionwise$Power <- sapply(MR_regionwise$r2Exp_calculated, calculate_mr_power) * 100
MR_global$Power <- sapply(MR_global$r2Exp_calculated, calculate_mr_power) * 100
MR_reverse$Power <- sapply(MR_reverse$r2Exp_calculated, calculate_mr_power) * 100

# Round to 2 decimal places for better readability (optional)
MR_regionwise$Power <- round(MR_regionwise$Power, 2)
MR_global$Power <- round(MR_global$Power, 2)
MR_reverse$Power <- round(MR_reverse$Power, 2)









#################### Merge significant values ####################
all_columns <- union(names(MR_global), names(MR_regionwise))

# Add missing columns to MR_global with NA values
for (col in all_columns) {
  if (!col %in% names(MR_global)) {
    MR_global[[col]] <- NA
  }
}

MR_main <- bind_rows(MR_global, MR_regionwise)
MR_main <- MR_main %>% filter(IVW_sig == 1) %>% relocate("Inverse_variance_weighted_pval_fdr",
                                                         "Weighted_median_pval_fdr",
                                                         "MR_Egger_pval_fdr",
                                                         "Simple_mode_pval_fdr",
                                                         "Weighted_mode_pval_fdr",
                                                         .after = "MR_error")









#################### Keep important columns only ####################
main_cols <- c("GWAS",
                "Inverse_variance_weighted_nsnp",
                "Inverse_variance_weighted_pval",
                "Inverse_variance_weighted_pval_fdr",
                "Weighted_median_pval",
                "Weighted_median_pval_fdr",
                "MR_Egger_pval",
                "MR_Egger_pval_fdr",
                "Simple_mode_pval",
                "Simple_mode_pval_fdr",
                "Weighted_mode_pval",
                "Weighted_mode_pval_fdr",
                "Pleiotropy_pval",
                "Pleiotropy_pval_fdr",
                "Heterogeneity_Inverse_variance_weighted_Q_pval",
                "Heterogeneity_Inverse_variance_weighted_Q_pval_fdr",
                'r2Exp_calculated',
                'Fvalue_from_r2Exp',
                'Power')
MR_main_summary <- MR_main %>% select(all_of(main_cols))













#################### Create a more understandable summary ####################
MR_total <- MR_main_summary %>% 
  mutate(IVW = case_when( is.na(Inverse_variance_weighted_pval_fdr) & Inverse_variance_weighted_pval <= 0.05 ~ 1,
                          !is.na(Inverse_variance_weighted_pval_fdr) & Inverse_variance_weighted_pval_fdr <= 0.05 ~ 1,
                          TRUE ~ NA ),
         Weighted_median = case_when( is.na(Weighted_median_pval_fdr) & Weighted_median_pval <= 0.05 ~ 1,
                                      !is.na(Weighted_median_pval_fdr) & Weighted_median_pval_fdr <= 0.05 ~ 1,
                                      TRUE ~ NA ),
         MR_Egger = case_when( is.na(MR_Egger_pval_fdr) & MR_Egger_pval <= 0.05 ~ 1,
                               !is.na(MR_Egger_pval_fdr) & MR_Egger_pval_fdr <= 0.05 ~ 1,
                               TRUE ~ NA ),
         Simple_mode = case_when( is.na(Simple_mode_pval_fdr) & Simple_mode_pval <= 0.05 ~ 1,
                                  !is.na(Simple_mode_pval_fdr) & Simple_mode_pval_fdr <= 0.05 ~ 1,
                                  TRUE ~ NA ),
         Weighted_mode = case_when( is.na(Weighted_mode_pval_fdr) & Weighted_mode_pval <= 0.05 ~ 1,
                                    !is.na(Weighted_mode_pval_fdr) & Weighted_mode_pval_fdr <= 0.05 ~ 1,
                                    TRUE ~ NA ),
         Pleiotropy = case_when( is.na(Pleiotropy_pval_fdr) & Pleiotropy_pval <= 0.05 ~ 1,
                                 !is.na(Pleiotropy_pval_fdr) & Pleiotropy_pval_fdr <= 0.05 ~ 1,
                                 TRUE ~ NA ),
         Heterogeneity_IVW = case_when( is.na(Heterogeneity_Inverse_variance_weighted_Q_pval_fdr) & Heterogeneity_Inverse_variance_weighted_Q_pval <= 0.05 ~ 1,
                                        !is.na(Heterogeneity_Inverse_variance_weighted_Q_pval_fdr) & Heterogeneity_Inverse_variance_weighted_Q_pval_fdr <= 0.05 ~ 1,
                                        TRUE ~ NA )) %>% 
  select(GWAS, IVW, Weighted_median, MR_Egger, Simple_mode, Weighted_mode, Pleiotropy, Heterogeneity_IVW) %>% 
  rowwise() %>% mutate(total_MR_sig = sum(IVW, Weighted_median, MR_Egger, Simple_mode, Weighted_mode, na.rm = TRUE))






#### Save
write.csv(MR_main, 'Outputs/result_summary_allMeasures/summary_tables/sig_Brain_vs_PD_MR_results.csv', row.names = FALSE)
write.csv(MR_main_summary, 'Outputs/result_summary_allMeasures/summary_tables/sig_Brain_vs_PD_MR_results_summary.csv', row.names = FALSE)
write.csv(MR_total, 'Outputs/result_summary_allMeasures/summary_tables/sig_Brain_vs_PD_MR_results_summary2.csv', row.names = FALSE)

write.csv(MR_regionwise, 'Outputs/result_summary_allMeasures/summary_tables/Regionwise_Brain_vs_PD_MR_results_complete.csv', row.names = FALSE)
write.csv(MR_global, 'Outputs/result_summary_allMeasures/summary_tables/Global_Brain_vs_PD_MR_results_complete.csv', row.names = FALSE)
write.csv(MR_reverse, 'Outputs/result_summary_allMeasures/summary_tables/Reverse_Brain_vs_PD_MR_results_complete.csv', row.names = FALSE)
