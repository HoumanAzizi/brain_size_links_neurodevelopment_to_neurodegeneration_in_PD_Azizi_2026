### Notes for final result
# r values are calculated as mean
# r2 values are calculate as sum of r2 across SNPs
###

rm(list = ls())
cat("\014")

#install.packages('remotes')
Sys.setenv(GITHUB_PAT = "token") # Personal github token
#remotes::install_github("MRCIEU/TwoSampleMR")
#remotes::install_github("MRCIEU/genetics.binaRies")
library(TwoSampleMR)
library(MRPRESSO)
library(R.utils)
library(ggplot2)
library(ieugwasr)
library(genetics.binaRies)
library(Cairo)
library(markdown)
library(dplyr)
library(data.table)

# get OPENGWAS token from https://cran.r-project.org/web/packages/ieugwasr/vignettes/guide.html
# usethis::edit_r_environ() -> add OPENGWAS_JWT=your_token_here
ieugwasr::get_opengwas_jwt()
ieugwasr::user() # Making sure the token works

setwd('/Users/HoumanAzizi/Library/Mobile\ Documents/com~apple~CloudDocs/Education/GitHub/Mendelian_Randomization_UKB/')



###################### Helper functions ######################
extract_mr_results_table <- function(file_content) {
  # Find the start of the Results section
  start_line <- which(grepl("### Results from two sample MR:", file_content))
  
  # Extract table lines (including header, but excluding separator)
  table_lines <- file_content[(start_line + 2):(start_line + 9)]  # Increased to capture all rows
  table_lines <- table_lines[table_lines != ""]  # Remove empty lines
  
  # Process each line
  data <- lapply(table_lines, function(line) {
    # Split the line by '|' and remove empty entries
    parts <- strsplit(line, "\\|")[[1]]
    parts <- parts[parts != ""]
    
    # Trim whitespace from each part
    parts <- trimws(parts)
    
    return(parts)
  })
  
  # Extract column names from the first row
  col_names <- data[[1]]
  
  # Remove the header and separator rows
  data <- data[-(1:2)]
  
  # Convert to dataframe
  df <- as.data.frame(do.call(rbind, data), stringsAsFactors = FALSE)
  colnames(df) <- col_names
  
  # Convert numeric columns
  for(col in c("nsnp", "b", "se", "pval")) {
    if(col %in% colnames(df)) {
      df[[col]] <- as.numeric(df[[col]])
    } else {
      warning(paste("Column", col, "not found in dataframe"))
    }
  }
  
  return(df)
}

widen_mr_results <- function(mr_results) {
  # Create an empty list to store our results
  wide_data <- list()
  
  # Loop through each row of the original data
  for (i in 1:nrow(mr_results)) {
    method <- gsub(" ", "_", mr_results$method[i])  # Replace spaces with underscores
    
    # Add each statistic to our list
    wide_data[[paste0(method, "_nsnp")]] <- mr_results$nsnp[i]
    wide_data[[paste0(method, "_b")]] <- mr_results$b[i]
    wide_data[[paste0(method, "_se")]] <- mr_results$se[i]
    wide_data[[paste0(method, "_pval")]] <- mr_results$pval[i]
  }
  
  # Convert the list to a data frame
  wide_df <- as.data.frame(wide_data)
  
  return(wide_df)
}

extract_heterogeneity_tests <- function(file_content) {
  # Find the start of the Heterogeneity tests section
  start_line <- which(grepl("### Heterogeneity tests", file_content))
  if (length(start_line) == 0) {
    stop("Heterogeneity tests section not found")
  }
  
  # Find the start of the table (look for a line starting with '|')
  table_start <- start_line
  while (!grepl("^\\|", file_content[table_start]) && table_start < length(file_content)) {
    table_start <- table_start + 1
  }
  
  if (table_start == length(file_content)) {
    stop("Table not found in Heterogeneity tests section")
  }
  
  # Extract table lines
  table_lines <- file_content[table_start:length(file_content)]
  table_end <- which(!grepl("^\\|", table_lines))[1] - 1
  if (is.na(table_end)) {
    table_end <- length(table_lines)
  }
  table_lines <- table_lines[1:table_end]
  
  # Process each line
  data <- lapply(table_lines, function(line) {
    # Split the line by '|' and remove empty entries
    parts <- strsplit(line, "\\|")[[1]]
    parts <- parts[parts != ""]
    
    # Trim whitespace from each part
    parts <- trimws(parts)
    
    return(parts)
  })
  
  # Extract column names from the first row
  col_names <- data[[1]]
  
  # Remove the header and separator rows
  data <- data[-(1:2)]
  
  # Convert to dataframe
  df <- as.data.frame(do.call(rbind, data), stringsAsFactors = FALSE)
  if (length(col_names) > 0) {
    colnames(df) <- col_names
  }
  
  # Convert numeric columns (use actual column names)
  for(col in colnames(df)) {
    if (col != "method") {  # Assume 'method' is always non-numeric
      df[[col]] <- as.numeric(df[[col]])
    }
  }
  
  # Add "Heterogeneity_" prefix to column names
  colnames(df) <- paste0("Heterogeneity_", colnames(df))
  
  return(df)
}

wide_heterogeneity_results <- function(heterogeneity_results) {
  # Create an empty list to store our results
  wide_data <- list()
  
  # Loop through each row of the original data
  for (i in 1:nrow(heterogeneity_results)) {
    # Extract the method name, removing the "Heterogeneity_" prefix
    method <- gsub("Heterogeneity_", "", heterogeneity_results$Heterogeneity_method[i])
    method <- gsub(" ", "_", method)  # Replace spaces with underscores
    
    # Add each statistic to our list
    wide_data[[paste0("Heterogeneity_", method, "_Q")]] <- heterogeneity_results$Heterogeneity_Q[i]
    wide_data[[paste0("Heterogeneity_", method, "_Q_df")]] <- heterogeneity_results$Heterogeneity_Q_df[i]
    wide_data[[paste0("Heterogeneity_", method, "_Q_pval")]] <- heterogeneity_results$Heterogeneity_Q_pval[i]
  }
  
  # Convert the list to a data frame
  wide_df <- as.data.frame(wide_data)
  
  return(wide_df)
}

extract_pleiotropy_test <- function(file_content) {
  # Find the start of the pleiotropy test section
  start_line <- which(grepl("### Test for directional horizontal pleiotropy", file_content))
  if (length(start_line) == 0) {
    stop("Test for directional horizontal pleiotropy section not found")
  }
  
  # Find the start of the table (look for a line starting with '|')
  table_start <- start_line
  while (!grepl("^\\|", file_content[table_start]) && table_start < length(file_content)) {
    table_start <- table_start + 1
  }
  
  if (table_start == length(file_content)) {
    stop("Table not found in pleiotropy test section")
  }
  
  # Extract table lines
  table_lines <- file_content[table_start:length(file_content)]
  table_end <- which(!grepl("^\\|", table_lines))[1] - 1
  if (is.na(table_end)) {
    table_end <- length(table_lines)
  }
  table_lines <- table_lines[1:table_end]
  
  # Process each line
  data <- lapply(table_lines, function(line) {
    # Split the line by '|' and remove empty entries
    parts <- strsplit(line, "\\|")[[1]]
    parts <- parts[parts != ""]
    
    # Trim whitespace from each part
    parts <- trimws(parts)
    
    return(parts)
  })
  
  # Extract column names from the first row
  col_names <- data[[1]]
  
  # Remove the header and separator rows
  data <- data[-(1:2)]
  
  # Convert to dataframe
  df <- as.data.frame(do.call(rbind, data), stringsAsFactors = FALSE)
  if (length(col_names) > 0) {
    colnames(df) <- col_names
  }
  
  # Convert all columns to numeric
  df[] <- lapply(df, as.numeric)
  
  # Add "Pleiotropy_" prefix to column names
  colnames(df) <- paste0("Pleiotropy_", colnames(df))
  
  return(df)
}

extract_upstream_test <- function(file_content) {
  # Find the start of the upstream test section
  start_line <- which(grepl("### Test that the exposure is upstream of the outcome", file_content))
  if (length(start_line) == 0) {
    stop("Test that the exposure is upstream of the outcome section not found")
  }
  
  # Find the start of the table (look for a line starting with '|')
  table_start <- start_line
  while (!grepl("^\\|", file_content[table_start]) && table_start < length(file_content)) {
    table_start <- table_start + 1
  }
  
  if (table_start == length(file_content)) {
    stop("Table not found in upstream test section")
  }
  
  # Extract table lines
  table_lines <- file_content[table_start:length(file_content)]
  table_end <- which(!grepl("^\\|", table_lines))[1] - 1
  if (is.na(table_end)) {
    table_end <- length(table_lines)
  }
  table_lines <- table_lines[1:table_end]
  
  # Process each line
  data <- lapply(table_lines, function(line) {
    # Split the line by '|' and remove empty entries
    parts <- strsplit(line, "\\|")[[1]]
    parts <- parts[parts != ""]
    
    # Trim whitespace from each part
    parts <- trimws(parts)
    
    return(parts)
  })
  
  # Extract column names from the first row
  col_names <- data[[1]]
  
  # Remove the header and separator rows
  data <- data[-(1:2)]
  
  # Convert to dataframe
  df <- as.data.frame(do.call(rbind, data), stringsAsFactors = FALSE)
  if (length(col_names) > 0) {
    colnames(df) <- col_names
  }
  
  # Convert numeric columns to numeric
  for(col in colnames(df)) {
    if (col != "correct_causal_direction") {  # Assume this column is not numeric
      df[[col]] <- as.numeric(df[[col]])
    }
  }
  
  # Add "Upstream_" prefix to column names and replace dots with underscores
  colnames(df) <- paste0("Upstream_", gsub("\\.", "_", colnames(df)))
  
  return(df)
}

preprocess_gwas_Smith2021 <- function(gwas_path) {
  # NOTE: a2 is the effect allele in BIG40 GWAS
  # Read GWAS file
  gwas_data <- fread(gwas_path)
  
  # Add/remove requires columns
  gwas_data[, P := 10^(-`pval(-log10)`)]
  gwas_data[, `pval(-log10)` := NULL]
  gwas_data[, N := 33224] # N for smith is 33224
  
  # Merge with allele frequency
  gwas_data <- merge(gwas_data, af_Smith, by=c('rsid','chr','a1','a2','pos'), all.x=T)
  
  # Create a temporary file
  temp_file <- tempfile(fileext = ".txt")
  fwrite(gwas_data, temp_file, sep="\t")
  
  # Return the temporary file
  return(temp_file)
}

preprocess_gwas_Smith2021_QSM <- function(gwas_path) {
  # NOTE: a2 is the effect allele in BIG40 GWAS
  # Read GWAS file
  gwas_data <- fread(gwas_path)
  
  # Add/remove requires columns
  gwas_data[, P := 10^(-`pval(-log10)`)]
  gwas_data[, `pval(-log10)` := NULL]
  gwas_data[, N := 19720] # N for smith is 19720 for discovery only
  
  # Merge with allele frequency
  gwas_data <- merge(gwas_data, af_Smith, by=c('rsid','chr','a1','a2','pos'), all.x=T)
  
  # Create a temporary file
  temp_file <- tempfile(fileext = ".txt")
  fwrite(gwas_data, temp_file, sep="\t")
  
  # Return the temporary file
  return(temp_file)
}

preprocess_gwas_Warrier2023 <- function(gwas_path) {
  # Note: A1 is the effect allele in Warrier
  # Read GWAS file
  gwas_data <- fread(gwas_path)
  
  # Fix sample size
  # NOTE: In Warrier GWAS, if N = 1, then sample size is 31,797 (UKB only). If N = 2, sample size = 36,663 (UKB + ABCD)
  gwas_data <- gwas_data %>% mutate(N = case_when( N == 1 ~ 31797,
                                                   N == 2 ~ 36663,
                                                   TRUE ~ NA))
  
  # Create a temporary file
  temp_file <- tempfile(fileext = ".txt")
  fwrite(gwas_data, temp_file, sep="\t")
  
  # Return the temporary file
  return(temp_file)
}

preprocess_gwas_Zhao2020 <- function(gwas_path) {
  # Note: A1 is the effect allele in Zhao & N is present
  # Read GWAS file
  gwas_data <- fread(gwas_path)
  
  # Create a temporary file
  temp_file <- tempfile(fileext = ".txt")
  fwrite(gwas_data, temp_file, sep="\t")
  
  # Return the temporary file
  return(temp_file)
}


###################### Define GWASes and frequency ######################
## Exposure and outcome GWAS files
# PD_GWAS <- "Data/Nalls2019_allSamples_allVariants_and23andMe/PD_GWAS_2019.no23.tab"
PD_GWAS <- "Data/GP2_PD_GWAS/GP2_ALL_EUR_ALL_DATASET_HG38_12162024.rsid_withN.txt"
GWAS_files <- c(list.files('/Volumes/WorkSSD-Houman/GWAS/Smith2021_GWAS/subcortical_volume', pattern = "\\.txt$", full.names = TRUE),
                list.files('/Volumes/WorkSSD-Houman/GWAS/Smith2021_GWAS/subcortical_volume_aseg', pattern = "\\.txt$", full.names = TRUE),
                list.files('/Volumes/WorkSSD-Houman/GWAS/Smith2021_GWAS/surface_area/DKTatlas_SA', pattern = "\\.txt$", full.names = TRUE),
                list.files('/Volumes/WorkSSD-Houman/GWAS/Smith2021_GWAS/DKTatlas_CT', pattern = "\\.txt$", full.names = TRUE),
                list.files('/Volumes/WorkSSD-Houman/GWAS/Smith2021_GWAS/white_matter/WM_FA', pattern = "\\.txt$", full.names = TRUE),
                list.files('/Volumes/WorkSSD-Houman/GWAS/Smith2021_GWAS/white_matter/WM_MD', pattern = "\\.txt$", full.names = TRUE),
                '/Volumes/WorkSSD-Houman/GWAS/Warrier2023_GWAS/SA_meta.txt',
                '/Volumes/WorkSSD-Houman/GWAS/Warrier2023_GWAS/CT_meta.txt', 
                list.files('/Volumes/WorkSSD-Houman/GWAS/Smith2021_GWAS/SWI/QSM', pattern = "\\.txt$", full.names = TRUE),
                '/Volumes/WorkSSD-Houman/GWAS/Zhao2020_GWAS/FA_average_allTracts_ukb_phase1to3_dti441_dec21_2019_pheno7.fastGWA',
                '/Volumes/WorkSSD-Houman/GWAS/Zhao2020_GWAS/MD_average_allTracts_ukb_phase1to3_dti441_dec21_2019_pheno259.fastGWA')
## Allele frequency
af_Smith <- read.table('/Volumes/WorkSSD-Houman/GWAS/Smith2021_GWAS/variants_freq.txt', header=T, sep=' ') # Read the allele frequencies
af_Warrier <- read.table('/Volumes/WorkSSD-Houman/GWAS/Warrier2023_GWAS/SNP_afreq.txt', header=T, sep=' ') # Read the allele frequencies




###################### Read exposure variable (PD) ######################
if (!file.exists('Data/Exposure_MR_files_generated/PD_GWAS_GP2_2025_MR.csv')) {
  ## Read exposure data
  exp_sum=read.table(PD_GWAS, header=T, sep='\t') # GWAS summary statistics for exposure
  
  ## Prepare exposure data
  exp_filtered <- exp_sum[ exp_sum$P<5e-08, ] # filter only significant SNPs
  rm(exp_sum)

  ## Write final exposure data exposure data
  write.table(exp_filtered, 'Data/Exposure_MR_files_generated/PD_GWAS_GP2_2025_MR.csv', 
              sep=',', row.names=FALSE) # Save to a csv for further processing
  rm(exp_filtered)
}
## Read final exposure data
exp_data <- read_exposure_data('Data/Exposure_MR_files_generated/PD_GWAS_GP2_2025_MR.csv', 
                               sep=",", 
                               snp_col = "rsID", beta_col = "BETA", 
                               eaf_col = "EAF", se_col="SE", effect_allele_col = "EA", 
                               other_allele_col= "NEA", pval_col = "P", 
                               ncase_col = "N_case_proxy", ncontrol_col = "N_control",
                               clump=TRUE)






###################### Loop over each GWAS file and run MR ######################
all_results <- data.frame(matrix(ncol = 10, nrow = length(GWAS_files)))
colnames(all_results) <- c("GWAS", "MR_error", "PRESSO_outlierCorrected_p", "PRESSO_outlier_SNP_rows", 'rExp_calculated', "Fvalue_from_rExp", 'r2Exp_calculated', "Fvalue_from_r2Exp", 'rOutcome_calculated', 'r2Outcome_calculated')
all_results$GWAS <- sapply(GWAS_files, function(x) {strsplit(basename(x), "\\.txt")[[1]]})


## Loop over each file
for (file in GWAS_files) {
  current_GWAS <- strsplit(basename(file), "\\.txt")[[1]]
  # set the GWAS type
  if (current_GWAS == 'SA_meta' || current_GWAS == 'CT_meta') {GWAS_type <- 'Warrier'} else if (current_GWAS == 'FA_average_allTracts_ukb_phase1to3_dti441_dec21_2019_pheno7.fastGWA' || current_GWAS == 'MD_average_allTracts_ukb_phase1to3_dti441_dec21_2019_pheno259.fastGWA') {GWAS_type <- 'Zhao'} else if (substr(current_GWAS, 1, 3) == 'QSM') {GWAS_type <- 'Smith_QSM'} else {GWAS_type <- 'Smith'}
  print(current_GWAS)
  print(GWAS_type)
  

  # If an error comes up, ignore and move to the next iteration
  tryCatch({
    
    ############ Outcome data ############
    ## Prepare outcome variables - Smith and Warrier GWASes
    if (GWAS_type == 'Smith') {
      current_GWAS_preprocessed <- preprocess_gwas_Smith2021(file)
      out_data <- read_outcome_data(snps = exp_data$SNP,
                                    filename = current_GWAS_preprocessed, 
                                    sep="\t",  snp_col = "rsid", beta_col = "beta", se_col = "se", 
                                    eaf_col = "af", effect_allele_col = "a2", 
                                    other_allele_col = "a1", pval_col = "P",
                                    samplesize_col = "N")
      file.remove(current_GWAS_preprocessed) # remove the temporary file created by the preprocessing function
      
    } else if (GWAS_type == 'Smith_QSM') {
      current_GWAS_preprocessed <- preprocess_gwas_Smith2021_QSM(file)
      out_data <- read_outcome_data(snps = exp_data$SNP,
                                    filename = current_GWAS_preprocessed, 
                                    sep="\t",  snp_col = "rsid", beta_col = "beta", se_col = "se", 
                                    eaf_col = "af", effect_allele_col = "a2", 
                                    other_allele_col = "a1", pval_col = "P",
                                    samplesize_col = "N")
      file.remove(current_GWAS_preprocessed) # remove the temporary file created by the preprocessing function
      
    } else if (GWAS_type == 'Warrier') {
      current_GWAS_preprocessed <- preprocess_gwas_Warrier2023(file)
      out_data <- read_outcome_data(snps = exp_data$SNP,
                                    filename = current_GWAS_preprocessed, 
                                    sep="\t",  snp_col = "SNP", beta_col = "BETA", se_col = "SE", 
                                    eaf_col = "AF1", effect_allele_col = "A1", 
                                    other_allele_col = "A2", pval_col = "P",
                                    samplesize_col = "N")
      file.remove(current_GWAS_preprocessed) # remove the temporary file created by the preprocessing function
    } else if (GWAS_type == 'Zhao') {
      current_GWAS_preprocessed <- preprocess_gwas_Zhao2020(file)
      out_data <- read_outcome_data(snps = exp_data$SNP,
                                    filename = current_GWAS_preprocessed, 
                                    sep="\t",  snp_col = "SNP", beta_col = "BETA", se_col = "SE", 
                                    eaf_col = "AF1", effect_allele_col = "A1", 
                                    other_allele_col = "A2", pval_col = "P",
                                    samplesize_col = "N")
      file.remove(current_GWAS_preprocessed) # remove the temporary file created by the preprocessing function
    }
    
    # add R values and harmonize
    out_data$r.outcome <- get_r_from_pn(out_data$pval.outcome, out_data$samplesize.outcome)
    dat <- harmonise_data(exposure_dat=exp_data, outcome_dat=out_data, action=2)
    
    
    
    
    ############ Running MR ############
    dat$units.outcome <-"mm2"
    dat$units.exposure <-"prob"
    dat$r.exposure <- get_r_from_lor(dat$beta.exposure, dat$eaf.exposure, dat$ncase.exposure, dat$ncontrol.exposure,
                                     prevalence=0.01, model = "logit")
    
    steiger <- steiger_filtering(dat)
    sig <- subset(steiger, steiger$steiger_dir==TRUE) # Only take SNPs remaining after Steiger filtering
    presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                        SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                        OUTLIERtest = F, DISTORTIONtest = T, data = sig, 
                        NbDistribution = 1000,  SignifThreshold = 0.05)
    print(presso)
    # save presso to a text file
    dir.create(paste0('Outputs/reverse_MR/PD_vs_', current_GWAS), recursive = TRUE)
    sink(paste0('Outputs/reverse_MR/PD_vs_', current_GWAS, '/', current_GWAS, '_pressoResults.txt'))
    print(presso)
    sink()
    
    all_results$PRESSO_outlierCorrected_p[which(all_results$GWAS == current_GWAS)] <- presso$`Main MR results`$`P-value`[2]
    if (paste(which(presso$`MR-PRESSO results`$`Outlier Test`[2] < 0.05), collapse = "/") != '') {
      all_results$PRESSO_outlier_SNP_rows[which(all_results$GWAS == current_GWAS)] <- paste(which(presso$`MR-PRESSO results`$`Outlier Test`[2] < 0.05), collapse = "/")
    }
    mr(sig)
    mr_report(sig, output_path=paste0('Outputs/reverse_MR/PD_vs_', current_GWAS))
    singlesnp_plot=mr_forest_plot(mr_singlesnp(sig))
    ggsave(singlesnp_plot[[1]], file= "plot.jpg", 
           path=paste0('Outputs/reverse_MR/PD_vs_', current_GWAS), 
           width=7, height=7)
    
    
    ### get and save all the MR results
    report_content <- readLines(paste0('Outputs/reverse_MR/PD_vs_', current_GWAS, '/mr_report.md'))
    mr_results <- extract_mr_results_table(report_content)
    write.csv(mr_results, paste0('Outputs/reverse_MR/PD_vs_', current_GWAS, '/', current_GWAS, '_mrResults.csv'))
    wide_mr_results <- widen_mr_results(mr_results)
    # get and save heterogeneity test
    heterogeneity_results <- extract_heterogeneity_tests(report_content)
    write.csv(heterogeneity_results, paste0('Outputs/reverse_MR/PD_vs_', current_GWAS, '/', current_GWAS, '_heterogeneityResults.csv'))
    wide_het_results <- wide_heterogeneity_results(heterogeneity_results)
    # get and save pleiotropy test
    pleiotropy_results <- extract_pleiotropy_test(report_content)
    write.csv(pleiotropy_results, paste0('Outputs/reverse_MR/PD_vs_', current_GWAS, '/', current_GWAS, '_pleiotropyResults.csv'))
    # get exposure upstream of outcome test
    upstream_results <- extract_upstream_test(report_content)
    write.csv(upstream_results, paste0('Outputs/reverse_MR/PD_vs_', current_GWAS, '/', current_GWAS, '_exposureUpstreamResults.csv'))
    
    # add them all to the all_results file
    current_gwas_row <- which(all_results$GWAS == current_GWAS)
    for (current_col in colnames(wide_mr_results)) { all_results[current_gwas_row, current_col] <- wide_mr_results[current_col] }
    for (current_col in colnames(wide_het_results)) { all_results[current_gwas_row, current_col] <- wide_het_results[current_col] }
    for (current_col in colnames(pleiotropy_results)) { all_results[current_gwas_row, current_col] <- pleiotropy_results[current_col] }
    for (current_col in colnames(upstream_results)) { all_results[current_gwas_row, current_col] <- upstream_results[current_col] }
    
    # Calculate F-value both ways and add it
    R2_rExp <- mean(sig$r.exposure)
    R2_r2Exp <- sum(sig$rsq.exposure)
    k <- nrow(sig)
    n <- mean(sig$samplesize.exposure)
    F_value_r <- (R2_rExp*(n-1-k))/((1-R2_rExp)*k)
    F_value_r2 <- (R2_r2Exp*(n-1-k))/((1-R2_r2Exp)*k)
    all_results$Fvalue_from_rExp[which(all_results$GWAS == current_GWAS)] <- F_value_r
    all_results$Fvalue_from_r2Exp[which(all_results$GWAS == current_GWAS)] <- F_value_r2
    # Save other values
    all_results$rExp_calculated[which(all_results$GWAS == current_GWAS)] <- mean(sig$r.exposure)
    all_results$r2Exp_calculated[which(all_results$GWAS == current_GWAS)] <- sum(sig$rsq.exposure)
    all_results$rOutcome_calculated[which(all_results$GWAS == current_GWAS)] <- mean(sig$r.outcome)
    all_results$r2Outcome_calculated[which(all_results$GWAS == current_GWAS)] <- sum(sig$rsq.outcome)
    
    
    
    
    
    
    
    
    
    
    
  }, error = function(e) {
    # Skip this iteration and move to the next
    all_results$MR_error[which(all_results$GWAS == current_GWAS)] <<- e$message
    message(paste("Skipping", current_GWAS, "due to error:", e$message))
    return()
  })
}

write.csv(all_results, paste0('Outputs/result_summary_allMeasures/reverseMR_PD_vsBrain_MR_Results_complete_QSM_Zhao.csv'))









