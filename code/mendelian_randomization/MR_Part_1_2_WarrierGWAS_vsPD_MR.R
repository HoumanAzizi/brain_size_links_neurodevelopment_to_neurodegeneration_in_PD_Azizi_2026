rm(list = ls())
cat("\014")

#install.packages('remotes')
Sys.setenv(GITHUB_PAT = "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") # Personal github token
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


###################### Define GWASes and frequency ######################
# Location of GWAS files: ~/projects/rrg-adagher/public_data/CM_gwas_0920
af=read.table('/Volumes/WorkSSD-Houman/GWAS/Warrier2023_GWAS/SNP_afreq.txt', header=T, sep=' ') # Read the allele frequencies
GWAS_files <- c('/Volumes/WorkSSD-Houman/GWAS/Warrier2023_GWAS/SA_meta.txt','/Volumes/WorkSSD-Houman/GWAS/Warrier2023_GWAS/CT_meta.txt')








###################### Loop over each GWAS file and run MR ######################
## Define the table to keep all results
all_results <- data.frame(matrix(ncol = 10, nrow = length(GWAS_files)))
colnames(all_results) <- c("GWAS", "MR_error", "PRESSO_outlierCorrected_p", "PRESSO_outlier_SNP_rows", 'rExp_calculated', "Fvalue_from_rExp", 'r2Exp_calculated', "Fvalue_from_r2Exp", 'rOutcome_calculated', 'r2Outcome_calculated')
all_results$GWAS <- sapply(GWAS_files, function(x) {strsplit(basename(x), "\\.txt")[[1]]})

## Loop over each file
for (file in GWAS_files) {
  current_GWAS <- strsplit(basename(file), "\\.txt")[[1]]
  print(current_GWAS)
  
  
  # If an error comes up, ignore and move to the next iteration
  tryCatch({
    ############ Exposure Data ############
    ## If have previously read exposure data, do not re-read it again
    if (!file.exists(paste0('Data/Exposure_MR_files_generated/', current_GWAS, '_Warrier_MR.csv'))) {
      ## Read exposure data
      exp_sum=read.table(file, header=T, sep=' ') # GWAS summary statistics for exposure
      
      ## Prepare exposure data
      exp_filtered = exp_sum[ exp_sum$P<5e-08, ] # filter only significant SNPs
      rm(exp_sum)
      # Fix sample size
      # NOTE: In Warrier GWAS, if N = 1, then sample size is 31,797 (UKB only). If N = 2, sample size = 36,663 (UKB + ABCD)
      exp_filtered <- exp_filtered %>% mutate(N = case_when( N == 1 ~ 31797,
                                                             N == 2 ~ 36663,
                                                             TRUE ~ NA))
      
      exp_filtered=merge(exp_filtered, af, by=c('SNP','CHR','A1','A2'), all.x=T) # Merge with allele frequencies
      
      ## Write final exposure data exposure data
      write.table(exp_filtered, paste0('Data/Exposure_MR_files_generated/', current_GWAS, '_Warrier_MR.csv'), 
                  sep=',', row.names=FALSE) # Save to a csv for further processing
      rm(exp_filtered)
    }
    ## Read final exposure data
    exp_data <- read_exposure_data(paste0('Data/Exposure_MR_files_generated/', current_GWAS, '_Warrier_MR.csv'), 
                                   sep=",", 
                                   snp_col = "SNP", beta_col = "BETA", 
                                   eaf_col = "AF1", se_col="SE", effect_allele_col = "A1", 
                                   other_allele_col= "A2", pval_col = "P", 
                                   samplesize_col = "N",  clump=TRUE)
    
    
    
    
    
    ############ Outcome Data ############
    ## Prepare outcome variables - PD here
    out_data <- read_outcome_data(snps = exp_data$SNP,
                                  filename = "Data/GP2_PD_GWAS/GP2_ALL_EUR_ALL_DATASET_HG38_12162024.rsid_withN.txt", # Using the NEW GP2 PD GWAS
                                  sep="\t",  snp_col = "rsID", beta_col = "BETA", se_col = "SE", 
                                  eaf_col = "EAF", effect_allele_col = "EA", 
                                  other_allele_col = "NEA", pval_col = "P",
                                  ncase_col = "N_case_proxy", ncontrol_col = "N_control") 
    # for previous PD GWAS Nalls
    # out_data <- read_outcome_data(snps = exp_data$SNP,
    #                               filename = "Data/Nalls2019_allSamples_allVariants_and23andMe/PD_GWAS_2019.no23.tab", 
    #                               sep="\t",  snp_col = "SNP", beta_col = "b", se_col = "se", 
    #                               eaf_col = "FREQ", effect_allele_col = "A1", 
    #                               other_allele_col = "A2", pval_col = "p",
    #                               ncase_col = "N_cases", ncontrol_col = "N_controls") 
    
    out_data$r.outcome <- get_r_from_lor(out_data$beta.outcome, out_data$eaf.outcome, out_data$ncase.outcome, out_data$ncontrol.outcome, 
                                         prevalence=0.01, model = "logit")
    dat <- harmonise_data(exposure_dat=exp_data, outcome_dat=out_data, action=2)
    
    ## EXCLUDE PLEIOTROPIC SNP rs56323304 - done post-hoc ###
    # dat=dat[!dat$SNP=='rs56323304',] # for SA_meta
    # dat=dat[!dat$SNP=='rs73802707',] # for SA_meta
    
    
    
    
    ############ Running MR ############
    dat$units.outcome <-"prob"
    dat$units.exposure <-"mm2"
    dat$r.exposure <- get_r_from_pn(dat$pval.exposure, dat$samplesize.exposure)

    steiger <- steiger_filtering(dat)
    sig <- subset(steiger, steiger$steiger_dir==TRUE) # Only take SNPs remaining after Steiger filtering
    presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                        SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                        OUTLIERtest = F, DISTORTIONtest = T, data = sig, 
                        NbDistribution = 1000,  SignifThreshold = 0.05)
    print(presso)
    # save presso to a text file
    dir.create(paste0('Outputs/', current_GWAS, '_vs_PD'), recursive = TRUE)
    sink(paste0('Outputs/', current_GWAS, '_vs_PD/', current_GWAS, '_pressoResults.txt'))
    print(presso)
    sink()
    
    all_results$PRESSO_outlierCorrected_p[which(all_results$GWAS == current_GWAS)] <- presso$`Main MR results`$`P-value`[2]
    if (paste(which(presso$`MR-PRESSO results`$`Outlier Test`[2] < 0.05), collapse = "/") != '') {
      all_results$PRESSO_outlier_SNP_rows[which(all_results$GWAS == current_GWAS)] <- paste(which(presso$`MR-PRESSO results`$`Outlier Test`[2] < 0.05), collapse = "/")
    }
    mr(sig)
    mr_report(sig, output_path=paste0('Outputs/', current_GWAS, '_vs_PD'))
    singlesnp_plot=mr_forest_plot(mr_singlesnp(sig))
    ggsave(singlesnp_plot[[1]], file= "plot.jpg", 
           path=paste0('Outputs/', current_GWAS, '_vs_PD'), 
           width=7, height=7)
    
    ### Save better figures separately
    ### Generate and save publication-ready figures
    results <- mr_report(dat)
    m <- list(
      mr = mr(dat), 
      mr_singlesnp = mr_singlesnp(dat), 
      mr_leaveoneout = mr_leaveoneout(dat)
    )
    
    # Create output directory if it doesn't exist
    output_dir <- paste0('Outputs/', current_GWAS, '_vs_PD', '/figure/')
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Generate base plots
    p <- list(
      scatter = mr_scatter_plot(m$mr, dat),
      forest = mr_forest_plot(m$mr_singlesnp),
      funnel = mr_funnel_plot(m$mr_singlesnp),
      leaveoneout = mr_leaveoneout_plot(m$mr_leaveoneout)
    )
    
    # Customize and save scatter plot
    # Get the base scatter plot
    scatter_plot <- p$scatter[[1]]
    
    # Store the regression line layer (layer 4)
    regression_layer <- scatter_plot$layers[[4]]
    
    # Remove the regression line layer temporarily
    scatter_plot$layers[[4]] <- NULL
    
    # Customize the existing error bars and points
    scatter_plot$layers[[1]]$aes_params$colour <- "grey90"
    scatter_plot$layers[[2]]$aes_params$colour <- "grey90"
    scatter_plot$layers[[3]]$aes_params$colour <- "grey40"
    scatter_plot$layers[[3]]$aes_params$size <- 2
    
    # Make regression lines thicker
    regression_layer$aes_params$size <- 0.8
    
    # Re-add the regression line layer on top
    scatter_plot$layers[[4]] <- regression_layer
    
    # Apply your custom theme
    scatter_plot <- scatter_plot +
      theme_minimal() +
      theme(
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.direction = "vertical",
        legend.position = "right",
        legend.background = element_rect(color = "black", fill = "white"),
        legend.box.background = element_rect(color = "black", fill = "white"),
        axis.line = element_line(color = "black", size = 0.5)
      ) +
      guides(colour = guide_legend(ncol = 1))
    
    ggsave(paste0(output_dir, "scatter_plot.svg"), scatter_plot, 
           width = 8, height = 6, units = "in", device = "svg")
    
    # Customize and save forest plot
    forest_plot <- p$forest[[1]] +
      theme_minimal() +
      theme(
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "white", size = 0.5),
        legend.position = "none"
      )
    
    ggsave(paste0(output_dir, "forest_plot.svg"), forest_plot, 
           width = 6, height = 6.5, units = "in", device = "svg")
    
    # Save funnel plot
    funnel_plot <- p$funnel[[1]] +
      theme_minimal() +
      theme(
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 0.5)
      )
    
    ggsave(paste0(output_dir, "funnel_plot.svg"), funnel_plot, 
           width = 8, height = 6, units = "in", device = "svg")
    
    # Save leave-one-out plot
    leaveoneout_plot <- p$leaveoneout[[1]] +
      theme_minimal() +
      theme(
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 0.5)
      )
    
    ggsave(paste0(output_dir, "leaveoneout_plot.svg"), leaveoneout_plot, 
           width = 6, height = 6.5, units = "in", device = "svg")
    ### end of new plots ###
    
    
    ### get and save all the MR results
    report_content <- readLines(paste0('Outputs/', current_GWAS, '_vs_PD/mr_report.md'))
    mr_results <- extract_mr_results_table(report_content)
    write.csv(mr_results, paste0('Outputs/', current_GWAS, '_vs_PD/', current_GWAS, '_mrResults.csv'))
    wide_mr_results <- widen_mr_results(mr_results)
    # get and save heterogeneity test
    heterogeneity_results <- extract_heterogeneity_tests(report_content)
    write.csv(heterogeneity_results, paste0('Outputs/', current_GWAS, '_vs_PD/', current_GWAS, '_heterogeneityResults.csv'))
    wide_het_results <- wide_heterogeneity_results(heterogeneity_results)
    # get and save pleiotropy test
    pleiotropy_results <- extract_pleiotropy_test(report_content)
    write.csv(pleiotropy_results, paste0('Outputs/', current_GWAS, '_vs_PD/', current_GWAS, '_pleiotropyResults.csv'))
    # get exposure upstream of outcome test
    upstream_results <- extract_upstream_test(report_content)
    write.csv(upstream_results, paste0('Outputs/', current_GWAS, '_vs_PD/', current_GWAS, '_exposureUpstreamResults.csv'))
    
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

write.csv(all_results, paste0('Outputs/result_summary_allMeasures/Warrier_GWAS_MR_Results_complete.csv'))



