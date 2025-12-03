library(tidyr)
library(dplyr)
library(lme4)
library(readr)
#library(ggplot2)

rm(list = ls())
cat("\014")

setwd("/Users/houmanazizi/Library/Mobile\ Documents/com~apple~CloudDocs/Education/GitHub/UKB_Analysis/")




############################# Setting up Confounds + Go Through Each Version #############################
n_PCA <- 15
PCAs <- "PCA_1"
for (i in 2:n_PCA) {
  PCAs <- paste0(PCAs," + PCA_",i)
}


#### SET POSSIBLE OPTIONS HERE FOR THE MODEL
model_sex_options <- c('all', 'male', 'female')


# Create main output path
main_output_folder <- 'Outputs/PD_PRS_Regression_Results/with_all_pathwayPRS_Behavioral_fullUKB/'
dir.create(main_output_folder, recursive = FALSE, showWarnings = FALSE)


#### Loop through all combinations
for (model_sex in model_sex_options) {
  #### Automatic model selection (updated model)
  if (model_sex == 'all') {
    Confounds <- paste0('Age + Age2 + Sex + Age:Sex + Age2:Sex + Center_Number + Genotype_batch + ', PCAs)
    sex_to_include <- c(1,0)
  } else if (model_sex == 'male') {
    sex_to_include <- 1 # Male (will be 1)
    Confounds <- paste0('Age + Age2 + Center_Number + Genotype_batch + ', PCAs)
  } else if (model_sex == 'female') {
    sex_to_include <- 0 # Female (will be 0)
    Confounds <- paste0('Age + Age2 + Center_Number + Genotype_batch + ', PCAs)
  } else { print('WRONG SEX SELECTED') }
  
  print(Confounds)
  
  
  
  
  ############################# Read the data or recreate them if not ready #############################
  if (!file.exists('Data/Current_UKB/UKB_Behavioral/UKB_BehavioralMeasures_vs_PDPRS_Data_Nov2025_fullUKB.csv')) { # Sept2025 version does not have sleep
    ############################# Reading Data #############################
    # UKB Data
    # UKB_BL <- as.data.frame(read_delim('Data/Current_UKB/UKB_Behavioral/subset_2025_May21_BehFigure_wide.tsv', delim = "\t", col_names = TRUE, na = "", quote = ""))
    UKB_BL <- as.data.frame(read_delim('Data/Current_UKB/UKB_Behavioral/subset_2025_Sept4_Beh_newMeasures_wide.tsv', delim = "\t", col_names = TRUE, na = "", quote = ""))
    
    # Read and add PCA
    UKB_PCA <- read.csv("Data/Genetic/PCA/UKB_PCAs_Lang.csv")
    UKB_BL <- UKB_BL %>% left_join(UKB_PCA, by = 'SubjectID')
    rm(UKB_PCA)
    
    # Read Disorders and exclude
    UKB_only_Disorders <- read.csv("Data/Exclusions/UKB_Disorders_Counts_fullUKB.csv")
    ids_to_remove <- UKB_only_Disorders %>% filter(All_Neuro_Disorder_ID_Grouped != 0)
    ids_to_remove <- ids_to_remove %>% select(SubjectID)
    ids_to_remove <- unique(ids_to_remove$SubjectID)
    UKB_BL <- UKB_BL %>% filter(!(SubjectID %in% ids_to_remove))
    rm(UKB_only_Disorders,ids_to_remove)
    
    
    # Separate UKB_BL from UKB_beh
    UKB_beh <- UKB_BL %>% select("SubjectID", "InstanceID", "ArrayID", "Sex_31", "Average total household income before tax_738", 
                                 "Alcohol intake frequency._1558", "Qualifications_6138", "Fluid intelligence score_20016", "Coffee intake_1498",
                                 "Ever addicted to any substance or behaviour_20401", "Amount of alcohol drunk on a typical drinking day_20403", "Frequency of drinking alcohol_20414", 
                                 "Body mass index (BMI)_21001", "Age when attended assessment centre_21003", "Genetic sex_22001", "Index of Multiple Deprivation (England)_26410", 
                                 "Pack years of smoking_20161", "Ever addicted to a behaviour or miscellanous_20431", 
                                 "Ever addicted to illicit or recreational drugs_20456", "Ever addicted to prescription or over-the-counter medication_20503", "Behavioural and miscellaneous addictions_20552", "Substance of prescription or over-the-counter medication addiction_20551",
                                 "Hand grip strength (left)_46", "Hand grip strength (right)_47", 
                                 "Frequency of depressed mood in last 2 weeks_2050", "Frequency of unenthusiasm / disinterest in last 2 weeks_2060", 
                                 "Diagnoses - main ICD10_41202",
                                 "Sleep duration_1160")
  
    UKB_BL <- UKB_BL %>% filter(InstanceID == 0) %>% filter(ArrayID == 0) %>% select(c("SubjectID", "Sex_31", "Age when attended assessment centre_21003", "UK Biobank assessment centre_54", "PCA_1", "PCA_2", "PCA_3", "PCA_4", "PCA_5",
                                                                                       "PCA_6", "PCA_7", "PCA_8", "PCA_9", "PCA_10", "PCA_11", "PCA_12", "PCA_13", "PCA_14", "PCA_15"), starts_with('PD_PRS'), "Genotype measurement batch_22000")
    # Read and add new PD-PRS
    PRS <- read.csv('Data/Genetic/PD_PRS/PD_PRS_allVersions_Final_mainPathways.csv') # has all PRS values
    UKB_BL <- UKB_BL %>% left_join(PRS, by = "SubjectID") %>% arrange(SubjectID)
    UKB_BL <- UKB_BL %>% filter(!is.na(PD_PRScs_Nalls))
    rm(PRS)
    
    
    
    
    
    ############################# Fix multi-instance UKB_beh data #############################
    ## Fixing Constipation, Urinary incontinence, and Orthostatic Hypotension
    UKB_beh <- UKB_beh %>% mutate(Constipation = case_when(`Diagnoses - main ICD10_41202` == 'K59.0 Constipation' ~ 1, TRUE ~ 0))
    UKB_beh <- UKB_beh %>% mutate(Orthostatic_hypotension = case_when(`Diagnoses - main ICD10_41202` == 'I95.1 Orthostatic hypotension' ~ 1, TRUE ~ 0))
    UKB_beh <- UKB_beh %>% mutate(Urinary_incontinence = case_when(`Diagnoses - main ICD10_41202` %in% c('N39.3 Stress incontinence', 'N39.4 Other specified urinary incontinence', 'R32 Unspecified urinary incontinence') ~ 1, TRUE ~ 0))
    UKB_beh <- UKB_beh %>% group_by(SubjectID) %>% mutate( Constipation = max(Constipation, na.rm = TRUE),
                                                           Urinary_incontinence = max(Urinary_incontinence, na.rm = TRUE),
                                                           Orthostatic_hypotension = max(Orthostatic_hypotension, na.rm = TRUE)) %>% ungroup()
    UKB_beh <- UKB_beh %>% select(-`Diagnoses - main ICD10_41202`)
    
    
    
    ############################# Fix the UKB_beh data #############################
    UKB_beh <- UKB_beh %>% filter(InstanceID == 0)
    ## Fixing income - using midpoint instead of 1-5
    UKB_beh <- UKB_beh %>% mutate(Household_income = case_when(`Average total household income before tax_738` == "Less than 18,000" ~ 9000,
                                                               `Average total household income before tax_738` == "18,000 to 30,999" ~ 24500,
                                                               `Average total household income before tax_738` == "31,000 to 51,999" ~ 41500,
                                                               `Average total household income before tax_738` == "52,000 to 100,000" ~ 76000,
                                                               `Average total household income before tax_738` == "Greater than 100,000" ~ 125000,
                                                               TRUE ~ NA)) %>% select(-`Average total household income before tax_738`)
    
    ## Fixing alcohol
    UKB_beh <- UKB_beh %>% mutate(Alcohol_drinks_per_month = case_when(`Frequency of drinking alcohol_20414` == "Never" ~ 0,
                                                                       `Frequency of drinking alcohol_20414` == "Monthly or less" ~ 1,
                                                                       `Frequency of drinking alcohol_20414` == "2 to 4 times a month" ~ 3,
                                                                       `Frequency of drinking alcohol_20414` == "2 to 3 times a week" ~ 10,
                                                                       `Frequency of drinking alcohol_20414` == "4 or more times a week" ~ 20,
                                                                       TRUE ~ NA)) %>% select(-`Frequency of drinking alcohol_20414`)
    UKB_beh <- UKB_beh %>% mutate(Alcohol_amount_per_drink = case_when(`Amount of alcohol drunk on a typical drinking day_20403` == "1 or 2" ~ 1,
                                                                       `Amount of alcohol drunk on a typical drinking day_20403` == "3 or 4" ~ 3,
                                                                       `Amount of alcohol drunk on a typical drinking day_20403` == "5 or 6" ~ 5,
                                                                       `Amount of alcohol drunk on a typical drinking day_20403` == "7, 8 or 9" ~ 7,
                                                                       `Amount of alcohol drunk on a typical drinking day_20403` == "10 or more" ~ 10,
                                                                       TRUE ~ NA)) %>% select(-`Amount of alcohol drunk on a typical drinking day_20403`)
    UKB_beh <- UKB_beh %>% mutate(Alcohol_usage_total = Alcohol_amount_per_drink*Alcohol_drinks_per_month) %>% 
      select(-Alcohol_drinks_per_month, -Alcohol_amount_per_drink, -`Alcohol intake frequency._1558`)
    
    ## Fixing Education based on ISCED number of years (check https://doi.org/10.1038/s41588-018-0147-3)
    UKB_beh <- UKB_beh %>% mutate(Education_level = case_when(Qualifications_6138 == "None of the above" ~ 7,
                                                              Qualifications_6138 == "NVQ or HND or HNC or equivalent" ~ 19,
                                                              Qualifications_6138 == "CSEs or equivalent" ~ 10,
                                                              Qualifications_6138 == "O levels/GCSEs or equivalent" ~ 10,
                                                              Qualifications_6138 == "A levels/AS levels or equivalent" ~ 13,
                                                              Qualifications_6138 == "College or University degree" ~ 20,
                                                              Qualifications_6138 == "Other professional qualifications eg: nursing, teaching" ~ 15,
                                                              TRUE ~ NA))
    UKB_beh <- UKB_beh %>% group_by(SubjectID) %>% mutate(Education_level = ifelse(all(is.na(Education_level)), NA, max(Education_level, na.rm = TRUE))) %>%
      select(-Qualifications_6138)
    UKB_beh <- as.data.frame(UKB_beh)
    
    ## Fixing Sex
    UKB_beh <- UKB_beh %>% mutate(Sex = case_when(Sex_31 == "Male" ~ 1,
                                                  Sex_31 == "Female" ~ 2,
                                                  TRUE ~ NA)) %>% select(-`Genetic sex_22001`, -Sex_31)
    
    ## Fixing addiction
    # For columns with possible multiple answers, do SUM
    UKB_beh <- UKB_beh %>% mutate(field_20551 = case_when(!is.na(`Substance of prescription or over-the-counter medication addiction_20551`) ~ 1,
                                                          TRUE ~ 0)) %>% select(-`Substance of prescription or over-the-counter medication addiction_20551`)
    UKB_beh <- UKB_beh %>% mutate(field_20552 = case_when(!is.na(`Behavioural and miscellaneous addictions_20552`) ~ 1,
                                                          TRUE ~ 0)) %>% select(-`Behavioural and miscellaneous addictions_20552`)
    UKB_beh <- UKB_beh %>% group_by(SubjectID) %>% mutate(field_20551 = sum(field_20551), field_20552=sum(field_20552))
    UKB_beh <- as.data.frame(UKB_beh)
    # Remove all rows with NA in 20401 - here NAs mean no response
    UKB_beh <- UKB_beh %>% filter(!is.na(`Ever addicted to any substance or behaviour_20401`))
    # Fix column 20431, accounting for 20552
    UKB_beh <- UKB_beh %>% mutate(addiction_beh = case_when(`Ever addicted to a behaviour or miscellanous_20431`=='Yes' & field_20552==0 ~ 1,
                                                            `Ever addicted to a behaviour or miscellanous_20431`=='Yes' & field_20552>0 ~ field_20552,
                                                            TRUE ~ 0)) %>% select(-field_20552, -`Ever addicted to a behaviour or miscellanous_20431`)
    # Fix column 20503, accounting for 20551
    UKB_beh <- UKB_beh %>% mutate(addiction_med = case_when(`Ever addicted to prescription or over-the-counter medication_20503`=='Yes' & field_20551==0 ~ 1,
                                                            `Ever addicted to prescription or over-the-counter medication_20503`=='Yes' & field_20551>0 ~ field_20551,
                                                            TRUE ~ 0)) %>% select(-field_20551, -`Ever addicted to prescription or over-the-counter medication_20503`)
    # Fix column 20456
    UKB_beh <- UKB_beh %>% mutate(addiction_drug = case_when(`Ever addicted to illicit or recreational drugs_20456`=='Yes' ~ 1,
                                                             TRUE ~ 0)) %>% select(-`Ever addicted to illicit or recreational drugs_20456`)
    # Calculate a final addiction score
    UKB_beh <- UKB_beh %>% mutate(addiction_sum = addiction_beh+addiction_med+addiction_drug) %>% select(-addiction_beh,-addiction_med,-addiction_drug)
    UKB_beh <- UKB_beh %>% mutate(Addiction_score_total = case_when(`Ever addicted to any substance or behaviour_20401` == 'Yes' & addiction_sum==0 ~ 1,
                                                                    `Ever addicted to any substance or behaviour_20401` == 'Yes' & addiction_sum>0 ~ addiction_sum,
                                                                    TRUE ~ 0)) %>% select(-`Ever addicted to any substance or behaviour_20401`, -addiction_sum)
    
    ## Fix coffee consumption
    UKB_beh$`Coffee intake_1498`[UKB_beh$`Coffee intake_1498` == "Less than one"] <- 0.5
    UKB_beh$`Coffee intake_1498` <- as.integer(UKB_beh$`Coffee intake_1498`)
    
    ## Fix hand grip strength
    UKB_beh <- UKB_beh %>% mutate(Hand_grip_strength = (`Hand grip strength (left)_46` + `Hand grip strength (right)_47`)/2) %>% select(-`Hand grip strength (left)_46`, -`Hand grip strength (right)_47`)
    
    ## Fix frequency of depression and apathy in past 2 weeks
    UKB_beh <- UKB_beh %>% mutate(Apathy_past2weeks = case_when(`Frequency of unenthusiasm / disinterest in last 2 weeks_2060` == 'Not at all' ~ 0,
                                                                `Frequency of unenthusiasm / disinterest in last 2 weeks_2060` == 'Several days' ~ 5,
                                                                `Frequency of unenthusiasm / disinterest in last 2 weeks_2060` == 'More than half the days' ~ 10,
                                                                `Frequency of unenthusiasm / disinterest in last 2 weeks_2060` == 'Nearly every day' ~ 14,
                                                                TRUE ~ NA)) %>% select(-`Frequency of unenthusiasm / disinterest in last 2 weeks_2060`)
    UKB_beh <- UKB_beh %>% mutate(Depressed_mood_past2weeks = case_when(`Frequency of depressed mood in last 2 weeks_2050` == 'Not at all' ~ 0,
                                                                        `Frequency of depressed mood in last 2 weeks_2050` == 'Several days' ~ 5,
                                                                        `Frequency of depressed mood in last 2 weeks_2050` == 'More than half the days' ~ 10,
                                                                        `Frequency of depressed mood in last 2 weeks_2050` == 'Nearly every day' ~ 14,
                                                                        TRUE ~ NA)) %>% select(-`Frequency of depressed mood in last 2 weeks_2050`)
    
    
    
    ## Set column names
    UKB_beh <- UKB_beh %>% select(-starts_with('Volume'), -ArrayID, -InstanceID)
    colnames(UKB_beh)
    new_col_names <- c("SubjectID", "Fluid_intelligence_score", "Coffee_intake", "BMI", "Age", 
                       "Multiple_deprivation_index", "Smoking_PackYears", "Sleep_duration",
                       "Constipation", "Orthostatic_hypotension", "Urinary_incontinence",
                       "Household_income", "Alcohol_usage_total", "Education_level", "Sex", "Addiction_score_total",
                       "Hand_grip_strength", "Apathy_past2weeks", "Depressed_mood_past2weeks" )
    colnames(UKB_beh) <- new_col_names
    UKB_beh <- UKB_beh %>% filter(SubjectID %in% UKB_BL$SubjectID)
    
    
    
    
    
    
    
    
    
    
    
    
    ############################# Keeping Confounds for all #############################
    UKB_data_final <- UKB_BL %>% select(c("SubjectID", "Sex_31", "Age when attended assessment centre_21003", "UK Biobank assessment centre_54", "PCA_1", "PCA_2", "PCA_3", "PCA_4", "PCA_5",
                                          "PCA_6", "PCA_7", "PCA_8", "PCA_9", "PCA_10", "PCA_11", "PCA_12", "PCA_13", "PCA_14", "PCA_15"), starts_with('PD_PRS'))
    
    colnames(UKB_data_final)[3:4] <- c("Age", "cn")
    
    UKB_data_final <- UKB_data_final %>% mutate(Sex = case_when(Sex_31 == "Male" ~ 1,
                                                                Sex_31 == "Female" ~ 0)) %>% select(-c("Sex_31"))
    
    UKB_data_final$Age2 <- (UKB_data_final$Age)^2
    
    
    # UKB_data_final$gn <- substr(UKB_BL$`Genotype measurement batch_22000`,start=1,stop=3)
    # UKB_data_final <- UKB_data_final %>%
    #   mutate(Genotype_batch = case_when(gn == "Bat" ~ 2,
    #                                     gn == "UKB" ~ 1)) %>%
    #   select(-c(gn)) #UKBiLEVEAX=1
    
    # Create a mapping of each unique batch to a number
    batches <- unique(UKB_BL$`Genotype measurement batch_22000`)
    batch_map <- setNames(1:length(batches), batches)
    # Apply the mapping to create the Genotype_batch variable
    UKB_data_final <- UKB_data_final %>% mutate(Genotype_batch = batch_map[UKB_BL$`Genotype measurement batch_22000`])
    
    
    UKB_data_final <- UKB_data_final %>% 
      mutate(Center_Number = case_when(
        cn == "Bury" ~ 1,
        cn == "Birmingham" ~ 2,
        cn == "Oxford" ~ 3,
        cn == "Liverpool" ~ 4,
        cn == "Barts" ~ 5,
        cn == "Sheffield" ~ 6,
        cn == "Leeds" ~ 7,
        cn == "Stoke" ~ 8,
        cn == "Manchester" ~ 9,
        cn == "Reading" ~ 10,
        cn == "Nottingham" ~ 11,
        cn == "Croydon" ~ 12,
        cn == "Newcastle" ~ 13,
        cn == "Bristol" ~ 14,
        cn == "Hounslow" ~ 15,
        cn == "Edinburgh" ~ 16,
        cn == "Middlesborough" ~ 17,
        cn == "Glasgow" ~ 18,
        cn == "Cardiff" ~ 19,
        cn == "Swansea" ~ 20,
        cn == "Wrexham" ~ 21,
        cn == "Stockport (pilot)" ~ 22 )) %>% select(-c(cn))
    
    
    UKB_BL <- UKB_data_final
    rm(UKB_data_final)
    
    
    
    
    
    
    ############################# Add Behavioral Data #############################
    colnames(UKB_BL)
    colnames(UKB_beh)
    UKB_beh <- UKB_beh %>% select(-Age, -Sex)
    UKB_BL <- UKB_BL %>% left_join(UKB_beh, by = 'SubjectID')
    rm(UKB_beh)
    
    
    ############################# Save to reuse later #############################
    write.csv(UKB_BL, 'Data/Current_UKB/UKB_Behavioral/UKB_BehavioralMeasures_vs_PDPRS_Data_Nov2025_fullUKB.csv', row.names = FALSE)
  }
  UKB_BL <- read.csv('Data/Current_UKB/UKB_Behavioral/UKB_BehavioralMeasures_vs_PDPRS_Data_Nov2025_fullUKB.csv')
  UKB_BL <- UKB_BL %>% filter(Sex %in% sex_to_include)
  PD_PRScs_list <- colnames(UKB_BL)[grep("^PD_PRScs_GP2", colnames(UKB_BL))] # or only PD_PRScs_list <- 'PD_PRScs_GP2'
  
  
  
  
  
  
  
  
  
  
  ############################# Run the models for behavioral measures #############################
  # Create list of behavioral measures to analyze
  behavioral_measures <- c("Coffee_intake", "Fluid_intelligence_score", "BMI",
                           "Sleep_duration",
                           "Multiple_deprivation_index", "Smoking_PackYears", 
                           "Household_income", "Alcohol_usage_total", 
                           "Education_level", "Addiction_score_total",
                           "Constipation", "Orthostatic_hypotension", "Urinary_incontinence", 
                           "Hand_grip_strength", "Apathy_past2weeks", "Depressed_mood_past2weeks"  )
  
  # Loop over all PRS versions
  for (current_prs in PD_PRScs_list) {
    print(paste0("Running behavioral measure models for ", current_prs))
    
    # Create current output path
    current_output_folder <- paste0(main_output_folder, 'with_all_pathwayPRS_Behavioral_', model_sex, '/')
    dir.create(current_output_folder, recursive = FALSE, showWarnings = FALSE)
    
    # Create independent variable list, including the PRS we are looking at
    ind_variables <- paste0(current_prs, ' + ', Confounds)
    
    # Create dataframe to store results for all behavioral measures
    behavioral_results <- data.frame(
      Measure = behavioral_measures,
      p_value = NA,
      t_value = NA,
      r2 = NA,
      r2_adj = NA,
      effect_size = NA,
      CI_low = NA,
      CI_high = NA,
      SE = NA,
      N = NA
    )
    
    # Prepare UKB data for modeling
    UKB_BL_model <- UKB_BL
    UKB_BL_model$Sex <- as.factor(UKB_BL_model$Sex)
    UKB_BL_model$Age <- scale(UKB_BL_model$Age)[,1]
    UKB_BL_model$Age2 <- scale(UKB_BL_model$Age2)[,1]
    UKB_BL_model$Center_Number <- as.factor(UKB_BL_model$Center_Number)
    UKB_BL_model$Genotype_batch <- as.factor(UKB_BL_model$Genotype_batch)
    UKB_BL_model[[current_prs]] <- scale(UKB_BL_model[[current_prs]])[,1]
    
    # Run models for each behavioral measure
    for (i in 1:nrow(behavioral_results)) {
      measure_name <- behavioral_results$Measure[i]
      
      # Skip if the measure has no data
      if(!(measure_name %in% colnames(UKB_BL_model))) {
        behavioral_results$p_value[i] <- NA
        behavioral_results$t_value[i] <- NA
        behavioral_results$r2[i] <- NA
        behavioral_results$r2_adj[i] <- NA
        behavioral_results$effect_size[i] <- NA
        behavioral_results$CI_low[i] <- NA
        behavioral_results$CI_high[i] <- NA
        behavioral_results$SE[i] <- NA
        behavioral_results$N[i] <- 0
        next
      }
      
      # Create formula
      mdl_formula <- paste0(measure_name, ' ~ ', ind_variables)
      
      # Scale the behavioral measure
      UKB_BL_model[[measure_name]] <- scale(UKB_BL_model[[measure_name]])[,1]
      
      # Run linear model
      out_temp <- lm(formula = mdl_formula, data = UKB_BL_model)
      out <- summary(out_temp)
      out_ci <- confint(out_temp)
      
      # Store results
      behavioral_results$p_value[i] <- out$coefficients[2,4]
      behavioral_results$t_value[i] <- out$coefficients[2,3]
      behavioral_results$r2[i] <- out$r.squared
      behavioral_results$r2_adj[i] <- out$adj.r.squared
      behavioral_results$effect_size[i] <- out$coefficients[2,1]
      behavioral_results$CI_low[i] <- out_ci[2, 1]
      behavioral_results$CI_high[i] <- out_ci[2, 2]
      behavioral_results$SE[i] <- out$coefficients[2,2]
      behavioral_results$N[i] <- nrow(model.frame(out_temp))
    }
    
    # Add FDR correction
    behavioral_results$p_FDR <- p.adjust(behavioral_results$p_value, method = "fdr", n = nrow(behavioral_results))
    behavioral_results$isSignificant <- ifelse(behavioral_results$p_FDR <= 0.05, 1, NA)
    
    # Save results
    write.csv(behavioral_results, paste0(current_output_folder, 'Behavioral_Measures_vs_', current_prs, '_Regression_Results.csv'), row.names = FALSE)
  }
  
  # Create a combined dataframe with results from all PRS versions
  behavioral_results_all <- data.frame()
  
  for (current_prs in PD_PRScs_list) {
    # Read the results file
    file_path <- paste0(current_output_folder, 'Behavioral_Measures_vs_', current_prs, '_Regression_Results.csv')
    if (file.exists(file_path)) {
      results <- read.csv(file_path)
      results$PRS <- current_prs
      results <- results %>% relocate(PRS, .before = 'Measure')
      behavioral_results_all <- rbind(behavioral_results_all, results)
    }
  }
  
  # Save the combined results
  write.csv(behavioral_results_all, paste0(main_output_folder, 'Behavioral_Measures_All_PRS_Regression_Results_', model_sex, '.csv'), row.names = FALSE)
  
  
}






