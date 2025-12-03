library(tidyr)
library(dplyr)
library(lme4)
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
model_motion_options <- c('withMotion', 'noMotion')


# Create main output path
main_output_folder <- 'Outputs/PD_PRS_Regression_Results/with_all_pathwayPRS_Brain_GlobalMeasures/'
dir.create(main_output_folder, recursive = FALSE, showWarnings = FALSE)


#### Loop through all combinations
for (model_sex in model_sex_options) {
  for (model_motion in model_motion_options) {
    
    #### Automatic model selection (updated model)
    if (model_motion == 'withMotion') {
      
      if (model_sex == 'all') {
        Confounds <- paste0('Age + Age2 + Sex + Age:Sex + Age2:Sex + Days + Days2 + Motion_Average_fMRI + Bed_Position + Center_Number + Genotype_batch + ', PCAs) # version with motion
        sex_to_include <- c('Male','Female')
      } else if (model_sex == 'male') {
        sex_to_include <- 'Male' # Male (will be 1)
        Confounds <- paste0('Age + Age2 + Days + Days2 + Motion_Average_fMRI + Bed_Position + Center_Number + Genotype_batch + ', PCAs)
      } else if (model_sex == 'female') {
        sex_to_include <- 'Female' # Female (will be 0)
        Confounds <- paste0('Age + Age2 + Days + Days2 + Motion_Average_fMRI + Bed_Position + Center_Number + Genotype_batch + ', PCAs)
      } else { print('WRONG SEX SELECTED') }
      
    } else if (model_motion == 'noMotion') {
      
      if (model_sex == 'all') {
        Confounds <- paste0('Age + Age2 + Sex + Age:Sex + Age2:Sex + Days + Days2 + Center_Number + Genotype_batch + ', PCAs) # version with motion
        sex_to_include <- c('Male','Female')
      } else if (model_sex == 'male') {
        sex_to_include <- 'Male' # Male (will be 1)
        Confounds <- paste0('Age + Age2 + Days + Days2 + Center_Number + Genotype_batch + ', PCAs)
      } else if (model_sex == 'female') {
        sex_to_include <- 'Female' # Female (will be 0)
        Confounds <- paste0('Age + Age2 + Days + Days2 + Center_Number + Genotype_batch + ', PCAs)
      } else { print('WRONG SEX SELECTED') }
      
    } else { print('WRONG MOTION SELECTED') }
    
    print(Confounds)
    
    
    
    
    
    
    ############################# Reading Data #############################
    # UKB Current Data - First imaging session
    UKB_BL <- read.csv('Data/Current_UKB/UKB_Imaging1_withPCA_afterExclusion.csv', header = TRUE)
    UKB_BL <- UKB_BL %>% filter(Sex_31 %in% sex_to_include)
    # SA Data
    SA_DKT <- read.csv('Data/Surface_Area/SA_in_DKT_FreeSurfer_Yashar/UKB_SA_DKT_FreeSurfer_Yashar_QCed_3SDremoved.csv', header = TRUE)
    # CT date
    CT_DKT <- read.csv('Data/Cortical_Thickness/CT_in_DKT_FreeSurfer_Yashar/UKB_CT_DKT_FreeSurfer_Yashar_QCed_3SDremoved.csv', header = TRUE)
    # New PD-PRS
    PRS <- read.csv('Data/Genetic/PD_PRS/PD_PRS_allVersions_Final_mainPathways.csv') # has all PRS values
    UKB_BL <- UKB_BL %>% left_join(PRS, by = "SubjectID") %>% arrange(SubjectID)
    UKB_BL <- UKB_BL %>% filter(!is.na(PD_PRScs_Nalls))
    PRS <- PRS %>% filter(!is.na(PD_PRScs_Nalls))
    # Subcortical Volume Data
    Subcortical_Volumes <- read.csv('Data/Subcortical/SV_aseg_FreeSurfer_Yashar/UKB_Subcortical_Volumes_aseg_FreeSurfer_Yashar_QCed_3SDremoved.csv', header = TRUE)
    # White Matter Data
    # FA_tracts <- read.csv('Data/White_Matter/Tractwise_reExtracted/UKB_ORGTract_FA_reExtracted_QCed_3SDremoved.csv', header = TRUE)
    # MD_tracts <- read.csv('Data/White_Matter/Tractwise_reExtracted/UKB_ORGTract_MD_reExtracted_QCed_3SDremoved.csv', header = TRUE)
    FW_tracts <- read.csv('Data/White_Matter/Tractwise_reExtracted/UKB_ORGTract_FW_reExtracted_QCed_3SDremoved.csv', header = TRUE)
    # White Matter Global Data
    FA_global <- read.csv('Data/White_Matter/Global_Measures_BISON/UKB_Global_WM_FA_BISON_QCed_3SDremoved.csv', header = TRUE)
    MD_global <- read.csv('Data/White_Matter/Global_Measures_BISON/UKB_Global_WM_MD_BISON_QCed_3SDremoved.csv', header = TRUE)
    # T2start and SWI data
    SWI <- read.csv('Data/SWI/UKB_Magnetic_Susceptibility_3SDremoved.csv', header = TRUE)
    T2star <- read.csv('Data/SWI/UKB_T2star_3SDremoved.csv', header = TRUE)
    # # Atlases
    # DKT_Region_Names <- read.csv('Data/Atlases/CIVET_DKT_labels.csv')
    # subcortical_indices <- read.csv('Data/Subcortical/Subcortical_indices.csv', header = FALSE)
    # ORG_tracts_names <- read.csv('Data/White_Matter/ORG_tract_names_bilateral.csv')
    
    
    
    
    
    
    ############################# Keeping Confounds for all #############################
    UKB_data_final <- UKB_BL %>% select(c("SubjectID", "Sex_31", "Age.when.attended.assessment.centre_21003", "UK.Biobank.assessment.centre_54", "PCA_1", "PCA_2", "PCA_3", "PCA_4", "PCA_5",
                                          "PCA_6", "PCA_7", "PCA_8", "PCA_9", "PCA_10", "PCA_11", "PCA_12", "PCA_13", "PCA_14", "PCA_15"), starts_with('PD_PRS'))
    
    colnames(UKB_data_final)[3:4] <- c("Age", "cn")
    
    UKB_data_final <- UKB_data_final %>% mutate(Sex = case_when(Sex_31 == "Male" ~ 1,
                                                                Sex_31 == "Female" ~ 0)) %>% select(-c("Sex_31"))
    
    UKB_data_final$Age2 <- (UKB_data_final$Age)^2
    
    UKB_BL$Date.of.attending.assessment.centre_53 <- as.Date(UKB_BL$Date.of.attending.assessment.centre_53)
    UKB_data_final$Days <- as.numeric(difftime(UKB_BL$Date.of.attending.assessment.centre_53,min(UKB_BL$Date.of.attending.assessment.centre_53), units = 'days'))
    UKB_data_final$Days2 <- as.numeric(UKB_data_final$Days)^2
    
    
    # UKB_data_final$gn <- substr(UKB_BL$Genotype.measurement.batch_22000,start=1,stop=3)
    # UKB_data_final <- UKB_data_final %>%
    #   mutate(Genotype_batch = case_when(gn == "Bat" ~ 2,
    #                                     gn == "UKB" ~ 1)) %>%
    #   select(-c(gn)) #UKBiLEVEAX=1
    
    # Create a mapping of each unique batch to a number
    batches <- unique(UKB_BL$Genotype.measurement.batch_22000)
    batch_map <- setNames(1:length(batches), batches)
    # Apply the mapping to create the Genotype_batch variable
    UKB_data_final <- UKB_data_final %>% mutate(Genotype_batch = batch_map[UKB_BL$Genotype.measurement.batch_22000])
    
    UKB_data_final <- UKB_data_final %>% mutate(Center_Number = case_when(cn == "Cheadle (imaging)" ~ 1,
                                                                          cn == "Newcastle (imaging)" ~ 2,
                                                                          cn == "Reading (imaging)" ~ 3)) %>% select(-c(cn))
    
    UKB_data_final$Motion_rfmri <- UKB_BL$Mean.rfMRI.head.motion..averaged.across.space.and.time.points_25741
    UKB_data_final$Motion_tfmri <- UKB_BL$Mean.tfMRI.head.motion..averaged.across.space.and.time.points_25742
    UKB_data_final$Motion_Average_fMRI <- rowMeans(UKB_data_final[, c("Motion_rfmri", "Motion_tfmri")], na.rm = TRUE)
    UKB_data_final$Motion_Average_fMRI[UKB_data_final$Motion_Average_fMRI == 'NaN'] <- NA
    UKB_data_final$Bed_Position <- UKB_BL$Scanner.table.position_25759
    UKB_data_final$Head_Size <- UKB_BL$Volume.of.EstimatedTotalIntraCranial..whole.brain._26521
    UKB_data_final$Educational_Reform <- UKB_BL$Educational_Reform
    
    UKB_BL <- UKB_data_final
    rm(UKB_data_final)
    
    
    
    
    
    
    ############################# Calculate global measures by averaging or sum #############################
    total_SA_DKT <- data.frame(
      SubjectID = SA_DKT$SubjectID,
      Total_SA_DKT = rowSums(SA_DKT[, -which(names(SA_DKT) == "SubjectID")]))
    
    total_CT_DKT <- data.frame(
      SubjectID = CT_DKT$SubjectID,
      Total_CT_DKT = rowMeans(CT_DKT[, -which(names(CT_DKT) == "SubjectID")]))
    
    total_Subcortical_Volume <- data.frame(
      SubjectID = Subcortical_Volumes$SubjectID,
      Total_Subcortical_Volume = rowSums(Subcortical_Volumes[, -which(names(Subcortical_Volumes) == "SubjectID")]))
    
    total_FA <- FA_global
    colnames(total_FA) <- c('SubjectID', 'Total_FA')
    # total_FA <- data.frame(
    #   SubjectID = FA_tracts$SubjectID,
    #   Total_FA = rowMeans(FA_tracts[, -which(names(FA_tracts) == "SubjectID")]))
    
    total_MD <- MD_global
    colnames(total_MD) <- c('SubjectID', 'Total_MD')
    # total_MD <- data.frame(
    #   SubjectID = MD_tracts$SubjectID,
    #   Total_MD = rowMeans(MD_tracts[, -which(names(MD_tracts) == "SubjectID")]))
    
    total_FW <- data.frame(
      SubjectID = FW_tracts$SubjectID,
      Total_FW = rowMeans(FW_tracts[, -which(names(FW_tracts) == "SubjectID")]))
    
    total_SWI <- data.frame(
      SubjectID = SWI$SubjectID, 
      Total_SWI = rowMeans(SWI[, -which(names(SWI) == "SubjectID")]))
    
    total_T2Star <- data.frame(
      SubjectID = T2star$SubjectID,
      Total_T2Star = rowMeans(T2star[, -which(names(T2star) == "SubjectID")]))

    
    ### Merge values
    global_Measures <- total_SA_DKT %>% 
      full_join(total_CT_DKT, by = "SubjectID") %>%
      full_join(total_Subcortical_Volume, by = "SubjectID") %>%
      full_join(total_FA, by = "SubjectID") %>%
      full_join(total_MD, by = "SubjectID") %>%
      full_join(total_FW, by = "SubjectID") %>%
      full_join(total_SWI, by = "SubjectID") %>%
      full_join(total_T2Star, by = "SubjectID")
    
    global_Measures <- global_Measures %>% left_join(PRS, by = "SubjectID") %>% arrange(SubjectID) %>% filter(!is.na(PD_PRScs_Nalls))
    rm(PRS, SA_DKT, CT_DKT, Subcortical_Volumes, FA_tracts, MD_tracts, FW_tracts, SWI, T2star, FA_global, MD_global,
       total_SA_DKT, total_CT_DKT, total_Subcortical_Volume, total_FA, total_MD, total_FW, total_SWI, total_T2Star)
    
    
    
    
    
    
    ############################# Loop over all PRS versions and calculate LM for global measures #############################
    PD_PRScs_list <- colnames(UKB_BL)[grep("^PD_PRS", colnames(UKB_BL))]
    
    ### Run the models
    # Loop over all PRS versions
    for (current_prs in PD_PRScs_list) {
      print(paste0("Running global measure models for ", current_prs))
      
      # Create current output path (using the same folder structure you already have)
      current_output_folder <- paste0(main_output_folder, 'with_all_pathwayPRS_Brain_GlobalMeasures_', model_motion, '_', model_sex, '/')
      dir.create(current_output_folder, recursive = FALSE, showWarnings = FALSE)
      
      # Create independent variable list, including the PRS we are looking at
      ind_variables <- paste0(current_prs, ' + ', Confounds)
      # save these variables separatelt as well
      variables <- trimws(unlist(strsplit(ind_variables, "\\+"))) # Split variables by "+" and trim spaces
      variables <- variables[!grepl(":", variables)] # Remove interaction terms (those containing ":")
      
      # Create dataframe to store results for all global measures
      global_results <- data.frame(
        Measure = c("Total_SA_DKT", "Total_CT_DKT", 
                    "Total_Subcortical_Volume", 
                    "Total_FA", "Total_MD", "Total_FW", 
                    "Total_SWI", "Total_T2Star"),
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
      
      # Get common subjects between UKB_BL and global measures
      common_subjects <- intersect(UKB_BL[complete.cases(UKB_BL[, variables]), ]$SubjectID, global_Measures$SubjectID) # the first value allows only the rows that will be used in the model (i.e. have all the required variables)
      UKB_BL_common <- UKB_BL[UKB_BL$SubjectID %in% common_subjects, ] %>% arrange(SubjectID)
      global_Measures_common <- global_Measures[global_Measures$SubjectID %in% common_subjects, ] %>% arrange(SubjectID)
      
      # Check if subjects are in the same order
      identical(UKB_BL_common$SubjectID, global_Measures_common$SubjectID)
      
      # Prepare UKB data for modeling
      UKB_BL_common$Sex <- as.factor(UKB_BL_common$Sex)
      UKB_BL_common$Age <- scale(UKB_BL_common$Age)[,1]
      UKB_BL_common$Age2 <- scale(UKB_BL_common$Age2)[,1]
      UKB_BL_common$Days <- scale(UKB_BL_common$Days)[,1]
      UKB_BL_common$Days2 <- scale(UKB_BL_common$Days2)[,1]
      UKB_BL_common$Center_Number <- as.factor(UKB_BL_common$Center_Number)
      UKB_BL_common$Genotype_batch <- as.factor(UKB_BL_common$Genotype_batch)
      UKB_BL_common$Head_Size <- scale(UKB_BL_common$Head_Size)[,1]
      
      # Run models for each global measure
      for (i in 1:nrow(global_results)) {
        measure_name <- global_results$Measure[i]
        
        # Skip if the measure has no data
        if(!(measure_name %in% colnames(global_Measures_common))) {
          global_results$p_value[i] <- NA
          global_results$t_value[i] <- NA
          global_results$r2[i] <- NA
          global_results$r2_adj[i] <- NA
          global_results$effect_size[i] <- NA
          global_results$CI_low[i] <- NA
          global_results$CI_high[i] <- NA
          global_results$SE[i] <- NA
          global_results$N[i] <- 0
          next
        }
        
        # Create formula
        mdl_formula <- paste0(measure_name, ' ~ ', ind_variables)
        
        # Add the measure to UKB data
        UKB_BL_common[[measure_name]] <- scale(global_Measures_common[[measure_name]])[,1]
        
        # Run model
        out_temp <- lm(formula = mdl_formula, data = UKB_BL_common)
        out <- summary(out_temp)
        out_ci <- confint(out_temp)
        
        # Store results
        global_results$p_value[i] <- out$coefficients[2,4]
        global_results$t_value[i] <- out$coefficients[2,3]
        global_results$r2[i] <- out$r.squared
        global_results$r2_adj[i] <- out$adj.r.squared
        global_results$effect_size[i] <- out$coefficients[2,1]
        global_results$CI_low[i] <- out_ci[2, 1]
        global_results$CI_high[i] <- out_ci[2, 2]
        global_results$SE[i] <- out$coefficients[2,2]
        global_results$N[i] <- nrow(UKB_BL_common[!is.na(UKB_BL_common[[measure_name]]),])
      }
      
      # Add FDR correction
      global_results$p_FDR <- p.adjust(global_results$p_value, method = "fdr", n = nrow(global_results))
      global_results$isSignificant <- ifelse(global_results$p_FDR <= 0.05, 1, NA)
      
      # Save results
      write.csv(global_results, paste0(current_output_folder, 'Global_Measures_vs_', current_prs, '_Regression_Results.csv'), row.names = FALSE)
      
      # Clean up
      rm(common_subjects, UKB_BL_common, global_Measures_common)
    }
  
    
    
    
    ############################# Create a combined dataframe with results from all PRS versions - per motion and sex #############################
    global_results_all <- data.frame()
    
    for (current_prs in PD_PRScs_list) {
      # Read the results file
      file_path <- paste0(current_output_folder, 'Global_Measures_vs_', current_prs, '_Regression_Results.csv')
      if (file.exists(file_path)) {
        results <- read.csv(file_path)
        results$PRS <- current_prs
        results <- results %>% relocate(PRS, .before = 'Measure')
        global_results_all <- rbind(global_results_all, results)
      }
    }
    
    # Save the combined results
    write.csv(global_results_all, paste0(current_output_folder, 'Global_Measures_All_PRS_Regression_Results.csv'), row.names = FALSE)  
  }
}



