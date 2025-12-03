# This is the final script to use for Aim1 after Merge - Houman Azizi
library(tidyr)
library(dplyr)
library(lme4)
#library(ggplot2)

rm(list = ls())
cat("\014")

setwd("/Users/houmanazizi/Library/Mobile\ Documents/com~apple~CloudDocs/Education/GitHub/UKB_Analysis/")


############################# Create an empty dataframe for all samples across all combinations of runs #############################
sample_summary_all_merged <- data.frame(matrix(ncol = 12, nrow = 0))
colnames(sample_summary_all_merged) <- c('PRS', 'Measure', 'Cohort_sex', 'Motion_confound_status', 'Age_mean', 'Age_SD', 'N_women', 'N_men', 'N_total', 'Women_percentage', 'Sig_regions_N', 'Total_regions_N')





############################# Setting up Confounds + Go Through Each Version #############################
n_PCA <- 15
PCAs <- "PCA_1"
for (i in 2:n_PCA) {
  PCAs <- paste0(PCAs," + PCA_",i)
}

#### SET POSSIBLE OPTIONS HERE FOR THE MODEL
model_sex_options <- c('all', 'male', 'female')
model_motion_options <- c('withMotion', 'noMotion')




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
    # UKB Data - TIV and Euler Number
    UKB_additions <- read.csv('Data/Current_UKB/UKB_TIV_EulerNumber_QCed_FreeSurfer_Yashar.csv', header = TRUE)
    UKB_BL <- UKB_BL %>% left_join(UKB_additions, by = 'SubjectID')
    rm(UKB_additions)
    # SA Data
    SA_DKT <- read.csv('Data/Surface_Area/SA_in_DKT_FreeSurfer_Yashar/UKB_SA_DKT_FreeSurfer_Yashar_QCed_3SDremoved.csv', header = TRUE)
    # CT date
    CT_DKT <- read.csv('Data/Cortical_Thickness/CT_in_DKT_FreeSurfer_Yashar/UKB_CT_DKT_FreeSurfer_Yashar_QCed_3SDremoved.csv', header = TRUE)
    # New PD-PRS
    PRS <- read.csv('Data/Genetic/PD_PRS/PD_PRS_allVersions_Final_mainPathways.csv') # has main PRS values - or use 'Data/Genetic/PD_PRS/PD_PRS_allVersions_Final.csv'
    UKB_BL <- UKB_BL %>% left_join(PRS, by = "SubjectID") %>% arrange(SubjectID)
    UKB_BL <- UKB_BL %>% filter(!is.na(PD_PRScs_Nalls))
    rm(PRS)
    # Subcortical Volume Data
    Subcortical_Volumes <- read.csv('Data/Subcortical/SV_aseg_FreeSurfer_Yashar/UKB_Subcortical_Volumes_aseg_FreeSurfer_Yashar_QCed_3SDremoved.csv', header = TRUE)
    # White Matter Data
    FA_tracts <- read.csv('Data/White_Matter/Tractwise_reExtracted/UKB_ORGTract_FA_reExtracted_QCed_3SDremoved.csv', header = TRUE)
    MD_tracts <- read.csv('Data/White_Matter/Tractwise_reExtracted/UKB_ORGTract_MD_reExtracted_QCed_3SDremoved.csv', header = TRUE)
    FW_tracts <- read.csv('Data/White_Matter/Tractwise_reExtracted/UKB_ORGTract_FW_reExtracted_QCed_3SDremoved.csv', header = TRUE)
    # T2start and SWI data
    SWI <- read.csv('Data/SWI/UKB_Magnetic_Susceptibility_3SDremoved.csv', header = TRUE)
    T2star <- read.csv('Data/SWI/UKB_T2star_3SDremoved.csv', header = TRUE)
    # Atlases
    DKT_Region_Names <- read.csv('Data/Atlases/CIVET_DKT_labels.csv')
    subcortical_indices <- read.csv('Data/Subcortical/Subcortical_indices.csv', header = FALSE)
    ORG_tracts_names <- read.csv('Data/White_Matter/ORG_tract_names_bilateral.csv')
    
    
    
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
    
    UKB_data_final$Motion_Euler <- UKB_BL$Euler_Number
    UKB_data_final$Motion_rfmri <- UKB_BL$Mean.rfMRI.head.motion..averaged.across.space.and.time.points_25741
    UKB_data_final$Motion_tfmri <- UKB_BL$Mean.tfMRI.head.motion..averaged.across.space.and.time.points_25742
    UKB_data_final$Motion_Average_fMRI <- rowMeans(UKB_data_final[, c("Motion_rfmri", "Motion_tfmri")], na.rm = TRUE)
    UKB_data_final$Motion_Average_fMRI[UKB_data_final$Motion_Average_fMRI == 'NaN'] <- NA
    UKB_data_final$Bed_Position <- UKB_BL$Scanner.table.position_25759
    UKB_data_final$Head_Size <- UKB_BL$Volume.of.EstimatedTotalIntraCranial..whole.brain._26521 # or $TIV
    UKB_data_final$Educational_Reform <- UKB_BL$Educational_Reform
    
    UKB_BL <- UKB_data_final
    rm(UKB_data_final)
    
    
    
    
    
    
    
    
    
    ############################# Loop over all PRS versions and calculate region-wise #############################
    # list of all PD-PRS versions
    PD_PRScs_list <- colnames(UKB_BL)[grep("^PD_PRS", colnames(UKB_BL))]
    # # shorten PRS list - can use this to select
    # PD_PRScs_list <- c('PD_PRScs_Lang', 'PD_PRScs_Lang2025',
    #                    'PD_PRScs_Autophagy_direct_genes_4', 'PD_PRScs_Autophagy_direct_genes_4_exclude',
    #                    'PD_PRScs_Lysosomal_direct_genes_2', 'PD_PRScs_Lysosomal_direct_genes_2_exclude', 
    #                    'PD_PRScs_Mitochondrial_genes_1', 'PD_PRScs_Mitochondrial_genes_1_exclude',
    #                    'PD_PRScs_Mitochondrial_Lysosomal_Autophagy_direct_genes_9', 'PD_PRScs_Mitochondrial_Lysosomal_Autophagy_direct_genes_9_exclude',
    #                    'PD_PRSice_Lang',
    #                    'PD_PRSice_Autophagy_direct_genes_4', 'PD_PRSice_Autophagy_direct_genes_4_exclude', 
    #                    'PD_PRSice_Lysosomal_direct_genes_2', 'PD_PRSice_Lysosomal_direct_genes_2_exclude', 
    #                    'PD_PRSice_Mitochondrial_genes_1', 'PD_PRSice_Mitochondrial_genes_1_exclude',
    #                    'PD_PRSice_Mitochondrial_Lysosomal_Autophagy_direct_genes_9', 'PD_PRSice_Mitochondrial_Lysosomal_Autophagy_direct_genes_9_exclude')
    # 
    # Create an empty dataframe for all samples
    sample_summary_all <- data.frame(matrix(ncol = 10, nrow = 0))
    colnames(sample_summary_all) <- c('PRS', 'Measure', 'Age_mean', 'Age_SD', 'N_women', 'N_men', 'N_total', 'Women_percentage', 'Sig_regions_N', 'Total_regions_N')
    
    
    # looping over them
    dir.create(paste0('Outputs/PD_PRS_Regression_Results/with_all_pathwayPRS_', model_motion, '_', model_sex, '/'), recursive = FALSE, showWarnings = FALSE)
    for (current_prs in PD_PRScs_list) {
      print(current_prs)
      
      # create current output path
      current_output_folder <- paste0('Outputs/PD_PRS_Regression_Results/with_all_pathwayPRS_', model_motion, '_', model_sex, '/', current_prs, '/')
      dir.create(current_output_folder, recursive = FALSE, showWarnings = FALSE)
      
      # create independent variable list, including the PRS we are looking at
      ind_variables <- paste0(current_prs, ' + ', Confounds)
      # save these variables separatelt as well
      variables <- trimws(unlist(strsplit(ind_variables, "\\+"))) # Split variables by "+" and trim spaces
      variables <- variables[!grepl(":", variables)] # Remove interaction terms (those containing ":")
      
      # create an empty dataframe to keep N
      sample_summary <- data.frame(matrix(ncol = 9, nrow = 8))
      colnames(sample_summary) <- c('Measure', 'Age_mean', 'Age_SD', 'N_women', 'N_men', 'N_total', 'Women_percentage', 'Sig_regions_N', 'Total_regions_N')
      
      
      
      
      
      ############ Surface Area (SA ~ PRS + Confounds) - DKT FreeSurfer ############
      SA_results <- data.frame(matrix(NA, nrow = 62, ncol = 12))
      colnames(SA_results) <- c("Region", "Region_Klein", "Region_Number", "Region_Number_Klein", 'Region_Column_Name', "p_value", "p_FDR", "t_value", 'r2', 'effect_size', 'r2_adj', 'SE')
      SA_results$Region <- DKT_Region_Names$Region_name_UKB
      SA_results$Region_Klein <- DKT_Region_Names$Region_name_Klein
      SA_results$Region_Number <- DKT_Region_Names$CIVET_DKT_number
      SA_results$Region_Number_Klein <- DKT_Region_Names$Klein_original_number
      mdl_formula <- paste0('SA_region ~ ', ind_variables)
      # find common subjects between UKB_BL and SA data to match
      common_subjects <- intersect(UKB_BL[complete.cases(UKB_BL[, variables]), ]$SubjectID, SA_DKT$SubjectID) # the first value allows only the rows that will be used in the model (i.e. have all the required variables)
      UKB_BL_common <- UKB_BL[UKB_BL$SubjectID %in% common_subjects, ] %>% arrange(SubjectID)
      SA_DKT <- SA_DKT[SA_DKT$SubjectID %in% common_subjects, ] %>% arrange(SubjectID)
      
      print(paste0('Age_mean=',mean(UKB_BL_common$Age), '   Age_SD=',sd(UKB_BL_common$Age), '   N_women=',sum(UKB_BL_common$Sex == 0) , '   Women%=',(sum(UKB_BL_common$Sex == 0)/nrow(UKB_BL_common))*100))
      current_common_subjects_age_mean <- mean(UKB_BL_common$Age)
      current_common_subjects_age_sd <- sd(UKB_BL_common$Age)
      print(paste0( 'Age range = ', min(UKB_BL_common$Age), ' to ', max(UKB_BL_common$Age) ))
      
      UKB_BL_common$Sex <- as.factor(UKB_BL_common$Sex)
      UKB_BL_common$Age <- scale(UKB_BL_common$Age)[,1]
      UKB_BL_common$Age2 <- scale(UKB_BL_common$Age2)[,1]
      UKB_BL_common$Days <- scale(UKB_BL_common$Days)[,1]
      UKB_BL_common$Days2 <- scale(UKB_BL_common$Days2)[,1]
      UKB_BL_common$Center_Number <- as.factor(UKB_BL_common$Center_Number)
      UKB_BL_common$Genotype_batch <- as.factor(UKB_BL_common$Genotype_batch)
      UKB_BL_common$Head_Size <- scale(UKB_BL_common$Head_Size)[,1]
      identical(UKB_BL_common$SubjectID, SA_DKT$SubjectID)
      # run the model
      for (i in 1:nrow(SA_results)) {
        
        UKB_BL_common$SA_region <- scale(SA_DKT[,(i+1)])[,1]
        SA_results$Region_Column_Name[i] <- colnames(SA_DKT)[i+1]
        
        out_temp <-  lm(formula = mdl_formula, data = UKB_BL_common)
        out <- summary(out_temp)
        
        SA_results$p_value[i] <- out$coefficients[2,4]
        SA_results$t_value[i] <- out$coefficients[2,3]
        SA_results$r2[i] <- out$r.squared
        SA_results$effect_size[i] <- out$coefficients[2,1]
        SA_results$r2_adj[i] <- out$adj.r.squared
        SA_results$SE[i] <- out$coefficients[2,2]
      }
      SA_results$p_FDR <- p.adjust(SA_results$p_value, method = "fdr", n = nrow(SA_results))
      sum(SA_results$p_FDR < 0.05, na.rm = TRUE)
      SA_results$isSignificant <- ifelse(SA_results$p_FDR <= 0.05, 1, NA)
      # save and remove extra variables
      sample_summary[1, ] <- c('SA_DKT', current_common_subjects_age_mean, current_common_subjects_age_sd, sum(UKB_BL_common$Sex == 0), sum(UKB_BL_common$Sex == 1), nrow(UKB_BL_common), (sum(UKB_BL_common$Sex == 0)/nrow(UKB_BL_common))*100, sum(SA_results$p_FDR < 0.05, na.rm = TRUE), nrow(SA_results))
      write.csv(SA_results, paste0(current_output_folder, 'SA_DKT_vs_', current_prs, '_Regression_Results.csv'), row.names = FALSE)
      rm(common_subjects, UKB_BL_common, current_common_subjects_age_mean, current_common_subjects_age_sd, SA_results)
      
      
      
      
      
      
      
      ############ Cortical Thickness (CT ~ PRS + Confounds) - DKT FreeSurfer ############
      CT_results <- data.frame(matrix(NA, nrow = 62, ncol = 12))
      colnames(CT_results) <- c("Region", "Region_Klein", "Region_Number", "Region_Number_Klein", 'Region_Column_Name', "p_value", "p_FDR", "t_value", 'r2', 'effect_size', 'r2_adj', 'SE')
      CT_results$Region <- DKT_Region_Names$Region_name_UKB
      CT_results$Region_Klein <- DKT_Region_Names$Region_name_Klein
      CT_results$Region_Number <- DKT_Region_Names$CIVET_DKT_number
      CT_results$Region_Number_Klein <- DKT_Region_Names$Klein_original_number
      mdl_formula <- paste0('CT_region ~ ', ind_variables)
      # find common subjects between UKB_BL and CT data to match
      common_subjects <- intersect(UKB_BL[complete.cases(UKB_BL[, variables]), ]$SubjectID, CT_DKT$SubjectID)
      UKB_BL_common <- UKB_BL[UKB_BL$SubjectID %in% common_subjects, ] %>% arrange(SubjectID)
      CT_DKT <- CT_DKT[CT_DKT$SubjectID %in% common_subjects, ] %>% arrange(SubjectID)
      
      print(paste0('Age_mean=',mean(UKB_BL_common$Age), '   Age_SD=',sd(UKB_BL_common$Age), '   N_women=',sum(UKB_BL_common$Sex == 0) , '   Women%=',(sum(UKB_BL_common$Sex == 0)/nrow(UKB_BL_common))*100))
      current_common_subjects_age_mean <- mean(UKB_BL_common$Age)
      current_common_subjects_age_sd <- sd(UKB_BL_common$Age)
      print(paste0( 'Age range = ', min(UKB_BL_common$Age), ' to ', max(UKB_BL_common$Age) ))
      
      UKB_BL_common$Sex <- as.factor(UKB_BL_common$Sex)
      UKB_BL_common$Age <- scale(UKB_BL_common$Age)[,1]
      UKB_BL_common$Age2 <- scale(UKB_BL_common$Age2)[,1]
      UKB_BL_common$Days <- scale(UKB_BL_common$Days)[,1]
      UKB_BL_common$Days2 <- scale(UKB_BL_common$Days2)[,1]
      UKB_BL_common$Center_Number <- as.factor(UKB_BL_common$Center_Number)
      UKB_BL_common$Genotype_batch <- as.factor(UKB_BL_common$Genotype_batch)
      UKB_BL_common$Head_Size <- scale(UKB_BL_common$Head_Size)[,1]
      identical(UKB_BL_common$SubjectID, CT_DKT$SubjectID)
      # run the model
      for (i in 1:nrow(CT_results)) {
        
        UKB_BL_common$CT_region <- scale(CT_DKT[,(i+1)])[,1]
        CT_results$Region_Column_Name[i] <- colnames(CT_DKT)[i+1]
        
        out_temp <-  lm(formula = mdl_formula, data = UKB_BL_common)
        out <- summary(out_temp)
        
        CT_results$p_value[i] <- out$coefficients[2,4]
        CT_results$t_value[i] <- out$coefficients[2,3]
        CT_results$r2[i] <- out$r.squared
        CT_results$effect_size[i] <- out$coefficients[2,1]
        CT_results$r2_adj[i] <- out$adj.r.squared
        CT_results$SE[i] <- out$coefficients[2,2]
      }
      CT_results$p_FDR <- p.adjust(CT_results$p_value, method = "fdr", n = nrow(CT_results))
      sum(CT_results$p_FDR < 0.05, na.rm = TRUE)
      CT_results$isSignificant <- ifelse(CT_results$p_FDR <= 0.05, 1, NA)
      # save and remove extra variables
      sample_summary[2, ] <- c('CT_DKT', current_common_subjects_age_mean, current_common_subjects_age_sd, sum(UKB_BL_common$Sex == 0), sum(UKB_BL_common$Sex == 1), nrow(UKB_BL_common), (sum(UKB_BL_common$Sex == 0)/nrow(UKB_BL_common))*100, sum(CT_results$p_FDR < 0.05, na.rm = TRUE), nrow(CT_results))
      write.csv(CT_results, paste0(current_output_folder, 'CT_DKT_vs_', current_prs, '_Regression_Results.csv'), row.names = FALSE)
      rm(common_subjects, UKB_BL_common, current_common_subjects_age_mean, current_common_subjects_age_sd, CT_results)
      
      
      
      
      
      
      
      ############ Subcortical Volume (Sub_Vol ~ PRS + Confounds) - aseg ############
      Sub_results <- data.frame(matrix(NA, nrow = 14, ncol = 9))
      colnames(Sub_results) <- c("Region", "Region_Number", "p_value", "p_FDR", "t_value", 'r2', 'effect_size', 'r2_adj', 'SE')
      Sub_results$Region <- subcortical_indices$V1
      Sub_results$Region_Number <- subcortical_indices$V2
      mdl_formula <- paste0('Subcortical_Volume ~ ', ind_variables)
      # find common subjects between UKB_BL and subcortical data to match
      common_subjects <- intersect(UKB_BL[complete.cases(UKB_BL[, variables]), ]$SubjectID, Subcortical_Volumes$SubjectID)
      UKB_BL_common <- UKB_BL[UKB_BL$SubjectID %in% common_subjects, ] %>% arrange(SubjectID)
      Subcortical_Volumes <- Subcortical_Volumes[Subcortical_Volumes$SubjectID %in% common_subjects, ] %>% arrange(SubjectID)
      
      print(paste0('Age_mean=',mean(UKB_BL_common$Age), '   Age_SD=',sd(UKB_BL_common$Age), '   N_women=',sum(UKB_BL_common$Sex == 0) , '   Women%=',(sum(UKB_BL_common$Sex == 0)/nrow(UKB_BL_common))*100))
      current_common_subjects_age_mean <- mean(UKB_BL_common$Age)
      current_common_subjects_age_sd <- sd(UKB_BL_common$Age)
      print(paste0( 'Age range = ', min(UKB_BL_common$Age), ' to ', max(UKB_BL_common$Age) ))
      
      UKB_BL_common$Sex <- as.factor(UKB_BL_common$Sex)
      UKB_BL_common$Age <- scale(UKB_BL_common$Age)[,1]
      UKB_BL_common$Age2 <- scale(UKB_BL_common$Age2)[,1]
      UKB_BL_common$Days <- scale(UKB_BL_common$Days)[,1]
      UKB_BL_common$Days2 <- scale(UKB_BL_common$Days2)[,1]
      UKB_BL_common$Center_Number <- as.factor(UKB_BL_common$Center_Number)
      UKB_BL_common$Genotype_batch <- as.factor(UKB_BL_common$Genotype_batch)
      UKB_BL_common$Head_Size <- scale(UKB_BL_common$Head_Size)[,1]
      identical(UKB_BL_common$SubjectID, Subcortical_Volumes$SubjectID)
      # run the model
      for (i in 1:nrow(Sub_results)) {
        UKB_BL_common$Subcortical_Volume <- scale(Subcortical_Volumes[,i+1])[,1]
        
        out_temp <-  lm(formula = mdl_formula, data = UKB_BL_common)
        out <- summary(out_temp)
        Sub_results$t_value[which(Sub_results$Region==colnames(Subcortical_Volumes)[i+1])] <- out$coefficients[2,3]
        Sub_results$p_value[which(Sub_results$Region==colnames(Subcortical_Volumes)[i+1])] <- out$coefficients[2,4]
        Sub_results$r2[which(Sub_results$Region==colnames(Subcortical_Volumes)[i+1])] <- out$r.squared
        Sub_results$effect_size[which(Sub_results$Region==colnames(Subcortical_Volumes)[i+1])] <- out$coefficients[2,1]
        Sub_results$r2_adj[which(Sub_results$Region==colnames(Subcortical_Volumes)[i+1])] <- out$adj.r.squared
        Sub_results$SE[which(Sub_results$Region==colnames(Subcortical_Volumes)[i+1])] <- out$coefficients[2,2]
      }
      Sub_results$p_FDR <- p.adjust(Sub_results$p_value, method = "fdr", n = nrow(Sub_results))
      sum(Sub_results$p_FDR <= 0.05, na.rm = TRUE)
      Sub_results$isSignificant <- ifelse(Sub_results$p_FDR <= 0.05, 1, NA)
      # save and remove extra variables
      sample_summary[3, ] <- c('Subcortical_volume', current_common_subjects_age_mean, current_common_subjects_age_sd, sum(UKB_BL_common$Sex == 0), sum(UKB_BL_common$Sex == 1), nrow(UKB_BL_common), (sum(UKB_BL_common$Sex == 0)/nrow(UKB_BL_common))*100, sum(Sub_results$p_FDR < 0.05, na.rm = TRUE), nrow(Sub_results))
      write.csv(Sub_results, paste0(current_output_folder, 'SubcorticaVolume_vs_', current_prs, '_Regression_Results.csv'), row.names = FALSE)
      rm(common_subjects, UKB_BL_common, current_common_subjects_age_mean, current_common_subjects_age_sd, Sub_results)
      
      
      
      
      
      
      
      
      ############ FA (FA ~ PRS + Confounds) - reExtrated ############
      FA_results <- data.frame(matrix(NA, nrow = ncol(FA_tracts)-1, ncol = 9))
      colnames(FA_results) <- c("Tract", "Tract_file_name", "p_value", "p_FDR", "t_value", 'r2', 'effect_size', 'r2_adj', 'SE')
      FA_results$Tract <- ORG_tracts_names$Tract_real_name
      FA_results$Tract_file_name <- ORG_tracts_names$Tract_file_name
      mdl_formula <- paste0('FA_tract ~ ', ind_variables)
      # find common subjects between UKB_BL and FA data to match
      common_subjects <- intersect(UKB_BL[complete.cases(UKB_BL[, variables]), ]$SubjectID, FA_tracts$SubjectID)
      UKB_BL_common <- UKB_BL[UKB_BL$SubjectID %in% common_subjects, ] %>% arrange(SubjectID)
      FA_tracts <- FA_tracts[FA_tracts$SubjectID %in% common_subjects, ] %>% arrange(SubjectID)
      
      print(paste0('Age_mean=',mean(UKB_BL_common$Age), '   Age_SD=',sd(UKB_BL_common$Age), '   N_women=',sum(UKB_BL_common$Sex == 0) , '   Women%=',(sum(UKB_BL_common$Sex == 0)/nrow(UKB_BL_common))*100))
      current_common_subjects_age_mean <- mean(UKB_BL_common$Age)
      current_common_subjects_age_sd <- sd(UKB_BL_common$Age)
      print(paste0( 'Age range = ', min(UKB_BL_common$Age), ' to ', max(UKB_BL_common$Age) ))
      
      UKB_BL_common$Sex <- as.factor(UKB_BL_common$Sex)
      UKB_BL_common$Age <- scale(UKB_BL_common$Age)[,1]
      UKB_BL_common$Age2 <- scale(UKB_BL_common$Age2)[,1]
      UKB_BL_common$Days <- scale(UKB_BL_common$Days)[,1]
      UKB_BL_common$Days2 <- scale(UKB_BL_common$Days2)[,1]
      UKB_BL_common$Center_Number <- as.factor(UKB_BL_common$Center_Number)
      UKB_BL_common$Genotype_batch <- as.factor(UKB_BL_common$Genotype_batch)
      UKB_BL_common$Head_Size <- scale(UKB_BL_common$Head_Size)[,1]
      identical(UKB_BL_common$SubjectID, FA_tracts$SubjectID)
      # run the model
      for (i in 1:nrow(FA_results)) {
        UKB_BL_common$FA_tract <- scale(FA_tracts[,i+1])[,1]
        
        out_temp <-  lm(formula = mdl_formula, data = UKB_BL_common)
        out <- summary(out_temp)
        
        FA_results$t_value[which(FA_results$Tract==colnames(FA_tracts)[i+1])] <- out$coefficients[2,3]
        FA_results$p_value[which(FA_results$Tract==colnames(FA_tracts)[i+1])] <- out$coefficients[2,4]
        FA_results$r2[which(FA_results$Tract==colnames(FA_tracts)[i+1])] <- out$r.squared
        FA_results$effect_size[which(FA_results$Tract==colnames(FA_tracts)[i+1])] <- out$coefficients[2,1]
        FA_results$r2_adj[which(FA_results$Tract==colnames(FA_tracts)[i+1])] <- out$adj.r.squared
        FA_results$SE[which(FA_results$Tract==colnames(FA_tracts)[i+1])] <- out$coefficients[2,2]
      }
      FA_results$p_FDR <- p.adjust(FA_results$p_value, method = "fdr", n = nrow(FA_results))
      sum(FA_results$p_FDR <= 0.05, na.rm = TRUE)
      FA_results$isSignificant <- ifelse(FA_results$p_FDR <= 0.05, 1, NA)
      # save and remove extra variables
      sample_summary[4, ] <- c('FA', current_common_subjects_age_mean, current_common_subjects_age_sd, sum(UKB_BL_common$Sex == 0), sum(UKB_BL_common$Sex == 1), nrow(UKB_BL_common), (sum(UKB_BL_common$Sex == 0)/nrow(UKB_BL_common))*100, sum(FA_results$p_FDR < 0.05, na.rm = TRUE), nrow(FA_results))
      write.csv(FA_results, paste0(current_output_folder, 'FA_Tracts_vs_', current_prs, '_Regression_Results.csv'), row.names = FALSE)
      rm(common_subjects, UKB_BL_common, current_common_subjects_age_mean, current_common_subjects_age_sd, FA_results)
      
      
      
      
      
      
      
      
      ############ MD (MD ~ PRS + Confounds) - reExtracted ############
      MD_results <- data.frame(matrix(NA, nrow = ncol(MD_tracts)-1, ncol = 9))
      colnames(MD_results) <- c("Tract", "Tract_file_name", "p_value", "p_FDR", "t_value", 'r2', 'effect_size', 'r2_adj', 'SE')
      MD_results$Tract <- ORG_tracts_names$Tract_real_name
      MD_results$Tract_file_name <- ORG_tracts_names$Tract_file_name
      mdl_formula <- paste0('MD_tract ~ ', ind_variables)
      # find common subjects between UKB_BL and MD data to match
      common_subjects <- intersect(UKB_BL[complete.cases(UKB_BL[, variables]), ]$SubjectID, MD_tracts$SubjectID)
      UKB_BL_common <- UKB_BL[UKB_BL$SubjectID %in% common_subjects, ] %>% arrange(SubjectID)
      MD_tracts <- MD_tracts[MD_tracts$SubjectID %in% common_subjects, ] %>% arrange(SubjectID)
      
      print(paste0('Age_mean=',mean(UKB_BL_common$Age), '   Age_SD=',sd(UKB_BL_common$Age), '   N_women=',sum(UKB_BL_common$Sex == 0) , '   Women%=',(sum(UKB_BL_common$Sex == 0)/nrow(UKB_BL_common))*100))
      current_common_subjects_age_mean <- mean(UKB_BL_common$Age)
      current_common_subjects_age_sd <- sd(UKB_BL_common$Age)
      print(paste0( 'Age range = ', min(UKB_BL_common$Age), ' to ', max(UKB_BL_common$Age) ))
      
      UKB_BL_common$Sex <- as.factor(UKB_BL_common$Sex)
      UKB_BL_common$Age <- scale(UKB_BL_common$Age)[,1]
      UKB_BL_common$Age2 <- scale(UKB_BL_common$Age2)[,1]
      UKB_BL_common$Days <- scale(UKB_BL_common$Days)[,1]
      UKB_BL_common$Days2 <- scale(UKB_BL_common$Days2)[,1]
      UKB_BL_common$Center_Number <- as.factor(UKB_BL_common$Center_Number)
      UKB_BL_common$Genotype_batch <- as.factor(UKB_BL_common$Genotype_batch)
      UKB_BL_common$Head_Size <- scale(UKB_BL_common$Head_Size)[,1]
      identical(UKB_BL_common$SubjectID, MD_tracts$SubjectID)
      # run the model
      for (i in 1:nrow(MD_results)) {
        UKB_BL_common$MD_tract <- scale(MD_tracts[,i+1])[,1]
        
        out_temp <-  lm(formula = mdl_formula, data = UKB_BL_common)
        out <- summary(out_temp)
        MD_results$t_value[which(MD_results$Tract==colnames(MD_tracts)[i+1])] <- out$coefficients[2,3]
        MD_results$p_value[which(MD_results$Tract==colnames(MD_tracts)[i+1])] <- out$coefficients[2,4]
        MD_results$r2[which(MD_results$Tract==colnames(MD_tracts)[i+1])] <- out$r.squared
        MD_results$effect_size[which(MD_results$Tract==colnames(MD_tracts)[i+1])] <- out$coefficients[2,1]
        MD_results$r2_adj[which(MD_results$Tract==colnames(MD_tracts)[i+1])] <- out$adj.r.squared
        MD_results$SE[which(MD_results$Tract==colnames(MD_tracts)[i+1])] <- out$coefficients[2,2]
      }
      MD_results$p_FDR <- p.adjust(MD_results$p_value, method = "fdr", n = nrow(MD_results))
      sum(MD_results$p_FDR <= 0.05, na.rm = TRUE)
      MD_results$isSignificant <- ifelse(MD_results$p_FDR <= 0.05, 1, NA)
      # save and remove extra variables
      sample_summary[5, ] <- c('MD', current_common_subjects_age_mean, current_common_subjects_age_sd, sum(UKB_BL_common$Sex == 0), sum(UKB_BL_common$Sex == 1), nrow(UKB_BL_common), (sum(UKB_BL_common$Sex == 0)/nrow(UKB_BL_common))*100, sum(MD_results$p_FDR < 0.05, na.rm = TRUE), nrow(MD_results))
      write.csv(MD_results, paste0(current_output_folder, 'MD_Tracts_vs_', current_prs, '_Regression_Results.csv'), row.names = FALSE)
      rm(common_subjects, UKB_BL_common, current_common_subjects_age_mean, current_common_subjects_age_sd, MD_results)
      
      
      
      
      
      
      
      
      ############ FW (FW ~ PRS + Confounds) - reExtrated ############ ##### COMPLETEEEEEEEEEEE THISSSSSSSSSSSSSSSS
      FW_results <- data.frame(matrix(NA, nrow = ncol(FW_tracts)-1, ncol = 9))
      colnames(FW_results) <- c("Tract", "Tract_file_name", "p_value", "p_FDR", "t_value", 'r2', 'effect_size', 'r2_adj', 'SE')
      FW_results$Tract <- ORG_tracts_names$Tract_real_name
      FW_results$Tract_file_name <- ORG_tracts_names$Tract_file_name
      mdl_formula <- paste0('FW_tract ~ ', ind_variables)
      # find common subjects between UKB_BL and FW data to match
      common_subjects <- intersect(UKB_BL[complete.cases(UKB_BL[, variables]), ]$SubjectID, FW_tracts$SubjectID)
      UKB_BL_common <- UKB_BL[UKB_BL$SubjectID %in% common_subjects, ] %>% arrange(SubjectID)
      FW_tracts <- FW_tracts[FW_tracts$SubjectID %in% common_subjects, ] %>% arrange(SubjectID)
      
      print(paste0('Age_mean=',mean(UKB_BL_common$Age), '   Age_SD=',sd(UKB_BL_common$Age), '   N_women=',sum(UKB_BL_common$Sex == 0) , '   Women%=',(sum(UKB_BL_common$Sex == 0)/nrow(UKB_BL_common))*100))
      current_common_subjects_age_mean <- mean(UKB_BL_common$Age)
      current_common_subjects_age_sd <- sd(UKB_BL_common$Age)
      
      UKB_BL_common$Sex <- as.factor(UKB_BL_common$Sex)
      UKB_BL_common$Age <- scale(UKB_BL_common$Age)[,1]
      UKB_BL_common$Age2 <- scale(UKB_BL_common$Age2)[,1]
      UKB_BL_common$Days <- scale(UKB_BL_common$Days)[,1]
      UKB_BL_common$Days2 <- scale(UKB_BL_common$Days2)[,1]
      UKB_BL_common$Center_Number <- as.factor(UKB_BL_common$Center_Number)
      UKB_BL_common$Genotype_batch <- as.factor(UKB_BL_common$Genotype_batch)
      UKB_BL_common$Head_Size <- scale(UKB_BL_common$Head_Size)[,1]
      identical(UKB_BL_common$SubjectID, FW_tracts$SubjectID)
      # run the model
      for (i in 1:nrow(FW_results)) {
        UKB_BL_common$FW_tract <- scale(FW_tracts[,i+1])[,1]
        
        out_temp <-  lm(formula = mdl_formula, data = UKB_BL_common)
        out <- summary(out_temp)
        
        FW_results$t_value[which(FW_results$Tract==colnames(FW_tracts)[i+1])] <- out$coefficients[2,3]
        FW_results$p_value[which(FW_results$Tract==colnames(FW_tracts)[i+1])] <- out$coefficients[2,4]
        FW_results$r2[which(FW_results$Tract==colnames(FW_tracts)[i+1])] <- out$r.squared
        FW_results$effect_size[which(FW_results$Tract==colnames(FW_tracts)[i+1])] <- out$coefficients[2,1]
        FW_results$r2_adj[which(FW_results$Tract==colnames(FW_tracts)[i+1])] <- out$adj.r.squared
        FW_results$SE[which(FW_results$Tract==colnames(FW_tracts)[i+1])] <- out$coefficients[2,2]
      }
      FW_results$p_FDR <- p.adjust(FW_results$p_value, method = "fdr", n = nrow(FW_results))
      sum(FW_results$p_FDR <= 0.05, na.rm = TRUE)
      FW_results$isSignificant <- ifelse(FW_results$p_FDR <= 0.05, 1, NA)
      # save and remove extra variables
      sample_summary[6, ] <- c('FW', current_common_subjects_age_mean, current_common_subjects_age_sd, sum(UKB_BL_common$Sex == 0), sum(UKB_BL_common$Sex == 1), nrow(UKB_BL_common), (sum(UKB_BL_common$Sex == 0)/nrow(UKB_BL_common))*100, sum(FW_results$p_FDR < 0.05, na.rm = TRUE), nrow(FW_results))
      write.csv(FW_results, paste0(current_output_folder, 'FW_Tracts_vs_', current_prs, '_Regression_Results.csv'), row.names = FALSE)
      rm(common_subjects, UKB_BL_common, current_common_subjects_age_mean, current_common_subjects_age_sd, FW_results)
      
      
      
      
      
      
      
      
      ############ SWI (SWI ~ PRS + Confounds) ############
      SWI_results <- data.frame(matrix(NA, nrow = ncol(SWI)-1, ncol = 8))
      colnames(SWI_results) <- c("Region", "p_value", "p_FDR", "t_value", 'r2', 'effect_size', 'r2_adj', 'SE')
      SWI_results$Region <- colnames(SWI)[-1]
      mdl_formula <- paste0('SWI_region ~ ', ind_variables)
      # find common subjects between UKB_BL and SWI data to match
      common_subjects <- intersect(UKB_BL[complete.cases(UKB_BL[, variables]), ]$SubjectID, SWI$SubjectID)
      UKB_BL_common <- UKB_BL[UKB_BL$SubjectID %in% common_subjects, ] %>% arrange(SubjectID)
      SWI <- SWI[SWI$SubjectID %in% common_subjects, ] %>% arrange(SubjectID)
      
      current_common_subjects_age_mean <- mean(UKB_BL_common$Age)
      current_common_subjects_age_sd <- sd(UKB_BL_common$Age)
      print(paste0( 'Age range = ', min(UKB_BL_common$Age), ' to ', max(UKB_BL_common$Age) ))
      
      UKB_BL_common$Sex <- as.factor(UKB_BL_common$Sex)
      UKB_BL_common$Age <- scale(UKB_BL_common$Age)[,1]
      UKB_BL_common$Age2 <- scale(UKB_BL_common$Age2)[,1]
      UKB_BL_common$Days <- scale(UKB_BL_common$Days)[,1]
      UKB_BL_common$Days2 <- scale(UKB_BL_common$Days2)[,1]
      UKB_BL_common$Center_Number <- as.factor(UKB_BL_common$Center_Number)
      UKB_BL_common$Genotype_batch <- as.factor(UKB_BL_common$Genotype_batch)
      UKB_BL_common$Head_Size <- scale(UKB_BL_common$Head_Size)[,1]
      identical(UKB_BL_common$SubjectID, SWI$SubjectID)
      # run the model
      for (i in 1:nrow(SWI_results)) {
        UKB_BL_common$SWI_region <- scale(SWI[,i+1])[,1]
        
        out_temp <-  lm(formula = mdl_formula, data = UKB_BL_common)
        out <- summary(out_temp)
        SWI_results$t_value[which(SWI_results$Region==colnames(SWI)[i+1])] <- out$coefficients[2,3]
        SWI_results$p_value[which(SWI_results$Region==colnames(SWI)[i+1])] <- out$coefficients[2,4]
        SWI_results$r2[which(SWI_results$Region==colnames(SWI)[i+1])] <- out$r.squared
        SWI_results$effect_size[which(SWI_results$Region==colnames(SWI)[i+1])] <- out$coefficients[2,1]
        SWI_results$r2_adj[which(SWI_results$Region==colnames(SWI)[i+1])] <- out$adj.r.squared
        SWI_results$SE[which(SWI_results$Region==colnames(SWI)[i+1])] <- out$coefficients[2,2]
      }
      SWI_results$p_FDR <- p.adjust(SWI_results$p_value, method = "fdr", n = nrow(SWI_results))
      sum(SWI_results$p_FDR <= 0.05, na.rm = TRUE)
      SWI_results$isSignificant <- ifelse(SWI_results$p_FDR <= 0.05, 1, NA)
      SWI_results <- SWI_results %>% arrange(Region)
      # save and remove extra variables
      sample_summary[7, ] <- c('SWI', current_common_subjects_age_mean, current_common_subjects_age_sd, sum(UKB_BL_common$Sex == 0), sum(UKB_BL_common$Sex == 1), nrow(UKB_BL_common), (sum(UKB_BL_common$Sex == 0)/nrow(UKB_BL_common))*100, sum(SWI_results$p_FDR < 0.05, na.rm = TRUE), nrow(SWI_results))
      write.csv(SWI_results, paste0(current_output_folder, 'SWI_vs_', current_prs, '_Regression_Results.csv'), row.names = FALSE)
      rm(common_subjects, UKB_BL_common, current_common_subjects_age_mean, current_common_subjects_age_sd, SWI_results)
      
      
      
      
      
      ############ T2star (T2star ~ PRS + Confounds) ############
      T2star_results <- data.frame(matrix(NA, nrow = ncol(T2star)-1, ncol = 8))
      colnames(T2star_results) <- c("Region", "p_value", "p_FDR", "t_value", 'r2', 'effect_size', 'r2_adj', 'SE')
      T2star_results$Region <- colnames(T2star)[-1]
      mdl_formula <- paste0('T2star_region ~ ', ind_variables)
      # find common subjects between UKB_BL and T2star data to match
      common_subjects <- intersect(UKB_BL[complete.cases(UKB_BL[, variables]), ]$SubjectID, T2star$SubjectID)
      UKB_BL_common <- UKB_BL[UKB_BL$SubjectID %in% common_subjects, ] %>% arrange(SubjectID)
      T2star <- T2star[T2star$SubjectID %in% common_subjects, ] %>% arrange(SubjectID)
      
      current_common_subjects_age_mean <- mean(UKB_BL_common$Age)
      current_common_subjects_age_sd <- sd(UKB_BL_common$Age)
      
      UKB_BL_common$Sex <- as.factor(UKB_BL_common$Sex)
      UKB_BL_common$Age <- scale(UKB_BL_common$Age)[,1]
      UKB_BL_common$Age2 <- scale(UKB_BL_common$Age2)[,1]
      UKB_BL_common$Days <- scale(UKB_BL_common$Days)[,1]
      UKB_BL_common$Days2 <- scale(UKB_BL_common$Days2)[,1]
      UKB_BL_common$Center_Number <- as.factor(UKB_BL_common$Center_Number)
      UKB_BL_common$Genotype_batch <- as.factor(UKB_BL_common$Genotype_batch)
      UKB_BL_common$Head_Size <- scale(UKB_BL_common$Head_Size)[,1]
      identical(UKB_BL_common$SubjectID, T2star$SubjectID)
      # run the model
      for (i in 1:nrow(T2star_results)) {
        UKB_BL_common$T2star_region <- scale(T2star[,i+1])[,1]
        
        out_temp <-  lm(formula = mdl_formula, data = UKB_BL_common)
        out <- summary(out_temp)
        T2star_results$t_value[which(T2star_results$Region==colnames(T2star)[i+1])] <- out$coefficients[2,3]
        T2star_results$p_value[which(T2star_results$Region==colnames(T2star)[i+1])] <- out$coefficients[2,4]
        T2star_results$r2[which(T2star_results$Region==colnames(T2star)[i+1])] <- out$r.squared
        T2star_results$effect_size[which(T2star_results$Region==colnames(T2star)[i+1])] <- out$coefficients[2,1]
        T2star_results$r2_adj[which(T2star_results$Region==colnames(T2star)[i+1])] <- out$adj.r.squared
        T2star_results$SE[which(T2star_results$Region==colnames(T2star)[i+1])] <- out$coefficients[2,2]
        
      }
      T2star_results$p_FDR <- p.adjust(T2star_results$p_value, method = "fdr", n = nrow(T2star_results))
      sum(T2star_results$p_FDR <= 0.05, na.rm = TRUE)
      T2star_results$isSignificant <- ifelse(T2star_results$p_FDR <= 0.05, 1, NA)
      T2star_results <- T2star_results %>% arrange(Region)
      # save and remove extra variables
      sample_summary[8, ] <- c('T2star', current_common_subjects_age_mean, current_common_subjects_age_sd, sum(UKB_BL_common$Sex == 0), sum(UKB_BL_common$Sex == 1), nrow(UKB_BL_common), (sum(UKB_BL_common$Sex == 0)/nrow(UKB_BL_common))*100, sum(T2star_results$p_FDR < 0.05, na.rm = TRUE), nrow(T2star_results))
      write.csv(T2star_results, paste0(current_output_folder, 'T2star_vs_', current_prs, '_Regression_Results.csv'), row.names = FALSE)
      rm(common_subjects, UKB_BL_common, current_common_subjects_age_mean, current_common_subjects_age_sd, T2star_results)
      
      
      ######## Save the sample summary ########
      # Add identifier columns and save it in the overall sample summary
      sample_summary$PRS <- current_prs
      sample_summary <- sample_summary %>% relocate(PRS, .before = 'Measure')
      write.csv(sample_summary, paste0(current_output_folder, 'sample_summary_', current_prs, '_Regression_Results.csv'), row.names = FALSE)
      
      sample_summary_all <- rbind(sample_summary_all, sample_summary)
      
      rm(variables)
    }
    
    write.csv(sample_summary_all, paste0('Outputs/PD_PRS_Regression_Results/', 'with_all_pathwayPRS_', model_motion, '_', model_sex, '/sample_summary_all_Regression_Results.csv'), row.names = FALSE)
    
    ## add sample_summary to the final large csv
    sample_summary_all <- sample_summary_all %>% 
      mutate(Cohort_sex = model_sex, Motion_confound_status = model_motion) %>% 
      relocate(Cohort_sex, .after = Measure) %>% relocate(Motion_confound_status, .after = Cohort_sex)
    sample_summary_all_merged <- rbind(sample_summary_all_merged, sample_summary_all)
    rm(sample_summary_all)
  }
}

## Save the final table
write.csv(sample_summary_all_merged, paste0('Outputs/PD_PRS_Regression_Results/', 'sample_summary_allPRS_allMeasures_allVersions_Regression_Results.csv'), row.names = FALSE)




















