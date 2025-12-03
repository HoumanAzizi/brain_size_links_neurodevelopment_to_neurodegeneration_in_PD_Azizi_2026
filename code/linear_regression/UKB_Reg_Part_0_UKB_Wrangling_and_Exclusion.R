library(tidyr)
library(dplyr)
library(readr)
library(data.table)
library(lme4)

rm(list = ls())
cat("\014")

setwd("/Users/houmanazizi/Library/Mobile\ Documents/com~apple~CloudDocs/Education/GitHub/UKB_Analysis/")

############################# Reading Data #############################
UKB_wide <- read_delim('Data/Current_UKB/2024March27/subset_2024March27_wide.tsv', delim = "\t", col_names = TRUE, na = "", quote = "")
UKB_wide_forColnames <- read.table('Data/Current_UKB/2024March27/subset_2024March27_wide.tsv', sep = "\t", header = TRUE, fill = TRUE)
colnames(UKB_wide) <- colnames(UKB_wide_forColnames)
UKB_BL <- UKB_wide %>% filter(InstanceID == 2)
# Fix family illnesses column by finding PD cases and making values 1 or 0 for each subject
UKB_BL <- UKB_BL %>% mutate(PD_History_1 = case_when(Illnesses.of.father_20107 == "Parkinson's disease" ~ 1, TRUE ~ 0))
UKB_BL <- UKB_BL %>% mutate(PD_History_2 = case_when(Illnesses.of.mother_20110 == "Parkinson's disease" ~ 1, TRUE ~ 0))
UKB_BL <- UKB_BL %>% mutate(PD_History_3 = case_when(Illnesses.of.siblings_20111 == "Parkinson's disease" ~ 1, TRUE ~ 0))
UKB_BL <- UKB_BL %>% mutate(PD_History_all = PD_History_1 + PD_History_2 + PD_History_3)
UKB_BL <- UKB_BL %>% group_by(SubjectID) %>% mutate(PD_Family_Hisotry = sum(PD_History_all))
UKB_BL <- UKB_BL %>% select(-c(Illnesses.of.father_20107,Illnesses.of.mother_20110,Illnesses.of.siblings_20111,
                               PD_History_1,PD_History_2,PD_History_3,PD_History_all))
UKB_BL <- as.data.frame(UKB_BL)
# Create a dummy variable for before or after educational reform
UKB_BL <- UKB_BL %>% mutate(Educational_Reform = case_when(Year.of.birth_34 > 1957 ~ 1,
                                                           Year.of.birth_34 == 1957 & Month.of.birth_52 == 'September' ~ 1,
                                                           Year.of.birth_34 == 1957 & Month.of.birth_52 == 'October' ~ 1,
                                                           Year.of.birth_34 == 1957 & Month.of.birth_52 == 'November' ~ 1,
                                                           Year.of.birth_34 == 1957 & Month.of.birth_52 == 'December' ~ 1,
                                                           TRUE ~ 0))



## NOTE: for PCA and ICD10 -> can just load the read.csv files ##
############################# Getting PCA's Per Subject #############################
# UKB_PCA <- UKB_BL %>% filter(!is.na(Genetic.principal.components_22009)) %>% select(c(SubjectID,ArrayID,Genetic.principal.components_22009))
# UKB_PCA <- UKB_PCA %>% mutate(PCA_Name = paste0("PCA_",ArrayID)) %>% select(-ArrayID)
# UKB_PCA <- UKB_PCA %>% pivot_wider(names_from = PCA_Name, values_from = Genetic.principal.components_22009)
# # Save PCAs separately as well
# write.csv(UKB_PCA, "Data/Genetic/PCA/UKB_PCAs.csv", row.names=FALSE)
UKB_PCA <- read.csv("Data/Genetic/PCA/UKB_PCAs.csv")




############################# Fixing Lang's PCA and PRS #############################
# PCA
PCA_Lang <- read.csv('Data/Genetic/PCA/ukbb_pc_unrelated_European_Lang.csv', fileEncoding = "UTF-7")
colnames(PCA_Lang) <- c("SubjectID", paste0("PCA_", 1:(ncol(PCA_Lang)-1)))
write.csv(PCA_Lang, "Data/Genetic/PCA/UKB_PCAs_Lang.csv", row.names=FALSE)
# PRSice
PRS_Lang <- read.csv('Data/Genetic/PD_PRS/PD_GWAS_2019_UKB_all_PRS_Lang.csv')
PRS_Lang <- PRS_Lang %>% filter(FID > 0) %>% select(-FID)
colnames(PRS_Lang) <- c('SubjectID', 'PD_PRScs')
write.csv(PRS_Lang, 'Data/Genetic/PD_PRS/PD_PRSice_Lang.csv', row.names = FALSE)
# PRScs
PRScs <- read.csv('Data/Genetic/PD_PRS/prscs_result_neurohub.csv')
PRScs <- PRScs %>% filter(NeuroHub_ID > 0) %>% select(-zscore)
colnames(PRScs) <- c('SubjectID', 'PD_PRScs')
write.csv(PRScs, 'Data/Genetic/PD_PRS/PD_PRScs_EricYu_Final.csv', row.names = FALSE)



############################# Fix PD Patient List to Use Later #############################
UKB_PD <- read_delim('Data/Current_UKB/PD_List/subset_PD_List_wide.tsv', delim = "\t", col_names = TRUE, na = "", quote = "")
UKB_PD <- as.data.frame(UKB_PD)
colnames(UKB_PD) <- c("SubjectID", "InstanceID", "ArrayID", "Date.of.attending.assessment.centre_53", "Date_PD_reported_131022")

UKB_PD <- UKB_PD %>% filter(!is.na(Date_PD_reported_131022)) %>% 
  filter(!is.na(Date.of.attending.assessment.centre_53)) %>% select(-ArrayID)
UKB_all_PD <- UKB_PD %>% select(-Date.of.attending.assessment.centre_53, -Date_PD_reported_131022) # this is to get a list of all PD subjects
UKB_PD <- UKB_PD %>% mutate(isPD = case_when(Date.of.attending.assessment.centre_53 <= Date_PD_reported_131022 ~ 1,
                                             TRUE ~ NA))
UKB_PD <- UKB_PD %>% select(-Date.of.attending.assessment.centre_53, -Date_PD_reported_131022)
UKB_PD <- UKB_PD %>% filter(InstanceID == 2 | InstanceID == 3) %>% filter(!is.na(isPD))
write.csv(UKB_PD, "Data/Current_UKB/PD_List/PD_List_GotPDAfterImaging.csv", row.names=FALSE)
write.csv(UKB_all_PD, "Data/Current_UKB/PD_List/PD_List_alltimes.csv", row.names=FALSE)







############################# Exclusion (Option 1 - Yashar): Removing ICD10_F or G > 1 #############################
# ## Getting ICD10 Per Subject
# # get only ICD10 columns and get related characters
# # Remove extra columns
# UKB_BL <- UKB_BL %>% filter(Age.when.attended.assessment.centre_21003 != '')
# UKB_BL <- as.data.frame(UKB_BL)
# rm(UKB_wide, UKB_wide_forColnames)
# # Actual exclusion
# UKB_only_ICD10 <- UKB_BL[,c("SubjectID","ArrayID","Diagnoses...ICD10_41270")]
# UKB_only_ICD10 <- UKB_only_ICD10 %>% filter(Diagnoses...ICD10_41270 != "")
# tmp <- UKB_only_ICD10 %>% mutate_all(function(x) substr(x,1,3)) # get the 3 related characters
# UKB_only_ICD10$Diagnoses...ICD10_41270 <- tmp$Diagnoses...ICD10_41270
# rm(tmp)
# ICDnames_full <- c("G00","G01","G02","G03","G04","G05","G06","G07","G08","G09","G10","G11","G12","G13","G14","G15","G16","G17","G18","G19",
#                    "G20","G21","G22","G23","G24","G25","G26","G27","G28","G29","G30","G31","G32","G33","G34","G35","G36","G37","G38","G39",
#                    "G40","G41","G42","G43","G44","G45","G46","G47","G48","G49","G50","G51","G52","G53","G54","G55","G56","G57","G58","G59",
#                    "G60","G61","G62","G63","G64","G65","G66","G67","G68","G69","G70","G71","G72","G73","G74","G75","G76","G77","G78","G79",
#                    "G80","G81","G82","G83","G84","G85","G86","G87","G88","G89","G90","G91","G92","G93","G94","G95","G96","G97","G98","G99",
#                    "F00","F01","F02","F03","F04","F05","F06","F07","F08","F09","F10","F11","F12","F13","F14","F15","F16","F17","F18","F19",
#                    "F20","F21","F22","F23","F24","F25","F26","F27","F28","F29","F30","F31","F32","F33","F34","F35","F36","F37","F38","F39",
#                    "F40","F41","F42","F43","F44","F45","F46","F47","F48","F49","F50","F51","F52","F53","F54","F55","F56","F57","F58","F59",
#                    "F60","F61","F62","F63","F64","F65","F66","F67","F68","F69","F70","F71","F72","F73","F74","F75","F76","F77","F78","F79",
#                    "F80","F81","F82","F83","F84","F85","F86","F87","F88","F89","F90","F91","F92","F93","F94","F95","F96","F97","F98","F99")
# # go over each eid and see if it has those variables
# UKB_only_ICD10 <- UKB_only_ICD10 %>% mutate(ICD10_GF = case_when(Diagnoses...ICD10_41270 %in% ICDnames_full ~ 1, TRUE ~ 0))
# UKB_only_ICD10 <- UKB_only_ICD10 %>% group_by(SubjectID) %>% mutate(total_icd10 = sum(ICD10_GF))
# UKB_only_ICD10 <- UKB_only_ICD10[duplicated(UKB_only_ICD10$SubjectID)==FALSE,c("SubjectID","total_icd10")]
# # Save ICD10 results separately
# write.csv(UKB_only_ICD10, "Data/Exclusions/UKB_ICD10_Counts.csv", row.names=FALSE)
# UKB_only_ICD10 <- read.csv("Data/Exclusions/UKB_ICD10_Counts.csv")
# # Apply the exclusion
# ids_to_keep_icd10 <- UKB_only_ICD10 %>% filter(total_icd10 == 0)
# UKB_BL <- UKB_BL %>% filter(SubjectID %in% ids_to_keep_icd10$SubjectID)
# rm(ICD10,ids_to_keep_icd10,UKB_only_ICD10,ICDnames_full)
# # Take only the main columsn from UKB_BL -> removing ICD10 and PCA data
# #UKB_BL <- UKB_BL %>% filter(Sex_31 != "") %>% select(-c(Diagnoses...ICD10_41270, Diagnoses...main.ICD10_41202, Genetic.principal.components_22009))
# UKB_BL <- UKB_BL %>% filter(Sex_31 != "")





############################# Exclusion (Option 2 - Nooshin): Removing Diseases using field 20002 #############################
## NOTE: read the csv OR use a file that has bipolar data to get the UKB_Disorders_Counts_fullUKB dataframe
list_of_exclusion_diseases <- c('parkinsons disease', 'dementia/alzheimers/cognitive impairment', 'chronic/degenerative neurological problem',
                                'acute infective polyneuritis/guillain-barre syndrome', 'multiple sclerosis',
                                'other demyelinating disease (not multiple sclerosis)', 'stroke', 'brain haemorrhage',
                                'brain abscess/intracranial abscess', 'cerebral aneurysm', 'cerebral palsy', 'encephalitis',
                                'epilepsy', 'head injury', 'infection of nervous system', 'ischaemic stroke',
                                'meningioma / benign meningeal tumour', 'meningitis', 'motor neurone disease', 'neurological injury/trauma',
                                'spina bifida', 'subdural haemorrhage/haematoma', 'subarachnoid haemorrhage', 'transient ischaemic attack (tia)')
UKB_wide <- UKB_wide %>% mutate(Exluded_Disease = case_when(Non.cancer.illness.code..self.reported_20002 %in% list_of_exclusion_diseases ~ 1, TRUE ~ 0))
UKB_wide <- UKB_wide %>% mutate(isBipolar = case_when(Bipolar.and.major.depression.status_20126 %in% c('Bipolar I Disorder', 'Bipolar II Disorder') ~ 1, TRUE ~ 0))
UKB_only_Disorders <- UKB_wide %>% group_by(SubjectID, InstanceID) %>% mutate(All_Neuro_Disorder_ID_Instance_Grouped = sum(Exluded_Disease, isBipolar))
UKB_only_Disorders <- as.data.frame(UKB_only_Disorders)
UKB_only_Disorders <- UKB_only_Disorders %>% group_by(SubjectID) %>% mutate(All_Neuro_Disorder_ID_Grouped = sum(Exluded_Disease, isBipolar)) %>% ungroup()
UKB_only_Disorders <- UKB_only_Disorders %>% select(SubjectID, InstanceID, All_Neuro_Disorder_ID_Grouped, All_Neuro_Disorder_ID_Instance_Grouped)
UKB_only_Disorders_Unique <- UKB_only_Disorders[!duplicated(UKB_only_Disorders[c("SubjectID", "InstanceID")]), ]
UKB_only_Disorders <- UKB_only_Disorders_Unique
# Save disorder results separately
#write.csv(UKB_only_Disorders_Unique, "Data/Exclusions/UKB_Disorders_Counts_fullUKB.csv", row.names=FALSE)
UKB_only_Disorders <- read.csv("Data/Exclusions/UKB_Disorders_Counts_fullUKB.csv")
# Apply the exclusion
ids_to_keep_disorders <- UKB_only_Disorders %>% filter(All_Neuro_Disorder_ID_Grouped == 0)
UKB_BL <- UKB_BL %>% filter(SubjectID %in% ids_to_keep_disorders$SubjectID)
rm(UKB_only_Disorders,ids_to_keep_disorders)


# Remove extra columns
UKB_BL <- UKB_BL %>% filter(Age.when.attended.assessment.centre_21003 != '')
UKB_BL <- as.data.frame(UKB_BL)








############################# Exclusion: BMI over 35 #############################
UKB_BL <- UKB_BL %>% filter(Body.mass.index..BMI._21001 <= 35)





############################# Exclusion: Family history of PD #############################
UKB_BL <- UKB_BL %>% filter(PD_Family_Hisotry == 0)






############################# Exclusion: Genetic and Self-Reported Sex Mismatch #############################
UKB_BL <- UKB_BL %>% filter(Sex_31 == Genetic.sex_22001)






############################# Exclusion: Removing non-European subjects #############################
# Using reported ethnicity - wrong method
# ethic_ids_to_keep <- UKB_wide %>% filter(!is.na(Ethnic.background_21000)) %>% select(SubjectID, Ethnic.background_21000)
# ethic_ids_to_keep <- ethic_ids_to_keep %>% filter(Ethnic.background_21000 == 'British')
# ethic_ids_to_keep <- ethic_ids_to_keep$SubjectID
# UKB_BL <- UKB_BL %>% filter(SubjectID %in% ethic_ids_to_keep)
# Using genetic ethnic grouping
UKB_BL <- UKB_BL %>% filter(Genetic.ethnic.grouping_22006 == "Caucasian")




############################# Exclusion (Option 1 - Yashar): Removing Related Subjects #############################
un_related <- read.csv("Data/Exclusions/ukbb_unrelated_European_covar.csv", sep=",", header = T)
UKB_BL <- UKB_BL %>% filter(SubjectID %in% un_related$ID)

############################# Exclusion (Option 2 - Nooshin): Removing Related Subjects #############################
un_related <- read.csv("Data/Exclusions/subjectsToExclude_kinship_Nooshin.csv", sep=",", header = T)
UKB_BL <- UKB_BL %>% filter(!(SubjectID %in% un_related$V1))

############################# Exclusion (Option 3 - Houman): Removing Related Subjects directly from UKB data #############################
# Find and remove one of the related subjects in the current remaining participants only
related_subjects <- UKB_wide %>% filter(InstanceID == 0 & ArrayID == 0) %>% 
  select(SubjectID, Genetic.relatedness.pairing_22011) %>% 
  filter(!is.na(Genetic.relatedness.pairing_22011))
related_subjects <- as.data.frame(related_subjects)
related_subjects <- related_subjects %>% filter(SubjectID %in% UKB_BL$SubjectID)
related_subjects <- related_subjects[duplicated(related_subjects$Genetic.relatedness.pairing_22011), 1]
UKB_BL <- UKB_BL %>% filter(!(SubjectID %in% related_subjects))
rm(UKB_wide, UKB_wide_forColnames)





# Get head size and add to the file
headsize <- read_delim('Data/Current_UKB/2023Nov23/subset_Nov23_NooshinPaper_Fig5_wide.tsv', delim = "\t", col_names = TRUE, na = "", quote = "")
headsize_forColnames <- read.table('Data/Current_UKB/2023Nov23/subset_Nov23_NooshinPaper_Fig5_wide.tsv', sep = "\t", header = TRUE, fill = TRUE)
colnames(headsize) <- colnames(headsize_forColnames)
headsize <- headsize %>% filter(InstanceID == 2, ArrayID == 0, !is.na(Volume.of.EstimatedTotalIntraCranial..whole.brain._26521)) %>% 
  select(SubjectID, Volume.of.EstimatedTotalIntraCranial..whole.brain._26521)
rm(headsize_forColnames)
UKB_BL <- UKB_BL %>% left_join(headsize, by = 'SubjectID')





############################# Add PCAs to end of file #############################
UKB_BL <- UKB_BL %>% left_join(UKB_PCA, by = "SubjectID")






############################# Separate UKB from Subcortical volume #############################
subcortical_vols <- UKB_BL %>% select(SubjectID, starts_with('Volume')) %>% select(-Volume.of.EstimatedTotalIntraCranial..whole.brain._26521)
subcortical_vols <- subcortical_vols %>% filter(!is.na(Volume.of.thalamus..left._25011))
UKB_BL <- UKB_BL %>% select(-ArrayID, -colnames(subcortical_vols)[2:ncol(subcortical_vols)])
# Fix subcortical colnames
subcortical_vols <- subcortical_vols %>% select(c("SubjectID",
                                                   "Volume.of.thalamus..left._25011", "Volume.of.thalamus..right._25012", "Volume.of.caudate..left._25013", "Volume.of.caudate..right._25014", 
                                                   "Volume.of.putamen..left._25015", "Volume.of.putamen..right._25016", "Volume.of.pallidum..left._25017", "Volume.of.pallidum..right._25018", 
                                                   "Volume.of.hippocampus..left._25019", "Volume.of.hippocampus..right._25020", "Volume.of.amygdala..left._25021", "Volume.of.amygdala..right._25022",
                                                   "Volume.of.accumbens..left._25023", "Volume.of.accumbens..right._25024"))
colnames(subcortical_vols) <- c("SubjectID",
                                 "Thalamus_L", "Thalamus_R", "Caudate_L", "Caudate_R", 
                                 "Putamen_L", "Putamen_R", "Pallidum_L", "Pallidum_R", 
                                 "Hippocampus_L", "Hippocampus_R", "Amygdala_L", "Amygdala_R",
                                 "Accumbens_L", "Accumbens_R")


############################# Save Final file for Analysis #############################
write.csv(UKB_BL, "Data/Current_UKB/UKB_Imaging1_withPCA_afterExclusion.csv", row.names=FALSE)






############################# Remove +-3SD from brain values and save #############################
SD_threshold <- 3

# Subcotical Volumes
average_zscore <- as.numeric(scale(rowSums(subcortical_vols[,-1], na.rm = TRUE)))
removed_subjects <- data.frame(SubjectID = subcortical_vols[abs(average_zscore) > SD_threshold, 1]) # keeping the list of subjects removed here
subcortical_vols <- subcortical_vols[average_zscore^2 <= SD_threshold^2,]
write.csv(subcortical_vols, "Data/Subcortical/UKB_Subcortical_Volumes_3SDremoved.csv", row.names=FALSE)
write.csv(removed_subjects, "Data/Subcortical/UKB_Subcortical_Subjects_Remove_by3SD.csv", row.names=FALSE)
rm(removed_subjects)

# Surface Area Schaefer100
SA_subjects <- read.csv('Data/Surface_Area/SA_in_Schaefer/SA_list.csv', header = TRUE)
SA100 <- read.csv('Data/Surface_Area/SA_vec_100_10mm.csv', header = FALSE)
SA100$SubjectID <- as.numeric(SA_subjects$eid_n)
SA100 <- SA100 %>% select(SubjectID, everything())
average_zscore <- as.numeric(scale(rowSums(SA100[,-1], na.rm = TRUE)))
removed_subjects <- data.frame(SubjectID = SA100[abs(average_zscore) > SD_threshold, 1])
SA100 <- SA100[average_zscore^2 <= SD_threshold^2,]
write.csv(SA100, "Data/Surface_Area/UKB_SA_Schaefer100_10mm_3SDremoved.csv", row.names=FALSE)
write.csv(removed_subjects, "Data/Surface_Area/UKB_SA_Schaefer100_10mm_Subjects_Remove_by3SD.csv", row.names=FALSE)
rm(removed_subjects)

# Surface Area Schaefer200
SA200 <- read.csv('Data/Surface_Area/SA_vec_200_10mm.csv', header = FALSE)
SA200$SubjectID <- as.numeric(SA_subjects$eid_n)
SA200 <- SA200 %>% select(SubjectID, everything())
average_zscore <- as.numeric(scale(rowSums(SA200[,-1], na.rm = TRUE)))
removed_subjects <- data.frame(SubjectID = SA200[abs(average_zscore) > SD_threshold, 1])
SA200 <- SA200[average_zscore^2 <= SD_threshold^2,]
write.csv(SA200, "Data/Surface_Area/UKB_SA_Schaefer200_10mm_3SDremoved.csv", row.names=FALSE)
write.csv(removed_subjects, "Data/Surface_Area/UKB_SA_Schaefer200_10mm_Subjects_Remove_by3SD.csv", row.names=FALSE)
rm(removed_subjects)

## Surface Area DKT
SA_subjects_DKT <- read.csv('Data/Surface_Area/SA_in_DKT/CIVET_in_DKT_IDs_UKB.csv', header = FALSE)
SA_DKT <- read.csv('Data/Surface_Area/SA_in_DKT/SA_CIVET_inDKT_UKB_10mm.csv', header = FALSE)
SA_DKT$SubjectID <- as.numeric(SA_subjects_DKT$V1)
SA_DKT <- SA_DKT %>% select(SubjectID, everything())
average_zscore <- as.numeric(scale(rowSums(SA_DKT[,-1], na.rm = TRUE)))
removed_subjects <- data.frame(SubjectID = SA_DKT[abs(average_zscore) > SD_threshold, 1])
SA_DKT <- SA_DKT[average_zscore^2 <= SD_threshold^2,]
# set column names
DKT_labels <- read.csv('Data/Atlases/CIVET_DKT_labels_withMiddleFilling.csv')
region_names <- paste0('Klein_', DKT_labels$Klein_original_number, '_CIVET_', DKT_labels$CIVET_DKT_number)
colnames(SA_DKT) <- c('SubjectID', region_names)
# save
write.csv(SA_DKT, "Data/Surface_Area/SA_in_DKT/UKB_SA_DKT_10mm_withMiddleFilling_3SDremoved.csv", row.names=FALSE)
write.csv(SA_DKT %>% select(-starts_with('Klein_NA_')), "Data/Surface_Area/SA_in_DKT/UKB_SA_DKT_10mm_3SDremoved.csv", row.names=FALSE)
write.csv(removed_subjects, "Data/Surface_Area/SA_in_DKT/UKB_SA_DKT_10mm_Subjects_Remove_by3SD.csv", row.names=FALSE)
rm(removed_subjects)

# Cortical Thickness Schaefer100
CT_subjects <- read.csv('Data/Cortical_Thickness/CT_list.csv', header = TRUE)
CT100 <- read.csv('Data/Cortical_Thickness/CT_vec_100_10mm.csv', header = FALSE)
CT100$SubjectID <- as.numeric(CT_subjects$eid_n)
CT100 <- CT100 %>% select(SubjectID, everything())
average_zscore <- as.numeric(scale(rowSums(CT100[,-1], na.rm = TRUE)))
removed_subjects <- data.frame(SubjectID = CT100[abs(average_zscore) > SD_threshold, 1])
CT100 <- CT100[average_zscore^2 <= SD_threshold^2,]
write.csv(CT100, "Data/Cortical_Thickness/UKB_CT_Schaefer100_10mm_3SDremoved.csv", row.names=FALSE)
write.csv(removed_subjects, "Data/Cortical_Thickness/UKB_CT_Schaefer100_10mm_Subjects_Remove_by3SD.csv", row.names=FALSE)
rm(removed_subjects)

## Cortical Thickness DKT
CT_subjects_DKT <- read.csv('Data/Cortical_Thickness/CT_in_DKT/CIVET_in_DKT_IDs_UKB.csv', header = FALSE)
CT_DKT <- read.csv('Data/Cortical_Thickness/CT_in_DKT/CT_CIVET_inDKT_UKB_10mm.csv', header = FALSE)
CT_DKT$SubjectID <- as.numeric(CT_subjects_DKT$V1)
CT_DKT <- CT_DKT %>% select(SubjectID, everything())
average_zscore <- as.numeric(scale(rowSums(CT_DKT[,-1], na.rm = TRUE)))
removed_subjects <- data.frame(SubjectID = CT_DKT[abs(average_zscore) > SD_threshold, 1])
CT_DKT <- CT_DKT[average_zscore^2 <= SD_threshold^2,]
# set column names
DKT_labels <- read.csv('Data/Atlases/CIVET_DKT_labels_withMiddleFilling.csv')
region_names <- paste0('Klein_', DKT_labels$Klein_original_number, '_CIVET_', DKT_labels$CIVET_DKT_number)
colnames(CT_DKT) <- c('SubjectID', region_names)
# save
write.csv(CT_DKT, "Data/Cortical_Thickness/CT_in_DKT/UKB_CT_DKT_10mm_withMiddleFilling_3SDremoved.csv", row.names=FALSE)
write.csv(CT_DKT %>% select(-starts_with('Klein_NA_')), "Data/Cortical_Thickness/CT_in_DKT/UKB_CT_DKT_10mm_3SDremoved.csv", row.names=FALSE)
write.csv(removed_subjects, "Data/Cortical_Thickness/CT_in_DKT/UKB_CT_DKT_10mm_Subjects_Remove_by3SD.csv", row.names=FALSE)
rm(removed_subjects)

# White Matter FA and MD
all_WM <- fread('Data/White_Matter/Tractwise/UKBB_WMA14Jun23.csv')
all_WM <- all_WM %>% group_by(eid) %>% mutate(dup = duplicated(Tract, fromLast = TRUE)) %>% filter(dup == FALSE)
all_WM <- all_WM %>% filter(Tract != '*')
FA_Tracts <- all_WM %>% select(c('eid', 'Tract', 'FA_MEAN')) %>% pivot_wider(names_from = Tract, values_from = FA_MEAN)
MD_Tracts <- all_WM %>% select(c('eid', 'Tract', 'MD_MEAN')) %>% pivot_wider(names_from = Tract, values_from = MD_MEAN)
rm(all_WM)
colnames(FA_Tracts)[1] <- 'SubjectID'
colnames(MD_Tracts)[1] <- 'SubjectID'
average_zscore <- as.numeric(scale(rowSums(FA_Tracts[,-1], na.rm = TRUE)))
removed_subjects_FA <- data.frame(SubjectID = FA_Tracts[abs(average_zscore) > SD_threshold, 1])
FA_Tracts <- FA_Tracts[average_zscore^2 <= SD_threshold^2,]
average_zscore <- as.numeric(scale(rowSums(MD_Tracts[,-1], na.rm = TRUE)))
removed_subjects_MD <- data.frame(SubjectID = MD_Tracts[abs(average_zscore) > SD_threshold, 1])
MD_Tracts <- MD_Tracts[average_zscore^2 <= SD_threshold^2,]
# fix tract names
tract_names <-read.csv('Data/White_Matter/ORG_tract_names_bilateral.csv', header = TRUE)
tract_names <- tract_names %>% arrange(match(Tract_file_name, colnames(FA_Tracts)[-1]))
colnames(FA_Tracts)[-1] <- tract_names$Tract_real_name
colnames(MD_Tracts)[-1] <- tract_names$Tract_real_name
identical(colnames(FA_Tracts),colnames(MD_Tracts))
write.csv(FA_Tracts, "Data/White_Matter/Tractwise/UKB_ORGTract_FA_3SDremoved.csv", row.names=FALSE)
write.csv(MD_Tracts, "Data/White_Matter/Tractwise/UKB_ORGTract_MD_3SDremoved.csv", row.names=FALSE)
write.csv(removed_subjects_FA, "Data/White_Matter/Tractwise/UKB_ORGTract_FA_Subjects_Remove_by3SD.csv", row.names=FALSE)
write.csv(removed_subjects_MD, "Data/White_Matter/Tractwise/UKB_ORGTract_MD_Subjects_Remove_by3SD.csv", row.names=FALSE)
rm(removed_subjects_FA, removed_subjects_MD)


# NODDI
NODDI <- read.csv('Data/NODDI/UKB_NODDI_ORG41_Tractwise_Data.csv', header = TRUE)
NODDI <- NODDI %>% filter(NODDI_Measure == 'ICVF') %>% select(-NODDI_Measure)
average_zscore <- as.numeric(scale(rowSums(NODDI[,-1], na.rm = TRUE)))
removed_subjects <- data.frame(SubjectID = NODDI[abs(average_zscore) > SD_threshold, 1])
NODDI <- NODDI[average_zscore^2 <= SD_threshold^2,]
write.csv(NODDI, "Data/NODDI/UKB_ORGTract_NODDI_ICVF_3SDremoved.csv", row.names=FALSE)
write.csv(removed_subjects, "Data/NODDI/UKB_ORGTract_NODDI_ICVF_Subjects_Remove_by3SD.csv", row.names=FALSE)
rm(removed_subjects)





############################# Fix SWI data, remove +-3SD, and save #############################
SWI <- read.csv('Data/SWI/ukb_category_109_bis.csv', header = TRUE)
# fix colnames
colnames(SWI)
column_mapping <- c("eid" = "SubjectID",
                    "X24483.2.0" = "Median_T2star_SubstantiaNigra_L",
                    "X24484.2.0" = "Median_T2star_SubstantiaNigra_R",
                    "X24479.2.0" = "Median_MagneticSusceptibility_Accumbens_L",
                    "X24480.2.0" = "Median_MagneticSusceptibility_Accumbens_R",
                    "X24477.2.0" = "Median_MagneticSusceptibility_Amygdala_L",
                    "X24478.2.0" = "Median_MagneticSusceptibility_Amygdala_R",
                    "X24469.2.0" = "Median_MagneticSusceptibility_Caudate_L",
                    "X24470.2.0" = "Median_MagneticSusceptibility_Caudate_R",
                    "X24475.2.0" = "Median_MagneticSusceptibility_Hippocampus_L",
                    "X24476.2.0" = "Median_MagneticSusceptibility_Hippocampus_R",
                    "X24473.2.0" = "Median_MagneticSusceptibility_Pallidum_L",
                    "X24474.2.0" = "Median_MagneticSusceptibility_Pallidum_R",
                    "X24471.2.0" = "Median_MagneticSusceptibility_Putamen_L",
                    "X24472.2.0" = "Median_MagneticSusceptibility_Putamen_R",
                    "X24481.2.0" = "Median_MagneticSusceptibility_SubstantiaNigra_L",
                    "X24482.2.0" = "Median_MagneticSusceptibility_SubstantiaNigra_R",
                    "X24467.2.0" = "Median_MagneticSusceptibility_Thalamus_L",
                    "X24468.2.0" = "Median_MagneticSusceptibility_Thalamus_R")
column_indices_onMap <- match(names(SWI), names(column_mapping))
column_indices_onSWI <- match(names(column_mapping), names(SWI))
colnames(SWI)[column_indices_onSWI] <- unname(column_mapping)
# select only the columns needed for first imaging
SWI <- SWI %>% select(-starts_with('X')) %>% filter(!is.na(Median_T2star_SubstantiaNigra_R))
# take the remaining T2start columns out
remaining_T2star <- SWI %>% select(SubjectID, starts_with("Median_T2star"))
SWI <- SWI %>% select(-starts_with("Median_T2star"))
# take out 3SD and save
average_zscore <- as.numeric(scale(rowSums(SWI[,-1], na.rm = TRUE)))
SWI <- SWI[average_zscore^2 <= SD_threshold^2,]
write.csv(SWI, "Data/SWI/UKB_Magnetic_Susceptibility_3SDremoved.csv", row.names=FALSE)


############################# Fix T2star data, remove +-3SD, and save #############################
T2star <- read.table('Data/SWI/subset_UKBB_SWI_IDPs_wide.tsv', sep = "\t", header = TRUE, fill = TRUE)
T2star <- T2star %>% filter(InstanceID == 2 & !is.na(Median.T2star.in.thalamus..left._25026)) %>% select(-ArrayID, -InstanceID, 
                                                                                                   -Susceptibility.weighted.brain.images...DICOM_20219, 
                                                                                                   -Susceptibility.weighted.brain.images...NIFTI_20251,
                                                                                                   -Discrepancy.between.SWI.brain.image.and.T1.brain.image_25738)
colnames(T2star)
colnames(T2star) <- c("SubjectID", 
                      "Median_T2star_Thalamus_L", "Median_T2star_Thalamus_R",
                      "Median_T2star_Caudate_L", "Median_T2star_Caudate_R", 
                      "Median_T2star_Putamen_L", "Median_T2star_Putamen_R", 
                      "Median_T2star_Pallidum_L", "Median_T2star_Pallidum_R", 
                      "Median_T2star_Hippocampus_L", "Median_T2star_Hippocampus_R",
                      "Median_T2star_Amygdala_L", "Median_T2star_Amygdala_R",
                      "Median_T2star_Accumbens_L", "Median_T2star_Accumbens_R")
# add remaining T2start columns
T2star <- T2star %>% left_join(remaining_T2star, by = 'SubjectID')
# remove 3SD and save
average_zscore <- as.numeric(scale(rowSums(T2star[,-1], na.rm = TRUE)))
T2star <- T2star[average_zscore^2 <= SD_threshold^2,]
write.csv(T2star, "Data/SWI/UKB_T2star_3SDremoved.csv", row.names=FALSE)




############################# Merge all PD-PRS (overall and pathway-specific PRScs and PRSice) #############################
### read all and rename columns
## patwhay PRS - all files together
# PRScs - Nalls GWAS
PRScs <- read.csv('Data/Genetic/PD_PRS/Pathway_PRS/full_pathway_PRScs-ice_withExclusions/UKB.PD.pathways.PRS.PRScs.csv')
colnames(PRScs)[1] <- 'SubjectID'
colnames(PRScs)[which(colnames(PRScs) == 'all')] <- 'PD_PRScs_Nalls'
colnames(PRScs)[2:(ncol(PRScs)-1)] <- paste0('PD_PRScs_', colnames(PRScs)[2:(ncol(PRScs)-1)])
# PRScs - GP2 2025 GWAS
PRScs_GP2 <- read.csv('Data/Genetic/PD_PRS/Pathway_PRS/full_pathway_PRScs-ice_withExclusions/UKB.PD.pathways.PRS.PRScs_GP2_2025.csv')
PRScs_GP2 <- PRScs_GP2 %>% select(-X.FID)
colnames(PRScs_GP2)[1] <- 'SubjectID'
colnames(PRScs_GP2)[which(colnames(PRScs_GP2) == 'all')] <- 'PD_PRScs_GP2'
colnames(PRScs_GP2)[2:(ncol(PRScs_GP2)-1)] <- paste0('PD_PRScs_GP2_', colnames(PRScs_GP2)[2:(ncol(PRScs_GP2)-1)])
# PRSice
PRSice <- read.csv('Data/Genetic/PD_PRS/Pathway_PRS/full_pathway_PRScs-ice_withExclusions/UKB.PD.pathways.PRS.PRSice.csv')
colnames(PRSice)[1] <- 'SubjectID'
colnames(PRSice)[2:ncol(PRSice)] <- paste0('PD_PRSice_', colnames(PRSice)[2:ncol(PRSice)])
# Overall PRSice
PRSice_Lang <- read.csv('Data/Genetic/PD_PRS/archived/Older_PRS/PD_PRSice_Lang.csv')
colnames(PRSice_Lang) <- c('SubjectID', 'PD_PRSice_Nalls')


### Merge all
all_PD_PRS <- PRScs %>% 
  left_join(PRSice_Lang, by = 'SubjectID') %>%
  left_join(PRScs_GP2, by = 'SubjectID') %>%
  left_join(PRSice, by = 'SubjectID')


## Remove extra variables
rm(PRScs, PRScs_GP2, PRSice, PRSice_Lang)

## Remove Autophagy Any genes because too much missing
all_PD_PRS <- all_PD_PRS %>% select(-PD_PRScs_Autophagy_any_genes_5)
# Remove incomplete rows
all_PD_PRS <- all_PD_PRS[complete.cases(all_PD_PRS),]

### z-score each PRS
z_score <- function(x) { (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) }
for (col in colnames(all_PD_PRS)[-1]) {
  if (is.numeric(all_PD_PRS[[col]])) {
    all_PD_PRS[[col]] <- z_score(all_PD_PRS[[col]])
  }
}

### Get correlation and plot
cor_matrix <- cor(all_PD_PRS[-1], use = "pairwise.complete.obs")
# Plot the correlation matrix
library(corrplot)
pdf("Data/Genetic/PD_PRS/PD_PRS_allVersions_Final_correlation.pdf", width = 20, height = 20)
corrplot(cor_matrix, method = "color", type = "upper", order = 'hclust',
         tl.col = "black", tl.srt = 45, # Text label color and rotation
         tl.cex = 0.7, # Reduce text size for labels
         addCoef.col = "black", # Add correlation coefficient values
         number.cex = 0.5, # Size of the correlation coefficient text
         col = colorRampPalette(c("royalblue4", "white", "red3"))(200), # Reverse color palette
         main = "Correlation Matrix of PD PRS") # Title of the plot
dev.off()

# save
write.csv(all_PD_PRS, 'Data/Genetic/PD_PRS/PD_PRS_allVersions_Final.csv', row.names = FALSE)

# save a smaller version with only important pathways
filtered_PD_PRS <- all_PD_PRS %>% select(SubjectID, PD_PRScs_GP2, PD_PRScs_Nalls, PD_PRSice_Nalls,
                                         contains("Autophagy_direct_genes_4"),
                                         contains("Lysosomal_direct_genes_2"),
                                         contains("Mitochondrial_genes_1"),
                                         contains("Mitochondrial_Lysosomal_Autophagy_direct_genes_9"))
write.csv(filtered_PD_PRS, 'Data/Genetic/PD_PRS/PD_PRS_allVersions_Final_mainPathways.csv', row.names = FALSE)



### Check Nalls Sig
# test <- read.table('Data/Genetic/PD_PRS/Pathway_PRS/pathway_PRSice/UKB.snp', sep = '\t', header = TRUE)
# nalls <- read.table('Data/Genetic/PD_PRS/Nalls_PD_SNPs/META5.1800.nochr.uniqID.summary.tab', sep = '\t', header = TRUE)
# nalls <- nalls %>% 
#   separate(SNP, into = c("CHR", "BP"), sep = ":", extra = "merge") %>% 
#   mutate(BP = sub("_.*", "", BP))  # Remove the ending part after the first underscore
# test$uniqueID <- paste0(test$CHR, '_', test$BP)
# nalls$uniqueID <- paste0(nalls$CHR, '_', nalls$BP)
# table(nalls$uniqueID %in% test$uniqueID)















############################# FreeSurfer Results - from Yashar #############################
### Keep QCed subjects only & related values
FS_QCed <- read.csv('Data/FreeSurfer_Results_Yashar/FreeSurfer_QC_Yashar.csv', header = TRUE)
FS_QCed <- FS_QCed %>% select(subject_id, tiv, euler_num)
colnames(FS_QCed) <- c('SubjectID', 'TIV', 'Euler_Number')
write.csv(FS_QCed, "Data/Current_UKB/UKB_TIV_EulerNumber_QCed_FreeSurfer_Yashar.csv", row.names=FALSE)


### Surface Area & Cortical Thickness FS DKT
FS_all <- read.csv('Data/FreeSurfer_Results_Yashar/FreeSurfer_SA_CT_Yashar.csv', header = TRUE)
FS_SA <- FS_all %>% select(subject_id, label, SurfArea)
colnames(FS_SA) <- c('SubjectID', 'Region', 'SurfArea')
FS_CT <- FS_all %>% select(subject_id, label, ThickAvg)
colnames(FS_CT) <- c('SubjectID', 'Region', 'ThickAvg')
rm(FS_all)
# Pivot wider
FS_SA <- FS_SA %>% pivot_wider(id_cols = SubjectID,
                               names_from = Region,
                               values_from = SurfArea)
FS_CT <- FS_CT %>% pivot_wider(id_cols = SubjectID,
                               names_from = Region,
                               values_from = ThickAvg)
# remove temporalpole_right as it does not exist in other places
FS_SA <- FS_SA %>% select(-temporalpole_right)
FS_CT <- FS_CT %>% select(-temporalpole_right)
# Keep only QC passed subjects
FS_SA <- FS_SA %>% filter(SubjectID %in% FS_QCed$SubjectID)
FS_CT <- FS_CT %>% filter(SubjectID %in% FS_QCed$SubjectID)
# Remove 3SD
average_zscore <- as.numeric(scale(rowSums(FS_SA[,-1], na.rm = TRUE)))
removed_subjects_SA <- data.frame(SubjectID = FS_SA[abs(average_zscore) > SD_threshold, 1])
FS_SA <- FS_SA[average_zscore^2 <= SD_threshold^2,]
average_zscore <- as.numeric(scale(rowSums(FS_CT[,-1], na.rm = TRUE)))
removed_subjects_CT <- data.frame(SubjectID = FS_CT[abs(average_zscore) > SD_threshold, 1])
FS_CT <- FS_CT[average_zscore^2 <= SD_threshold^2,]
# set column names
DKT_labels <- read.csv('Data/Atlases/CIVET_DKT_labels_withMiddleFilling.csv')
DKT_labels$Region_name_FS[which(DKT_labels$Region_name_FS == '')] <- NA
DKT_labels <- DKT_labels %>% filter(!is.na(Region_name_FS))
region_names <- paste0('Klein_', DKT_labels$Klein_original_number, '_CIVET_', DKT_labels$CIVET_DKT_number)
FS_SA <- FS_SA[, c("SubjectID", DKT_labels$Region_name_FS)]
FS_CT <- FS_CT[, c("SubjectID", DKT_labels$Region_name_FS)]
colnames(FS_SA)[-1] <- region_names
colnames(FS_CT)[-1] <- region_names
# save
write.csv(FS_SA, "Data/Surface_Area/SA_in_DKT_FreeSurfer_Yashar/UKB_SA_DKT_FreeSurfer_Yashar_QCed_3SDremoved.csv", row.names=FALSE)
write.csv(removed_subjects_SA, "Data/Surface_Area/SA_in_DKT_FreeSurfer_Yashar/UKB_SA_DKT_FreeSurfer_Yashar_Subjects_Remove_by3SD.csv", row.names=FALSE)
rm(removed_subjects_SA)
write.csv(FS_CT, "Data/Cortical_Thickness/CT_in_DKT_FreeSurfer_Yashar/UKB_CT_DKT_FreeSurfer_Yashar_QCed_3SDremoved.csv", row.names=FALSE)
write.csv(removed_subjects_CT, "Data/Cortical_Thickness/CT_in_DKT_FreeSurfer_Yashar/UKB_CT_DKT_FreeSurfer_Yashar_Subjects_Remove_by3SD.csv", row.names=FALSE)
rm(removed_subjects_CT)




### Subcortical Volume aseg FS DKT
FS_SV <- read.csv('Data/FreeSurfer_Results_Yashar/FreeSurfer_aseg_Yashar.csv', header = TRUE)
FS_SV <- FS_SV %>% filter(type == 'aseg') %>% select(subject_id, label, volume)
colnames(FS_SV) <- c('SubjectID', 'Region', 'volume')
# Pivot wider
FS_SV <- FS_SV %>% pivot_wider(id_cols = SubjectID,
                               names_from = Region,
                               values_from = volume)
# Select only the main regions
FS_SV <- FS_SV %>% select(SubjectID, `Left-Thalamus`, `Right-Thalamus`, `Left-Caudate`, `Right-Caudate`, 
                          `Left-Putamen`, `Right-Putamen`, `Left-Pallidum`, `Right-Pallidum`, 
                          `Left-Hippocampus`, `Right-Hippocampus`, `Left-Amygdala`, `Right-Amygdala`, 
                          `Left-Accumbens-area`, `Right-Accumbens-area`)
# Fix column names
colnames(FS_SV) <- c("SubjectID","Thalamus_L","Thalamus_R","Caudate_L","Caudate_R","Putamen_L","Putamen_R","Pallidum_L","Pallidum_R","Hippocampus_L","Hippocampus_R","Amygdala_L","Amygdala_R","Accumbens_L","Accumbens_R")
# Keep only QC passed subjects
FS_SV <- FS_SV %>% filter(SubjectID %in% FS_QCed$SubjectID)
# Remove 3SD
average_zscore <- as.numeric(scale(rowSums(FS_SV[,-1], na.rm = TRUE)))
removed_subjects <- data.frame(SubjectID = FS_SV[abs(average_zscore) > SD_threshold, 1])
FS_SV <- FS_SV[average_zscore^2 <= SD_threshold^2,]
# save
write.csv(FS_SV, "Data/Subcortical/SV_aseg_FreeSurfer_Yashar/UKB_Subcortical_Volumes_aseg_FreeSurfer_Yashar_QCed_3SDremoved.csv", row.names=FALSE)
write.csv(removed_subjects, "Data/Subcortical/SV_aseg_FreeSurfer_Yashar/UKB_Subcortical_aseg_FreeSurfer_Yashar_Subjects_Remove_by3SD.csv", row.names=FALSE)
rm(removed_subjects)




### TOTAL Subcortical Volume aseg FS DKT
FS_SV <- read.csv('Data/FreeSurfer_Results_Yashar/FreeSurfer_aseg_Yashar.csv', header = TRUE)
FS_SV <- FS_SV %>% filter(type == 'generals') %>% select(subject_id, label, volume)
colnames(FS_SV) <- c('SubjectID', 'Region', 'volume')
# Pivot wider
FS_SV <- FS_SV %>% pivot_wider(id_cols = SubjectID,
                               names_from = Region,
                               values_from = volume)
# Keep only QC passed subjects
FS_SV <- FS_SV %>% filter(SubjectID %in% FS_QCed$SubjectID)
# Remove 3SD
average_zscore <- as.numeric(scale(rowSums(FS_SV[,-1], na.rm = TRUE)))
removed_subjects <- data.frame(SubjectID = FS_SV[abs(average_zscore) > SD_threshold, 1])
FS_SV <- FS_SV[average_zscore^2 <= SD_threshold^2,]
# save
write.csv(FS_SV, "Data/Subcortical/SV_aseg_FreeSurfer_Yashar/UKB_aseg_Total_Values_FreeSurfer_Yashar_QCed_3SDremoved.csv", row.names=FALSE)
write.csv(removed_subjects, "Data/Subcortical/SV_aseg_FreeSurfer_Yashar/UKB_aseg_Total_Values_FreeSurfer_Yashar_Subjects_Remove_by3SD.csv", row.names=FALSE)
rm(removed_subjects)










############################# WM Results - re-extracted Houman #############################
all_WM <- read.csv('Data/White_Matter/Tractwise_reExtracted/UKB_WM_Regionwise_Data_WMMaskMultiplied.csv', header = TRUE)
FS_QCed <- read.csv('Data/Current_UKB/UKB_TIV_EulerNumber_QCed_FreeSurfer_Yashar.csv', header = TRUE)

# fix column names
colnames(all_WM) <- gsub("\\.nii\\.gz$", "", colnames(all_WM))
colnames(all_WM) <- gsub("\\.", "-", colnames(all_WM))
colnames(all_WM) <- gsub("IP", "I&P", colnames(all_WM))
colnames(all_WM)[-(1:2)] <- paste0(colnames(all_WM)[-(1:2)], ".nii")

# fix tract names
tract_names <- read.csv('Data/White_Matter/ORG_tract_names_bilateral.csv', header = TRUE)
setdiff(colnames(all_WM), tract_names$Tract_file_name)
tract_names <- tract_names %>% arrange(match(Tract_file_name, colnames(all_WM)[-(1:2)]))
colnames(all_WM)[-(1:2)] <- tract_names$Tract_real_name

# separate measures
FA_Tracts <- all_WM %>% filter(WM_Measure == 'FA') %>% select(-WM_Measure)
MD_Tracts <- all_WM %>% filter(WM_Measure == 'MD') %>% select(-WM_Measure)
FW_Tracts <- all_WM %>% filter(WM_Measure == 'FW') %>% select(-WM_Measure)
rm(all_WM)

colnames(FA_Tracts)[1] <- 'SubjectID'
colnames(MD_Tracts)[1] <- 'SubjectID'
colnames(FW_Tracts)[1] <- 'SubjectID'

# Keep only QC passed subjects
FA_Tracts <- FA_Tracts %>% filter(SubjectID %in% FS_QCed$SubjectID)
MD_Tracts <- MD_Tracts %>% filter(SubjectID %in% FS_QCed$SubjectID)
FW_Tracts <- FW_Tracts %>% filter(SubjectID %in% FS_QCed$SubjectID)

average_zscore <- as.numeric(scale(rowSums(FA_Tracts[,-1], na.rm = TRUE)))
removed_subjects_FA <- data.frame(SubjectID = FA_Tracts[abs(average_zscore) > SD_threshold, 1])
FA_Tracts <- FA_Tracts[average_zscore^2 <= SD_threshold^2,]
average_zscore <- as.numeric(scale(rowSums(MD_Tracts[,-1], na.rm = TRUE)))
removed_subjects_MD <- data.frame(SubjectID = MD_Tracts[abs(average_zscore) > SD_threshold, 1])
MD_Tracts <- MD_Tracts[average_zscore^2 <= SD_threshold^2,]
average_zscore <- as.numeric(scale(rowSums(FW_Tracts[,-1], na.rm = TRUE)))
removed_subjects_FW <- data.frame(SubjectID = FW_Tracts[abs(average_zscore) > SD_threshold, 1])
FW_Tracts <- FW_Tracts[average_zscore^2 <= SD_threshold^2,]

# save
write.csv(FA_Tracts, "Data/White_Matter/Tractwise_reExtracted/UKB_ORGTract_FA_reExtracted_QCed_3SDremoved.csv", row.names=FALSE)
write.csv(removed_subjects_FA, "Data/White_Matter/Tractwise_reExtracted/UKB_ORGTract_FA_reExtracted_Subjects_Remove_by3SD.csv", row.names=FALSE)
write.csv(MD_Tracts, "Data/White_Matter/Tractwise_reExtracted/UKB_ORGTract_MD_reExtracted_QCed_3SDremoved.csv", row.names=FALSE)
write.csv(removed_subjects_MD, "Data/White_Matter/Tractwise_reExtracted/UKB_ORGTract_MD_reExtracted_Subjects_Remove_by3SD.csv", row.names=FALSE)
write.csv(FW_Tracts, "Data/White_Matter/Tractwise_reExtracted/UKB_ORGTract_FW_reExtracted_QCed_3SDremoved.csv", row.names=FALSE)
write.csv(removed_subjects_FW, "Data/White_Matter/Tractwise_reExtracted/UKB_ORGTract_FW_reExtracted_Subjects_Remove_by3SD.csv", row.names=FALSE)
rm(removed_subjects_FA, removed_subjects_MD, removed_subjects_FW)








############################# WM Results - Global Measures by BISON WM Mask #############################
all_WM <- read.csv('Data/White_Matter/Global_Measures_BISON/global_measures_BISON_reExtracted/UKB_WM_Regionwise_Data_BISON.csv', header = TRUE)
FS_QCed <- read.csv('Data/Current_UKB/UKB_TIV_EulerNumber_QCed_FreeSurfer_Yashar.csv', header = TRUE)

# separate measures
FA_Global <- all_WM %>% filter(WM_Measure == 'FA') %>% select(-WM_Measure, -BISON_RF_icbm_asym_CSF_mean, -BISON_RF_icbm_asym_CSF_median, -BISON_RF_icbm_asym_GM_mean, -BISON_RF_icbm_asym_GM_median, -BISON_RF_icbm_asym_WM_median) %>% filter(!is.na(BISON_RF_icbm_asym_WM_mean))
MD_Global <- all_WM %>% filter(WM_Measure == 'MD') %>% select(-WM_Measure, -BISON_RF_icbm_asym_CSF_mean, -BISON_RF_icbm_asym_CSF_median, -BISON_RF_icbm_asym_GM_mean, -BISON_RF_icbm_asym_GM_median, -BISON_RF_icbm_asym_WM_median) %>% filter(!is.na(BISON_RF_icbm_asym_WM_mean))

colnames(FA_Global) <- c('SubjectID', 'Total_FA_BISON_WM')
colnames(MD_Global) <- c('SubjectID', 'Total_MD_BISON_WM')
rm(all_WM)

# Keep only QC passed subjects
FA_Global <- FA_Global %>% filter(SubjectID %in% FS_QCed$SubjectID)
MD_Global <- MD_Global %>% filter(SubjectID %in% FS_QCed$SubjectID)

average_zscore <- as.numeric(scale(FA_Global$Total_FA_BISON_WM))
removed_subjects_FA <- data.frame(SubjectID = FA_Global[abs(average_zscore) > SD_threshold, 1])
FA_Global <- FA_Global[average_zscore^2 <= SD_threshold^2,]
average_zscore <- as.numeric(scale(MD_Global$Total_MD_BISON_WM))
removed_subjects_MD <- data.frame(SubjectID = MD_Global[abs(average_zscore) > SD_threshold, 1])
MD_Global <- MD_Global[average_zscore^2 <= SD_threshold^2,]

# save
write.csv(FA_Global, "Data/White_Matter/Global_Measures_BISON/UKB_Global_WM_FA_BISON_QCed_3SDremoved.csv", row.names=FALSE)
write.csv(removed_subjects_FA, "Data/White_Matter/Global_Measures_BISON/UKB_Global_WM_FA_BISON_Subjects_Remove_by3SD.csv", row.names=FALSE)
write.csv(MD_Global, "Data/White_Matter/Global_Measures_BISON/UKB_Global_WM_MD_BISON_QCed_3SDremoved.csv", row.names=FALSE)
write.csv(removed_subjects_MD, "Data/White_Matter/Global_Measures_BISON/UKB_Global_WM_MD_BISON_Subjects_Remove_by3SD.csv", row.names=FALSE)
rm(removed_subjects_FA, removed_subjects_MD)
