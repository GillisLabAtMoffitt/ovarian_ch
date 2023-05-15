# Import Library
library(tidyverse)
library(data.table)
library(lubridate)

################################################################################# II ### Data cleaning

# Ovarian samples----
ovarian_dna <- ovarian_DNA %>% 
  # Select sample type needed for CH detection
  filter(collection_site_anatomic_desc == "Blood", 
         str_detect(sample_type, "Buffy|Genomic|Unprocessed|CD138|MNC$")) %>% 
  # mutate(across(where(is.character), ~str_to_sentence(.))) %>% 
  select(patient_id, sample_family_id, sample_id,
         specimen_collection_dt) %>%
  # add same sample/same date on the same row
  arrange(patient_id, specimen_collection_dt) %>% 
  # Summarize to have 1 sample/day per row and not 1 row for each aliquot of the same sample 
  group_by(patient_id, sample_family_id, specimen_collection_dt) %>% 
  summarise_at(vars(sample_id), str_c, collapse = "; ") %>%
  # separate(col = sample_id, paste("sample_id", 1:3, sep="_"), sep = "; ", extra = "drop", fill = "right")
  ungroup() %>% 
  # Create deidentify IDs
  mutate(txt = "ovarian_study_") %>%
  group_by(patient_id) %>% 
  mutate(id = cur_group_id()) %>%
  ungroup() %>%
  mutate(zero = 6 - nchar(id)) %>%
  mutate(ii = stringi::stri_dup("0", zero)) %>%
  select(c(txt, ii, id, patient_id, everything())) %>%
  unite(deidentified_patient_id, txt:id, sep = "") %>% 
  select(c(deidentified_patient_id, patient_id, everything(), -zero))

write_rds(ovarian_dna, "ovarian_dna.rds")

# Cancer Characteristics----
subsequent_cancer <- ovarian_info %>%
  mutate(across(where(is.character), ~str_to_sentence(.))) %>% 
  filter(!str_detect(tumor_seq_num, "00")) %>%
  distinct(patient_id, dx_dt, .keep_all = TRUE) %>% 
  arrange(patient_id, dx_dt) %>% 
  select(patient_id, subsequent_tumor_id = tumor_id, 
         subsequent_tumor_seq_num = tumor_seq_num, 
         subsequent_tumor_dx_dt = dx_dt, 
         subsequent_tumor_primary_site_group = primary_site_group_desc, 
         subsequent_tumor_histology = histology_desc, 
         subsequent_tumor_stage_tnm = stage_tnm_cs_mixed_group_desc) %>% 
  mutate(had_subsequent_tumor = "had subsequent tumor")

ovarian_info <- ovarian_info %>%
  mutate(across(where(is.character), ~str_to_sentence(.))) %>% 
  filter(str_detect(tumor_seq_num, "00")) %>%
  # filter(str_detect(primary_site_group_desc, "OVARY")) %>%
  # Remove potential duplicated record
  distinct(patient_id, dx_dt, .keep_all = TRUE) %>% 
  arrange(patient_id, dx_dt) %>% 
  select(patient_id, tumor_id,
         last_contact_or_death_dt, dx_dt,
         tumor_behavior_desc, class_of_case_cd, grade_differentiation_cd,
         grade_differentiation_desc, grade_clinical_desc, grade_pathological_desc,
         grade_post_therapy_path_desc, comorbidity_01_desc,
         lymph_vascular_invasion_desc, regional_nodes_examined, 
         regional_nodes_positive, stage_clinical_tnm_group_desc, 
         stage_pathological_tnm_group_desc, stage_tnm_cs_mixed_group_desc,
         cancer_status_last_dt, type_of_first_recurrence_desc,
         metastatic_site_at_diagnosis_desc)

ovarian_info <- left_join(ovarian_info, 
                          subsequent_cancer,
                          by = "patient_id")

# Demographic----
Demographic <- Demographic %>%
  select(-c(is_moffitt_patient, is_active_tcc, is_approached_tcc_consent_oncore))


################################################################################# III ### Merge clinical data
ovarian_patients <- ovarian_dna %>% 
  # Merge with Cancer Char for patients with samples available
  left_join(., ovarian_info, by = c("patient_id")) %>% 
  # Merge with Demographic for patients with samples available
  left_join(., Demographic, 
            c("patient_id"))
write_rds(ovarian_patients, "ovarian_patients.rds")

ovarian_patients_id <- paste0(ovarian_patients$patient_id, collapse = "|")


rm(Demographic, ovarian_DNA, ovarian_dna, ovarian_info, 
   subsequent_cancer)
################################################################################# IV ### Clean treatments
# Beva
Beva <- Chemot %>% 
  # Limit to ovarian cancer patients
  filter(str_detect(patient_id, ovarian_patients_id)) %>% 
  select(-c(treatment_start_dt_src, treatment_end_dt_src)) %>% 
  mutate(across(where(is.character), ~str_to_lower(.))) %>% 
  # limit to chemotherapy 
  filter(treatment_type == "immuno") %>% 
  filter(str_detect(treatment_drug_desc, "beva")) %>% 
  select(patient_id, tumor_id, 
         bevacizumab_start_dt = treatment_start_dt, 
         bevacizumab_end_dt = treatment_end_dt) %>% 
  arrange(patient_id, bevacizumab_start_dt) %>% 
  # Combine drugs into regimen with same start and stop date
  group_by(patient_id, tumor_id) %>%
  summarise_at(vars(bevacizumab_start_dt, bevacizumab_end_dt), 
               str_c, collapse = "_") %>% 
  ungroup() %>% 
  mutate(regimen_count = sapply(strsplit(bevacizumab_start_dt, "_"), length)) %>% 
  separate_wider_delim(cols = bevacizumab_start_dt, delim = "_",
                       names = c(paste("bevacizumab_start_dt_regimen", 1:max(.$regimen_count), sep = "_")), 
                       too_few = "align_start", too_many = "drop", 
                       cols_remove = TRUE) %>% 
  separate_wider_delim(cols = bevacizumab_end_dt, delim = "_",
                       names = c(paste("bevacizumab_end_dt_regimen", 1:max(.$regimen_count), sep = "_")), 
                       too_few = "align_start", too_many = "drop", 
                       cols_remove = TRUE) %>% 
  select(-regimen_count) %>% 
  mutate(across(c(starts_with("bevacizumab_start_dt"), 
                  starts_with("bevacizumab_end_dt")), 
                ~ as.POSIXct(.))) %>% 
  mutate(received_beva = "Received Bevacizumab")


# Chemotherapy----
Chemot1 <- Chemot %>% 
  # Limit to ovarian cancer patients
  filter(str_detect(patient_id, ovarian_patients_id)) %>% 
  select(-c(treatment_start_dt_src, treatment_end_dt_src)) %>% 
  mutate(across(where(is.character), ~str_to_lower(.))) %>% 
  # limit to chemotherapy 
  filter(treatment_type %in% c("chemo", "other", "palliative")) %>% 
  
  # deeper cleaning
  # Fix 1 palliative data
  mutate(treatment_drug_desc = case_when(
    treatment_drug_desc == 
      "unknown; not stated"                  ~ NA_character_,
    treatment_type == "palliative" & 
      !is.na(treatment_start_dt) &
      str_detect(treatment_txt, "carbo")     ~ "carboplatin, paclitaxel",
    treatment_drug_desc == 
      "abraxane (form of taxol)"             ~ "paclitaxel",
    TRUE                                     ~ treatment_drug_desc
  )) %>% 
  # Select administered chemo
  filter(
    # condition chemo
    (treatment_type == "chemo" &
       !is.na(treatment_start_dt)
    ) |
      (treatment_type == "chemo" &
         is.na(treatment_start_dt) &
         treatment_rx_summary_desc == "multi-agent chemo"
      ) |
      # condition pal
      (treatment_type == "palliative" & 
         !is.na(treatment_start_dt) &
         !is.na(treatment_drug_desc)
      ) |
      # condition other
      (treatment_type == "other" &
         !is.na(treatment_start_dt)
      ) 
  ) %>% 
  # Create a really early date to add to NAs date to make sure to exclude these patients
  mutate(treatment_start_dt =
           coalesce(treatment_start_dt,
                    as.POSIXct("1700-01-01", tz = "GMT"))
  ) %>% 
  distinct(patient_id, treatment_drug_desc, treatment_start_dt, 
           treatment_end_dt, .keep_all = TRUE) 

Chemot <- Chemot1 %>%
  arrange(patient_id, treatment_start_dt, treatment_drug_desc, treatment_type) %>% 
  # Combine drugs into regimen with same start and stop date
  group_by(patient_id, tumor_id, treatment_start_dt, treatment_end_dt) %>%
  summarise_at(vars(treatment_drug_desc, treatment_type), 
               str_c, collapse = "; ") %>% 
  arrange(patient_id, treatment_start_dt, desc(treatment_end_dt), treatment_type) %>% 
  # Combine drugs into regimen with same start (still in the same regimen)
  # Do it in 2 step to keep an order in what was given together
  group_by(patient_id, tumor_id, treatment_start_dt) %>%
  summarise_at(vars(treatment_drug_desc, treatment_end_dt, treatment_type), 
               str_c, collapse = "; ") %>% 
  ungroup() %>% 
  mutate(treatment_drug_desc = case_when(
    is.na(treatment_drug_desc) &
      !is.na(treatment_start_dt)           ~ "Ukn what was given",
    TRUE                                   ~ treatment_drug_desc
  )) %>% 
  mutate(treatment_type = case_when(
    treatment_type == "chemo; chemo; chemo; chemo" |
      treatment_type == "chemo; chemo; chemo" |
      treatment_type == "chemo; chemo"               ~ "chemo",
    TRUE                                             ~ treatment_type
  )) %>% 
  mutate(enddate_count = sapply(strsplit(treatment_end_dt, ";"), length), .before = 5) %>% 
  separate_wider_delim(cols = treatment_end_dt, delim = "; ",
                       names = c(paste("treatment_end_dt", 1:max(.$enddate_count), sep = "")), 
                       too_few = "align_start", too_many = "drop", 
                       cols_remove = TRUE) %>% 
  mutate(across(starts_with("treatment_end_dt"), ~ as.POSIXct(.))) %>% 
  distinct(patient_id, treatment_drug_desc, treatment_start_dt, treatment_end_dt1, 
           .keep_all = TRUE) %>%
  group_by(patient_id) %>% 
  mutate(tumor_seq = dense_rank(tumor_id), .before = 3) %>% 
  group_by(patient_id, tumor_id) %>% 
  mutate(linenumber = row_number(), .before = 4)

write_rds(Chemot, "Chemotherapy.rds")

# pivot wider regimen but keep tumor on different rows
chemot <- Chemot %>% 
  arrange(patient_id, treatment_start_dt, treatment_end_dt1, treatment_type) %>% 
  group_by(patient_id, tumor_id, tumor_seq) %>% 
  summarise_at(vars(treatment_drug_desc, treatment_start_dt, treatment_end_dt1, treatment_type), 
               str_c, collapse = "_") %>% 
  ungroup() %>% 
  mutate(regimen_count = sapply(strsplit(treatment_start_dt, "_"), length), .before = 5) %>% 
  separate_wider_delim(cols = treatment_drug_desc, delim = "_",
                       names = c(paste("treatment_drug_regimen", 1:max(.$regimen_count), sep = "_")), 
                       too_few = "align_start", too_many = "drop", 
                       cols_remove = TRUE) %>% 
  separate_wider_delim(cols = treatment_start_dt, delim = "_",
                       names = c(paste("treatment_start_dt_regimen", 1:max(.$regimen_count), sep = "_")), 
                       too_few = "align_start", too_many = "drop", 
                       cols_remove = TRUE) %>% 
  separate_wider_delim(cols = treatment_end_dt1, delim = "_",
                       names = c(paste("treatment_end_dt1_regimen", 1:max(.$regimen_count), sep = "_")), 
                       too_few = "align_start", too_many = "drop", 
                       cols_remove = TRUE) %>% 
  mutate(across(c(starts_with("treatment_start_dt"), 
                  starts_with("treatment_end_dt")), 
                ~ as.POSIXct(.))) %>% 
  mutate(subsequent_tumor_treatment_tumor_id = case_when(
    tumor_seq == 2                    ~ tumor_id
  )) %>% 
  mutate(subsequent_tumor_treatment_start_dt = case_when(
    tumor_seq == 2                    ~ treatment_start_dt_regimen_1
  )) %>% 
  mutate(subsequent_tumor_treatment_drug = case_when(
    tumor_seq == 2                    ~ treatment_drug_regimen_1
  )) %>% 
  group_by(patient_id) %>% 
  fill(subsequent_tumor_treatment_tumor_id, subsequent_tumor_treatment_start_dt, 
       subsequent_tumor_treatment_drug, .direction = "updown") %>% 
  distinct(patient_id, .keep_all = TRUE)

write_rds(chemot, "chemotherapy_wide.rds")


################################################################################# IV ### Merge all
ovarian_data <- 
  left_join(ovarian_patients, chemot, by = c("patient_id", "tumor_id")) %>% 
  left_join(., Beva, by = c("patient_id", "tumor_id"))

write_rds(ovarian_data, "ovarian_data_with_blood.rds")
rm(Beva, Chemot1, ovarian_patients_id)

# END cleaning




