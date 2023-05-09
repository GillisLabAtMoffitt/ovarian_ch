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
  ungroup()

write_rds(ovarian_dna, "ovarian_dna.rds")

# Cancer Characteristics----
subsequent_cancer <- ovarian_info %>%
  mutate(across(where(is.character), ~str_to_sentence(.))) %>% 
  filter(!str_detect(tumor_seq_num, "00")) %>%
  distinct(patient_id, dx_dt, .keep_all = TRUE) %>% 
  arrange(patient_id, dx_dt)

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

# Demographic----
Demographic <- Demographic %>%
  select(-c(is_moffitt_patient, is_active_tcc, is_approached_tcc_consent_oncore))


################################################################################# III ### Merge data
ovarian_patients <- ovarian_dna %>% 
  # Merge with Cancer Char for patients with samples available
  left_join(., ovarian_info, by = c("patient_id")) %>% 
  # Merge with Demographic for patients with samples available
  left_join(., Demographic, 
            c("patient_id"))
write_rds(ovarian_patients, "ovarian_patients.rds")

ovarian_patients_id <- paste0(ovarian_patients$patient_id, collapse = "|")


rm(Demographic, ovarian_DNA, ovarian_dna, ovarian_info)
################################################################################# IV ### Clean treatments
# Beva
Beva <- Chemot %>% 
  # Limit to ovarian cancer patients
  filter(str_detect(patient_id, ovarian_patients_id)) %>% 
  select(-c(treatment_start_dt_src, treatment_end_dt_src)) %>% 
  mutate(across(where(is.character), ~str_to_lower(.))) %>% 
  # limit to chemotherapy 
  filter(treatment_type == "immuno")


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
  arrange(patient_id, tumor_id, treatment_start_dt, treatment_drug_desc, treatment_type) %>% 
  # Combine drugs into regimen
  group_by(patient_id, tumor_id, treatment_start_dt, treatment_end_dt) %>%
  summarise_at(vars(treatment_drug_desc, treatment_type), 
               str_c, collapse = "; ") %>% 
  arrange(patient_id, tumor_id, treatment_start_dt, desc(treatment_end_dt), treatment_type) %>% 
  # Combine drugs into regimen
  group_by(patient_id, tumor_id, treatment_start_dt) %>%
  summarise_at(vars(treatment_drug_desc, treatment_end_dt, treatment_type), 
               str_c, collapse = "; ") %>% 
  ungroup() %>% 
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
  mutate(linenumber = row_number())


# pivot wider, use dcast bc better to keep date class
chemot <- dcast(setDT(Chemot), patient_id+tumor_id ~ rowid(patient_id),
                value.var = c(
                  "treatment_drug_desc",
                  "treatment_start_dt",
                  "treatment_end_dt1",
                  "treatment_type"
                )) %>% 
  select(patient_id, starts_with("treatment_drug_desc"),
         starts_with("chemotherapy_start"),
         starts_with("chemotherapy_end"),
         treatment_type_1, treatment_type_2) %>% 
  mutate(treatment_drug_desc_1 = case_when(
    !is.na(treatment_drug_desc_2)             
    ~ coalesce(treatment_drug_desc_1, paste(treatment_drug_desc_2, "for 2nd tumor")),
    TRUE                 ~ treatment_drug_desc_1
  ) ) %>% 
  mutate(treatment_drug_desc_2 = case_when(
    str_detect(treatment_drug_desc_1, "for 2nd tumor")       ~ NA_character_,
    TRUE                                                     ~ treatment_drug_desc_2
  ) ) %>% 
  mutate(treatment_drug_desc_2 = case_when(
    !is.na(treatment_drug_desc_3)             
    ~ coalesce(treatment_drug_desc_2, paste(treatment_drug_desc_3, "for 2nd tumor")),
    TRUE                 ~ treatment_drug_desc_2
  ) ) %>% 
  mutate(treatment_drug_desc_3 = case_when(
    str_detect(treatment_drug_desc_2, "for 2nd tumor")       ~ NA_character_,
    TRUE                                                     ~ treatment_drug_desc_3
  ) )

# Need to do the same for dates......might be better to use summarize then posixct










