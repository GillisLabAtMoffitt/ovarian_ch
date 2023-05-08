################################################################################# II ### Data cleaning

# Ovarian samples----
table(ovarian_DNA$sample_type)

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
ovarian_info <- ovarian_info %>%
  mutate(across(where(is.character), ~str_to_sentence(.))) %>% 
  filter(str_detect(tumor_seq_num, "00")) %>%
  # filter(str_detect(primary_site_group_desc, "OVARY")) %>%
  distinct(patient_id, dx_dt, .keep_all = TRUE) %>% 
  arrange(patient_id, dx_dt)












