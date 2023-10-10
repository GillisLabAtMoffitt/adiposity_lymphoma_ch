# Import Library
library(tidyverse)
library(lubridate)


################################################################################# I ### Load data
path <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", "adiposity_lymphoma_ch")

scans_data <- 
  readxl::read_xlsx(paste0(here::here(), "/10R23000188_20230816_outfile.xlsx"),
                    sheet = "PET_CT SCANS"
  ) %>% 
  janitor::clean_names() %>% 
  mutate(mrn = as.character(mrn))

scans_type <- 
  readxl::read_xlsx(paste0(here::here(), "/Lymphoma obesity_Imaging types_08.22.23dj.xlsx")
  ) %>% 
  janitor::clean_names()

dna_data <- 
  readxl::read_xlsx(paste0(here::here(), "/10R23000188_20230816_outfile.xlsx"),
                    sheet = "BioBanking Data"
  ) %>% 
  janitor::clean_names() %>% 
  mutate(mrn = as.character(mrn))

clinical <- 
  readxl::read_xlsx(paste0(here::here(), "/10R23000188_20230816_outfile.xlsx"),
                    sheet = "Cancer Registry"
  ) %>% 
  janitor::clean_names() %>% 
  mutate(mrn = as.character(mrn))

weight <- 
  readxl::read_xlsx(paste0(here::here(), "/10R23000188_20230816_outfile.xlsx"),
                    sheet = "HT_WT "
  ) %>% 
  janitor::clean_names() %>% 
  mutate(mrn = as.character(mrn))


################################################################################# II ### Data cleaning
# Scans
scans_type_contrast <- scans_type %>% 
  filter(ideal_cohort_c_ts_post_contrast_only == 1)
scans_type_petct <- scans_type %>% 
  filter(ideal_cohort_pet_ct_only == 1)

scans_contrast <- scans_data %>% 
  filter(charge_service_desc %in% c(scans_type_contrast$charge_service_desc)) %>% 
  select(mrn, contrast_scans_date = service_dt,  contrast_imaging_types = charge_service_desc)
  
scans_petct <- scans_data %>% 
  filter(charge_service_desc %in% c(scans_type_petct$charge_service_desc)) %>% 
  select(mrn, petct_scans_date = service_dt, petct_imaging_types = charge_service_desc)

rm(scans_data, scans_type,
   scans_type_contrast, scans_type_petct)


### DNA
dna_data <- dna_data %>% 
  group_by(mrn, tissue_type, sample_family_id, specimen_collection_dt) %>% 
  summarise_at(
    vars(sample_id, sample_type), 
    str_c, collapse = "; ") %>%
  group_by(mrn, tissue_type, specimen_collection_dt) %>% 
  summarise_at(
    vars(sample_family_id, sample_id, sample_type), 
    str_c, collapse = " or ") %>%
  ungroup()

germline_dna <- dna_data %>% 
  select(mrn, specimen_collection_dt, sample_type, 
         sample_family_id, sample_id, tissue_type) %>% 
  left_join(., clinical %>% 
              select(mrn, first_treatment_dt),
            by = "mrn") %>% 
  mutate(blood_bf_tx = case_when(
    specimen_collection_dt <= first_treatment_dt        ~ "Blood before",
    TRUE                                                ~ "No"
  )) %>% 
  mutate(interval_blood_tx = case_when(
    blood_bf_tx == "Blood before"                       ~ interval(
      start = specimen_collection_dt, end = first_treatment_dt) /
      duration(n=1, unit = "days")
  )) %>% 
  arrange(mrn, interval_blood_tx) %>% 
  distinct(mrn, .keep_all = TRUE) %>% 
  filter(!is.na(interval_blood_tx)) %>% 
  select(collection_dt_germline = specimen_collection_dt,
         sample_type_germline = sample_type, 
         sample_family_id_germline = sample_family_id,
         sample_id_germline = sample_id,
         -first_treatment_dt)

# Clinical
clinical <- clinical %>% 
  mutate(across(where(is.character), ~str_to_sentence(.))) %>% 
  filter(primary_site_region_desc == "Lymphoma" &
           histology_desc == "Malignant lymphoma large b-cell diffuse nos") %>% 
  mutate(surgery_radiation_seq_num = na_if(surgery_radiation_seq_num, "Not appl")) %>% 
  mutate(across(where(is.character), ~ na_if(., "Na")))

# Weight and height - select value before and closest to tx 





