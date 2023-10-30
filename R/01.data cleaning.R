# Import Library
library(tidyverse)
library(lubridate)


################################################################################# I ### Load data
path <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", "adiposity_lymphoma_ch")

scans_data <- 
  readxl::read_xlsx(paste0(here::here(), "/10R23000188_20231010_outfile_updated.xlsx"),
                    sheet = "PET_CT scans"
  ) %>% 
  janitor::clean_names() %>% 
  mutate(mrn = as.character(mrn))

scans_type <- 
  readxl::read_xlsx(paste0(here::here(), "/Lymphoma obesity_Imaging types_08.22.23dj.xlsx")
  ) %>% 
  janitor::clean_names()

dna_data <- 
  readxl::read_xlsx(paste0(here::here(), "/10R23000188_20231010_outfile_updated.xlsx"),
                    sheet = "BioBanking Data"
  ) %>% 
  janitor::clean_names() %>% 
  mutate(mrn = as.character(mrn))

clinical <- 
  readxl::read_xlsx(paste0(here::here(), "/10R23000188_20231010_outfile_updated.xlsx"),
                    sheet = "Cancer Registry"
  ) %>% 
  janitor::clean_names() %>% 
  mutate(mrn = as.character(mrn))

weight <- 
  readxl::read_xlsx(paste0(here::here(), "/10R23000188_20231010_outfile_updated.xlsx"),
                    sheet = "HT_WT "
  ) %>% 
  janitor::clean_names() %>% 
  mutate(mrn = as.character(mrn))


################################################################################# II ### Data cleaning
# Get body area being scanned  
scans_type_contrast <- scans_type %>% 
  filter(ideal_cohort_c_ts_post_contrast_only == 1)
scans_type_petct <- scans_type %>% 
  filter(ideal_cohort_pet_ct_only == 1)

# Scans # Make a contrast and non contrast data to find samples for each dataset
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
  select(mrn, collection_dt_germline = specimen_collection_dt,
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
height <- weight %>% 
  filter(!is.na(vitals_value_height)) %>% 
  mutate(height_date = as.Date(vitals_dtm)) %>% 
  select(mrn, treatment_start_dt, height = vitals_value_height, 
         height_unit = vitals_units_in, height_date) %>% 
  mutate(time_height_tx_days = case_when(
    height_date <= treatment_start_dt                   ~ interval(
      start = height_date, end = treatment_start_dt) /
      duration(n=1, unit = "days")
  )) %>% 
  mutate(time_height_after_tx_days = case_when(
    treatment_start_dt < height_date                   ~ interval(
      start = height_date, end = treatment_start_dt) /
      duration(n=1, unit = "days")
  )) %>% 
  filter(!is.na(height)) %>% 
  arrange(mrn, time_height_tx_days, time_height_after_tx_days) %>% 
  distinct(mrn, .keep_all = TRUE)

weight <- weight %>% 
  filter(!is.na(vitals_value_weight)) %>% 
  mutate(weight_date = as.Date(vitals_dtm)) %>% 
  select(mrn, treatment_start_dt, weight = vitals_value_weight, 
         weight_unit = vitals_units_lbs, weight_date) %>% 
  mutate(time_weight_tx_days = case_when(
    weight_date <= treatment_start_dt                   ~ interval(
      start = weight_date, end = treatment_start_dt) /
      duration(n=1, unit = "days")
  )) %>% 
  mutate(time_weight_after_tx_days = case_when(
    treatment_start_dt < weight_date                   ~ interval(
      start = weight_date, end = treatment_start_dt) /
      duration(n=1, unit = "days")
  )) %>% 
  # select(-treatment_start_dt) %>% 
  filter(!is.na(weight)) %>% 
  arrange(mrn, time_weight_tx_days, time_weight_after_tx_days) %>% 
  distinct(mrn, .keep_all = TRUE)

bmi <- full_join(weight, height,
                 by = c("mrn", "treatment_start_dt")) %>% 
  filter(!is.na(weight) & !is.na(height))


################################################################################# III ### Merge data
lymphoma_data <- clinical %>% 
  inner_join(., germline_dna,
             by = "mrn") %>% 
  inner_join(., bmi,
             by = "mrn")

write_rds(lymphoma_data, "lymphoma data.rds")

rm(clinical, dna_data, germline_dna,
   height, weight, bmi)


################################################################################# III ### Create variables
lymphoma_data <- read_rds(paste0(here::here(), "/lymphoma data.rds"))

lymphoma_data <- lymphoma_data %>% 
  mutate(weight_kg = weight / 2.205,
         height_m = height / 39.37,
         bmi = weight_kg / (height_m * height_m)) %>%
  mutate(bmi_cat = case_when(
    bmi < 25                    ~ "Underweight and normal weight",
    bmi >= 25 &
      bmi < 30                  ~ "Overweight",
    bmi >= 30                   ~ "Obese"
  )) %>%
  mutate(bmi_cat = factor(bmi_cat, levels = c("Underweight and normal weight", "Overweight", "Obese"))) #%>% 
  # mutate(os_time_months = interval(
  #   start = first_treatment_dt, end = date) /
  #     duration(n=1, unit = "months"))
  

write_rds(lymphoma_data, "lymphoma data with new variables.rds")


################################################################################# IV ### Find patients fitting criteria
# PetCT
lymphoma_data_petct <- lymphoma_data %>% 
  inner_join(., scans_petct,
             by = "mrn")

lymphoma_data_petct <- lymphoma_data_petct %>% 
  mutate(petct_bf_tx = case_when(
    petct_scans_date <= first_treatment_dt              ~ "Scan before",
    TRUE                                                ~ "No"
  )) %>% 
  mutate(time_petct_tx_days = case_when(
    petct_bf_tx == "Scan before"                        ~ interval(
      start = petct_scans_date, end = treatment_start_dt) /
      duration(n=1, unit = "days")
  )) %>% 
  # select(mrn,petct_scans_date, petct_bf_tx, treatment_start_dt, time_petct_tx_days) %>% 
  filter(!is.na(time_petct_tx_days)) %>% 
  arrange(mrn, time_petct_tx_days) %>% 
  distinct(mrn, .keep_all = TRUE)

table(lymphoma_data_petct$bmi_cat)
write_rds(lymphoma_data_petct, "lymphoma_data_petct.rds")

# Contrast
lymphoma_data_contrast <- lymphoma_data %>% 
  inner_join(., scans_contrast,
             by = "mrn")

lymphoma_data_contrast <- lymphoma_data_contrast %>% 
  mutate(contrast_bf_tx = case_when(
    contrast_scans_date <= first_treatment_dt           ~ "Scan before",
    TRUE                                                ~ "No"
  )) %>% 
  mutate(time_contrast_tx_days = case_when(
    contrast_bf_tx == "Scan before"                        ~ interval(
      start = contrast_scans_date, end = treatment_start_dt) /
      duration(n=1, unit = "days")
  )) %>% 
  # select(mrn,contrast_scans_date, contrast_bf_tx, treatment_start_dt, time_contrast_tx_days) %>% 
  filter(!is.na(time_contrast_tx_days)) %>% 
  arrange(mrn, time_contrast_tx_days) %>% 
  distinct(mrn, .keep_all = TRUE)

write_rds(lymphoma_data_contrast, "lymphoma_data_contrast.rds")
table(lymphoma_data_contrast$bmi_cat)


# End cleaning
