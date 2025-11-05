# Import Library
library(tidyverse)
library(lubridate)


################################################################################# I ### Load data----
path <- fs::path("", "Volumes", "Gillis_Research",
                 "Lab_Data", "LymphomaObesity")

clinical <- 
  readxl::read_xlsx(paste0(
    # path, "/RawData",
    here::here(), "/data/raw data",
    "/Lymphoma_ClinicalData_10R24000003_20240131.xlsx"),
    sheet = "Cancer Registry"
  ) %>% 
  janitor::clean_names() %>% 
  mutate(mrn = as.character(mrn))

vitals <- 
  readxl::read_xlsx(paste0(
    # path, "/RawData",
    here::here(), "/data/raw data",
    "/Lymphoma_ClinicalData_10R24000003_20240131.xlsx"),
    sheet = "Vital Status"
  ) %>% 
  janitor::clean_names() %>% 
  mutate(mrn = as.character(mrn)) %>% 
  select(mrn, date_of_last_contact, vital_status,
         starts_with("cause_of_death_"))

samples_selected <- 
  read_csv(paste0(# path,
    here::here(), "/data/processed data",
    "/lymphoma data with new variables_01082024.csv")) %>% 
  select(mrn, tumor_id,
         collection_dt_germline, sample_type_germline,
         sample_family_id_germline, 
         sample_id_requested = sample_id_germline)

weight <- 
  readxl::read_xlsx(paste0(
    # path, "/RawData",
    here::here(), "/data/raw data",
    "/Lymphoma_ClinicalData_10R24000003_20240131.xlsx"),
    sheet = "HT_WT"
  ) %>% 
  janitor::clean_names() %>% 
  mutate(mrn = as.character(mrn))

body_comp <- 
  readxl::read_xlsx(paste0(# path,
    here::here(), "/data/raw data",
    "/Lymphoma_SliceomaticResults_20240216.xlsx"), n_max = 105) %>% 
  janitor::clean_names()


ch_calls <- 
  readxl::read_xlsx(paste0(# path,
    here::here(), "/data/raw data",
    "/Lymphoma_Shipping Manifest_EXP2052DNA_PHI_20240504.xlsx"), 
    sheet = "Shipping Manifest_w PHI") %>% 
  janitor::clean_names()
ch_mutation <- 
  readxl::read_xlsx(paste0(# path,
    here::here(), "/data/raw data",
    "/Lymphoma_CHcalls_20251021.xlsx"), 
    na = "NA") %>% 
  janitor::clean_names()


################################################################################# II ### Data cleaning----
samples_selected <- samples_selected %>% 
  mutate(sample_id_requested = str_replace(sample_id_requested, " or ", "; ")) %>% 
  separate_wider_delim(cols = sample_id_requested, delim = "; ",
                       names = c("sample_id_requested_1", "sample_id_requested_2"), 
                       too_few = "align_start", too_many = "merge", 
                       cols_remove = TRUE) %>% 
  pivot_longer(cols = c(sample_id_requested_1, sample_id_requested_2), 
               names_to = NULL, 
               values_to = "sample_id_requested", 
               values_drop_na = TRUE)

ch_mutation <- ch_mutation %>% 
  mutate(ch_status = factor(ch_status, levels = c("No CH", "CH")))
ch_calls <- ch_calls %>% 
  select(mrn : ncol(ch_calls)) %>% 
  full_join(ch_mutation %>% 
              select(source_sample_id, ch_status) %>% 
              distinct(), .,
            by = "source_sample_id") %>% 
  filter(str_detect(study_id,  "Lymphoma_")) %>% 
  select(mrn, study_id, ch_status, 
         sample_id_requested = source_sample_id, 
         sample_id_germline = dna_aliquot_number_1_lv_id) %>% 
  left_join(., samples_selected,
            by = c("mrn", "sample_id_requested")) %>% 
  mutate(mrn = as.character(mrn))

body_comp <- body_comp %>% 
  # Extract mrn and scan date
  mutate(mrn = str_extract(slice_name, "([:digit:]*)"), 
         .before = 1) %>% 
  mutate(scan_date = str_match(slice_name, 
                               "[:digit:]*\\.([:digit:]*\\.[:digit:]*\\.[:digit:]*).*$")[,2], 
         .after = mrn) %>% 
  mutate(scan_date = as.POSIXct(scan_date, format = "%m.%d.%Y")) %>% 
  # This is to remove the duplicated id for which we have 2 identical measurement at 2 different dates.
  # I am using the latest date which is related to the dx date.
  arrange(mrn, desc(scan_date)) %>% 
  distinct(mrn, .keep_all = TRUE)

# Weight and height - select value before and closest to tx 
# I am re-doing this step with the new updated data received
weight <- weight %>% 
  filter(str_detect(mrn, paste0(ch_calls$mrn, collapse = "|"))) %>% 
  left_join(., body_comp %>% 
            select(mrn, scan_date),
          by = "mrn") %>% 
  mutate(vitals_dtm = as.Date(vitals_dtm))

height <- weight %>% 
  filter(!is.na(vitals_value_height)) %>% 
  select(mrn, scan_date,
         height = vitals_value_height, 
         height_unit = vitals_units_in, 
         height_date = vitals_dtm) %>% 
  mutate(time_height_to_scan_days = case_when(
    height_date <= scan_date                   ~ interval(
      start = height_date, end = scan_date) /
      duration(n=1, unit = "days")
  )) %>% 
  mutate(time_height_after_scan_days = case_when(
    scan_date < height_date                   ~ interval(
      start = scan_date, end = height_date) /
      duration(n=1, unit = "days")
  )) %>% 
  filter(!is.na(height)) %>% 
  arrange(mrn, time_height_to_scan_days, time_height_after_scan_days) %>% 
  distinct(mrn, .keep_all = TRUE) %>% 
  select(-scan_date) %>% 
  mutate_at(c("time_height_to_scan_days",
              "time_height_after_scan_days"), ~ round(., 0))

weight <- weight %>% 
  filter(!is.na(vitals_value_weight)) %>% 
  rename(weight_date = vitals_dtm) %>% 
  select(mrn, treatment_start_dt, scan_date,
         weight = vitals_value_weight, 
         weight_unit = vitals_units_lbs, weight_date) %>% 
  mutate(time_weight_to_scan_days = case_when(
    weight_date <= scan_date                   ~ interval(
      start = weight_date, end = scan_date) /
      duration(n=1, unit = "days")
  )) %>% 
  mutate(time_weight_after_scan_days = case_when(
    scan_date < weight_date                   ~ interval(
      start = scan_date, end = weight_date) /
      duration(n=1, unit = "days")
  )) %>% 
  # select(-treatment_start_dt) %>% 
  filter(!is.na(weight)) %>% 
  arrange(mrn, time_weight_to_scan_days, time_weight_after_scan_days) %>% 
  distinct(mrn, .keep_all = TRUE) %>% 
  select(-scan_date) %>% 
  mutate_at(c("time_weight_to_scan_days",
              "time_weight_after_scan_days"), ~ round(., 0))


################################################################################# III ### Merge and create new varaibles----
bmi <- full_join(weight, height,
                 by = c("mrn")) %>% 
  filter(!is.na(weight) & !is.na(height))

lymphoma_data <- ch_calls %>% 
  left_join(., clinical,
            by = c("mrn", "tumor_id")) %>% 
  left_join(., vitals,
            by = "mrn") %>% 
  left_join(., body_comp,
            by = c("mrn"))%>% 
  left_join(., bmi,
            by = "mrn")

lymphoma_data <- lymphoma_data %>% 
  # General cleaning
  mutate_at(c("gender_src_desc", 
              "race_cr_src_desc_1"), 
            ~ str_to_sentence(.)) %>% 
  rename(sex = gender_src_desc) %>% 
  mutate(race = case_when(
    race_cr_src_desc_1 == "White"    ~ "White",
    !is.na(race_cr_src_desc_1)       ~ "Other"
  ), .before = race_cr_src_desc_1) %>% 
  # For BMI and BSA
  mutate(weight_kg = weight / 2.205,
         height_m = height / 39.37,
         height_cm = height * 2.54,
         bmi = weight_kg / (height_m * height_m)) %>%
  mutate(bmi_cat = case_when(
    bmi < 25                    ~ "<25",
    bmi >= 25 &
      bmi < 30                  ~ "≥25 - <30",
    bmi >= 30                   ~ "≥30"
  )) %>%
  mutate(bmi_cat = factor(bmi_cat, levels = c("<25", "≥25 - <30", "≥30"))) %>% 
  mutate(bsa_du_bois = 0.007184 * (height_cm^0.725) * (weight_kg^0.425)) %>% 
  mutate(bsa_cat = case_when(
    sex == "Male" & 
      bsa_du_bois > 1.9         ~ "> average",
    sex == "Female" & 
      bsa_du_bois > 1.6         ~ "> average",
    sex == "Male" & 
      bsa_du_bois <= 1.9        ~ "≤ average",
    sex == "Female" & 
      bsa_du_bois <= 1.6        ~ "≤ average"
  )) %>% 
  # Survival
  mutate(os_time_from_treatment_months = interval(start = first_treatment_dt, end = date_of_last_contact) /
           duration(n=1, unit = "months")) %>% 
  mutate(os_event = case_when(
    vital_status == "ALIVE"     ~ 0,
    vital_status == "DEAD"      ~ 1
  ))
  
  
# Save
write_csv(lymphoma_data,
          paste0(here::here(), 
                 "/data/processed data",
                 "/Lymphoma_Data_",
                 str_remove_all(today(), "-"), ".csv"))
write_rds(lymphoma_data, 
          paste0(here::here(), 
                 "/data/processed data",
                 "/Lymphoma_Data_",
                 str_remove_all(today(), "-"), ".rds"))

write_csv(lymphoma_data, 
          paste0(path, "/ProcessedData",
                 "/Lymphoma_Data_",
                 str_remove_all(today(), "-"), ".csv"))


# End data cleaning

