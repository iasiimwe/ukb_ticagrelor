# Load libraries
# --------------
library(data.table) # Reading tables faster
library(tidyverse) # Includes library(dplyr) and tidyr #data processing
library(lubridate) # Data processing
library(readxl) # Reading Excel
library(survival) # for Cox Proportional-Hazards Model i.e. for computing survival analyses
library(survminer) # for visualizing survival analysis results
source("helper_functions.R")

# Import primary care data
# -------------------------
gp_clinical.df <- fread("C:/Users/asiimwe/Documents/UKBiobank/gp_clinical.txt.gz") %>%
  select(eid, event_dt, read_2, read_3) %>%
  rename(f.eid = eid)

# Import converted UK Biobank dataset
# -----------------------------------
bd <- as_tibble(fread("C:/Users/asiimwe/Documents/UKBiobank/thiazides/bd_enc.csv.gz")) %>%
  filter(f.eid %in% gp_clinical.df$f.eid) # Limit to those with primary care data
length(unique(bd$f.eid)) # 230060
dict <- as_tibble(fread("C:/Users/asiimwe/Documents/UKBiobank/thiazides/Data_Dictionary_Showcase.csv")) # import dictionary

# Get patients with ischemic heart disease (IHD)
# ----------------------------------------------
icd10_codes <- c("I20", "I21", "I22", "I23", "I24", "I25") # provided
# with_ihd <- first_occurrence(icd10_codes, bd, gp_clinical.df, dict)
# write_rds(with_ihd, "with_ihd.rds")  # remove names and this intermediate step.
with_ihd <- read_rds("with_ihd.rds")
with_ihd <- with_ihd[[1]] %>%
  rbind(with_ihd[[2]]) %>%
  rbind(with_ihd[[3]]) %>%
  rbind(with_ihd[[4]]) %>%
  unique()
length(unique(with_ihd$f.eid)) # 28332

# DRUK Phenotype Library additional codes
# ----------------------------------------
# Download file from the caliber/HDRUK Phenotype Library  (Coronary Heart Disease (CHD), https://phenotypes.healthdatagateway.org/phenotypes/PH1027/version/2263/detail/)
used_codes <- fread(paste0("occurrence_results_", paste0(icd10_codes, collapse = "_"),".csv")) 
caliber <- fread("phenotype_PH1027_ver_2263_concepts_20230209T143624.csv") %>%
  select(code, description, coding_system)  # No read v3 or self reported
used_codes_list <- used_codes %>%
  select(icd9_codes_short, icd10_codes_short, read_2_codes_short)

# ICD-9
caliber_icd9 <- caliber %>% filter(coding_system == "ICD9 codes") %>% pull(code)
not_used_icd9 <- caliber_icd9[!caliber_icd9 %in% unlist(str_split(used_codes_list$icd9_codes_short[[2]], ", "))]
      # "410"  "411"  "414"  "4292" "4295" "4296" "4297"
caliber$description[caliber$coding_system == "ICD9 codes"][!caliber_icd9 %in% unlist(str_split(used_codes_list$icd9_codes_short[[2]], ", "))]
# related_codes <- code_lkup_conditions_2(not_used_icd9, "icd9") 
# additional_icd9 <- first_occurrence(related_codes$icd9[[1]], bd, gp_clinical.df, dict, codes_input_type = "icd9")
# Should capture even the read codes
# write_rds(additional_icd9, "additional_icd9.rds")
additional_icd9 <- read_rds("additional_icd9.rds")

# ICD-10
caliber_icd10 <- caliber %>% filter(coding_system == "ICD10 codes") %>% pull(code) %>% substr(1, 3)
not_used_icd10 <- caliber_icd10[!caliber_icd10 %in% unlist(str_split(used_codes_list$icd10_codes_short[[2]], ", "))] # 0

# Read
caliber_read <- caliber %>% filter(coding_system == "Read codes v2") %>% pull(code) %>% substr(1, 5)
not_used_read <- caliber_read[!caliber_read %in% unlist(str_split(used_codes_list$read_2_codes_short[[4]], " "))]
caliber$description[caliber$coding_system == "Read codes v2"][!caliber_read %in% unlist(str_split(used_codes_list$read_2_codes_short[[4]], " "))]
# related_codes <- code_lkup_conditions_2(not_used_read, "read_2") 
# additional_read <- first_occurrence(not_used_read, bd, gp_clinical.df, dict, codes_input_type = "read_2")
# Should capture even the read codes
# write_rds(additional_read, "additional_read.rds")
additional_read <- read_rds("additional_read.rds")

# Compare 'with_ihd' and the additional patients to see if there are new hits
length(unique(additional_icd9$f.eid[!additional_icd9$f.eid %in% with_ihd$f.eid]))  # 419 
length(unique(additional_read$f.eid[!additional_read$f.eid %in% with_ihd$f.eid]))  # 1567
additional <- unique(c(additional_icd9$f.eid, additional_read$f.eid))
length(additional[!additional %in% with_ihd$f.eid])  # 1733/28332 or 6.1%

# Participants total
with_ihd <- with_ihd %>%
  rbind(additional_icd9) %>%
  rbind(additional_read) %>%
  unique()
length(unique(with_ihd$f.eid))  # 30065 all participants

## Additional fields (https://pubmed.ncbi.nlm.nih.gov/35381005/)
# -------------------------------------------------------------
# Data-Field 41272 (Operative procedures (OPCS4)): CAD Operative Procedures 
# Data-Field 41282, Description:	Date of first operative procedure - OPCS4 (encoded using Data-Coding 240)
Codings <- fread("Codings.csv.gz") %>% #load Codings
  subset(Coding == 240) %>%
  select(Value, Meaning) %>% 
  filter(str_detect(tolower(.$Meaning), "coronary arter(y|ies)|arter(y|ies) into heart|coronary angioplasty|coronary thrombolysis|transluminal angioscopy")) %>%
  mutate(Meaning = "yes")
bd_data <- bd %>%
  select(contains(c("f.eid", "f.53.", "f.41272.", "f.41282."))) %>%
  ukb_reshape_long(dict) %>%
  rename(date = `Date of attending assessment centre`,
         Value = `Operative procedures - OPCS4`,
         date_first_diagnosed = `Date of first operative procedure - OPCS4`) %>%
  left_join(Codings) %>%
  filter(str_detect(.$Meaning, "yes")) %>%
  mutate(year_first_diagnosed = as.integer(year(date_first_diagnosed))) %>%
  select(f.eid, date_first_diagnosed, year_first_diagnosed)
length(unique(bd_data$f.eid[!bd_data$f.eid %in% with_ihd$f.eid])) # 80 additional participants
with_ihd <- with_ihd %>%
  rbind(bd_data) %>%
  unique()
length(unique(with_ihd$f.eid)) # 30145 all participants

# Data-Field 20004, Self-report operative procedures (verbal interview), 92 is the operation year/age first occurred
# Obtained Self-report CABG, Self-report PCI, Self-report angiogram	# Encoded using Data-Coding 5.
# 1070             coronary angioplasty (ptca) +/- stent (obtained by Codings %>% filter(str_detect(tolower(.$Meaning), "coronary")))
# The combination of coronary angioplasty with stenting is usually referred to as percutaneous coronary intervention (PCI).
# 1095      coronary artery bypass grafts (cabg)
# 1514      coronary angiogram
Codings <- fread("Codings.csv.gz") %>% #load Codings
  subset(Coding == 5) %>%
  select(Value, Meaning) %>%
  filter(str_detect(.$Value, "1070|1095|1514")) %>%
  mutate(Meaning = "yes",
         Value = as.integer(Value))
patient_dob <- select(bd, contains(c("f.eid", "f.53.", "f.34.")))  %>%
  ukb_reshape_long(dict) %>%
  filter(I == ".0.") %>%
  rename(year_dob = `Year of birth`) %>%
  select(f.eid, year_dob)
bd_data <- bd %>%
  select(contains(c("f.eid", "f.53.", "f.92.", "f.20004."))) %>%
  ukb_reshape_long(dict) %>%
  rename(date = `Date of attending assessment centre`,
         age_year_first_diagnosed = `Operation year/age first occurred`,
         Value = `Operation code`) %>%
  left_join(Codings) %>%
  filter(str_detect(.$Meaning, "yes")) %>%
  left_join(patient_dob) %>%
  mutate(year_first_diagnosed = ifelse(nchar(age_year_first_diagnosed) == 4,
                                       age_year_first_diagnosed,
                                       year_dob + age_year_first_diagnosed),
         date_first_diagnosed = as.Date(NA)) %>%
  select(f.eid, date_first_diagnosed, year_first_diagnosed)
length(unique(bd_data$f.eid[!bd_data$f.eid %in% with_ihd$f.eid])) # 260 additional participants
with_ihd <- with_ihd %>%
  rbind(bd_data) %>%
  unique()
length(unique(with_ihd$f.eid)) # 30405 all participants

# Data-Field 6150 (Vascular/heart problems diagnosed by doctor)
# Heart attack diagnosed by doctor (self-report), Angina diagnosed by doctor (self-report)
# Field 3627,	Age angina diagnosed; Field 3894,	Age heart attack diagnosed
bd_data <- bd %>%
  select(contains(c("f.eid", "f.53.", "f.6150.", "f.3627.", "f.3894.")))  %>%
  ukb_reshape_long(dict) %>%
  rename(date = `Date of attending assessment centre`,
         condition = `Vascular/heart problems diagnosed by doctor`,
         angina_age = `Age angina diagnosed`,
         heart_attack_age = `Age heart attack diagnosed`) %>%
  filter(condition == "Angina" | condition == "Heart attack") 
no_dates_but_with_ihd <- bd_data %>%
  filter(is.na(angina_age) == TRUE & is.na(heart_attack_age) == TRUE) %>% # To exclude from analysis (including only ihd patients so no need to worry about these anymore)
  select(f.eid)
bd_data <- bd_data %>%
  filter(!f.eid %in% no_dates_but_with_ihd$f.eid) %>%
  rowwise() %>%
  mutate(condition_age = ifelse(is.na(angina_age) == TRUE, min(angina_age, heart_attack_age, na.rm = TRUE), angina_age)) %>%
  group_by(f.eid) %>% 
  summarise(age_first_diagnosed = min(condition_age, na.rm = TRUE)) %>%
  left_join(patient_dob) %>%
  mutate(year_first_diagnosed = year_dob + age_first_diagnosed,
         date_first_diagnosed = as.Date(NA)) %>%
  select(f.eid, date_first_diagnosed, year_first_diagnosed)
length(unique(bd_data$f.eid[!bd_data$f.eid %in% with_ihd$f.eid])) # 156 additional participants
with_ihd <- with_ihd %>%
  rbind(bd_data) %>%
  unique()
length(unique(with_ihd$f.eid)) # 30561 all participants

# Data-Field 42000, Description:	Date of myocardial infarction
# Under Category 44, Health-related outcomes, Algorithmically-defined outcomes, Myocardial infarction outcomes
bd_data <- bd %>%
  select(contains(c("f.eid", "f.53.", "f.42000.")))  %>%
  ukb_reshape_long(dict) %>%
  rename(date = `Date of attending assessment centre`,
         event_dt = `Date of myocardial infarction`) %>%
  filter(is.na(event_dt) == FALSE) %>%
  group_by(f.eid) %>%
  summarise(date_first_diagnosed = min(event_dt, na.rm = TRUE)) %>%
  mutate(year_first_diagnosed = as.integer(year(date_first_diagnosed)))
length(unique(bd_data$f.eid[!bd_data$f.eid %in% with_ihd$f.eid])) # 3 additional participants
with_ihd <- with_ihd %>%
  rbind(bd_data) %>%
  unique()
length(unique(with_ihd$f.eid)) # 30564 all participants

# Get earliest dates/years
# ------------------------
with_ihd_dates <- with_ihd %>%
  filter(is.na(date_first_diagnosed) == FALSE) %>%
  group_by(f.eid) %>%
  summarise(date_first_diagnosed = min(date_first_diagnosed, na.rm = TRUE)) # 29,566 with valid dates
with_ihd <- with_ihd %>%
  group_by(f.eid) %>%
  summarise(year_first_diagnosed = min(year_first_diagnosed, na.rm = TRUE)) %>%
  left_join(with_ihd_dates) 

# Get prescriptions data
# ----------------------
gp_scripts.df <- fread("C:/Users/asiimwe/Documents/UKBiobank/gp_scripts.txt.gz") %>%
  select(-data_provider) %>%
  filter(eid %in% with_ihd$f.eid) %>%
  rename(f.eid = eid) %>%
  mutate(read_2 = tolower(substr(read_2, 1, 5)),
         drug_name = tolower(drug_name),
         bnf_code = tolower(bnf_code))
length(unique(gp_scripts.df$f.eid)) # 29630 participants

# Patients with at least one prescription
drugs <- c("ticagrelor", "clopidogrel", "prasugrel", "aspirin")
# drugs_ids <- number_on_drugs(drugs, gp_scripts.df) 
# write_rds(drugs_ids, "drugs_ids_tica.rds")
drugs_ids <- read_rds("drugs_ids_tica.rds")

# Those taking ticagrelor
with_ihd <- with_ihd %>%
  filter(f.eid %in% drugs_ids[[1]]$f.eid) 

# Check earliest ticagrelor prescription
earliest_drug_date <- drugs_ids[[1]] %>%
  group_by(f.eid) %>%
  summarize(earliest_drug_date = min(as.Date(issue_date, format="%d/%m/%Y"), na.rm = TRUE))
with_ihd <- with_ihd %>%
  left_join(earliest_drug_date) %>%
  filter(date_first_diagnosed <= earliest_drug_date) 

# Obtain registration records
gp_registrations.df <- fread("C:/Users/asiimwe/Documents/UKBiobank/gp_registrations.txt.gz") %>%
  rename(f.eid = eid) %>%
  filter(f.eid %in% with_ihd$f.eid) %>%
  mutate(reg_date = as.Date(reg_date, format="%d/%m/%Y")) %>%
  filter(is.na(reg_date) == FALSE) %>%
  group_by(f.eid) %>%
  summarise(earliest_reg_date = min(reg_date)) # 754 participants
with_ihd <- gp_registrations.df %>%
  left_join(with_ihd) %>%
  mutate(date_diff = as.numeric(earliest_drug_date - earliest_reg_date)) %>%
  filter(date_diff > 730.5) # n = 729, excluded those with ticagrelor purchase within 2 years of registration to ensure 2-year washout

# Get genotype data 
# ------------------
# Individuals to keep, when extracting genotype
bd_keep <- tibble(FID = unique(with_ihd$f.eid), IID = unique(with_ihd$f.eid))
names(bd_keep) <- NULL
write.table(bd_keep, "tica.individuals.txt", sep = " ", row.names = FALSE, quote = FALSE) # file to use on LINUX

# ON LINUX
# cd /pub36/iasiimwe/UKBB/GWAS/TICA
# /pub36/iasiimwe/plink2 --bgen /pub29/andrewt/ukbb_v2/ukb_imp_chr7_v3.bgen ref-first --sample /pub36/iasiimwe/ukbb.sample --keep tica.individuals.txt --snps rs35599367,rs776746 --recode vcf --out chr7_ticageclor_snps
# /pub36/iasiimwe/plink2 --bgen /pub29/andrewt/ukbb_v2/ukb_imp_chr19_v3.bgen ref-first --sample /pub36/iasiimwe/ukbb.sample --keep tica.individuals.txt --snp rs3093135 --recode vcf --out chr19_ticageclor_snps
# for i in 7 19; do /pub36/iasiimwe/plink1.9/plink --vcf chr${i}_ticageclor_snps.vcf --recode tab --out chr${i}_ticageclor_snps; done

## CYP4F2 (rs3093135)
snp_rs3093135 <- fread("chr19_ticageclor_snps.ped") %>%
  select(V1, V7) %>%
  filter(V7 != "0 0") %>%
  mutate(snp_rs3093135 = case_when(V7 == "A A" ~ 0, V7 == "T A" ~ 1, V7 == "T T" ~ 2)) %>%  # Dominant mode of inheritance
  select(-V7) %>%
  rename(f.eid = V1)
table(snp_rs3093135$snp_rs3093135)
# 0   1   2 
# 509 185  21

## CYP3A5 (rs776746, *3) and CYP3A4 (rs35599367, *22)
chr_7_snps <- fread("chr7_ticageclor_snps.ped") %>%
  select(V1, V7, V8) %>%
  mutate(snp_rs776746 = case_when(V7 == "C C" ~ 0, V7 == "T C" ~ 1, V7 == "T T" ~ 2),
         snp_rs35599367 = case_when(V8 == "G G" ~ 0, V8 == "A G" ~ 1, V8 == "A A" ~ 2)) %>%
  select(-V7, -V8) %>%
  rename(f.eid = V1)
table(chr_7_snps$snp_rs776746)
# 0   1   2 
# 605 105   7 
table(chr_7_snps$snp_rs35599367)
# 0   1   2 
# 641  75   1 

# Combine SNPs (exclude the two individuals without snp_rs3093135, so left_join)
ticageclor_snps <- left_join(snp_rs3093135, chr_7_snps) %>%
  as_tibble()

# Join with ihd dataset and keep index date
with_ihd <- with_ihd %>%
  rename(index_date = earliest_drug_date) %>%
  select(f.eid, index_date) 
with_ihd <- ticageclor_snps %>%
  left_join(with_ihd) # 715 participants with genotype data for the 3 SNPs

# Ticagrelor exposure
ticagrelor_prescriptions <- drugs_ids[[1]]
ticagrelor_prescriptions$quantity <- get_quantity_ticagrelor(ticagrelor_prescriptions$quantity)
ticagrelor_prescriptions <- ticagrelor_prescriptions %>%
  filter(quantity != 0)

# Other antiplatelets
length(unique(drugs_ids[[2]]$f.eid[drugs_ids[[2]]$f.eid %in% with_ihd$f.eid])) # 209 patients got clopidogrel
length(unique(drugs_ids[[3]]$f.eid[drugs_ids[[3]]$f.eid %in% with_ihd$f.eid])) # 14 patients got prasugrel
other_antiplatelets_df <- drugs_ids[[2]] %>%
  rbind(drugs_ids[[3]]) %>%
  filter(quantity != "-1.000") %>%
  filter(quantity != "-1 INVALID") %>%
  filter(quantity != "") %>%
  filter(quantity != "0") %>%
  filter(quantity != "0.000") %>%
  filter(quantity != "See Dosage For Quantity") %>%
  filter(quantity != "1.4E-45 tablet(s)") %>%
  select(f.eid, issue_date)

# Obtain follow-up days and exclude patients who took antiplatelets two years within 'index' date
patient_list <- unique(with_ihd$f.eid)
for (i in seq_along(patient_list)) {
  patient <- patient_list[[i]]
  
  # Aim is to get the exposure period (already have index date, default end_date is "2017-09-20")
  patient_df <- with_ihd %>%
    filter(f.eid == patient) %>%
    mutate(end_date = ymd("2017-09-20"))
  
  # ticagrelor prescriptions, given twice a day according to https://bnf.nice.org.uk/drugs/ticagrelor/#indications-and-dose
  ticagrelor_prescription <- ticagrelor_prescriptions %>%
    filter(f.eid == patient) %>%
    select(-drug_name) %>%
    mutate(issue_date = as.Date(issue_date, format="%d/%m/%Y"),
           lead_time = as.numeric(lead(issue_date) - issue_date),
           extra_tablets = quantity/2 - lead_time,
           cummulative_tabs = cumsum(extra_tablets)) %>%
    filter(is.na(cummulative_tabs) == FALSE) %>%
    arrange(issue_date)
  
  if (TRUE %in% (ticagrelor_prescription$cummulative_tabs < -31)) { # all dispensed packages were consumed and no new package was dispensed within 30 days
    patient_df$end_date <- ticagrelor_prescription$issue_date[ticagrelor_prescription$cummulative_tabs < -31][1]
  }
  
  # clopidogrel/prasugrel
  other_antiplatelets <- other_antiplatelets_df %>%
    filter(f.eid == patient) %>%
    mutate(issue_date = as.Date(issue_date, format="%d/%m/%Y")) %>%
    filter(issue_date > (patient_df$index_date - years(2)) & issue_date < patient_df$end_date) 
  # Antiplatelets should be issued at earliest within two years of index date and before current end date
  # A patient with antiplatelets issued before index date will be removed later.
  patient_df$before_index <- "no"
  if (nrow(other_antiplatelets) > 0) {
    if(min(other_antiplatelets$issue_date, na.rm = TRUE) < (patient_df$index_date)) patient_df$before_index <- "yes" 
    patient_df$end_date <- min(other_antiplatelets$issue_date, na.rm = TRUE)
  } 
  
  if (i == 1) with_ihd_2 <- patient_df else with_ihd_2 <- rbind(with_ihd_2, patient_df)
  message(paste(round(i/length(patient_list) * 100, 3), "% complete",sep = ""))
}
with_ihd_2 <- with_ihd_2 %>%
  filter(before_index == "no") %>% # exclude patients who took other antiplatelets at least two years before 'index' date
  filter(end_date - index_date != 0) # exclude patients who have zero days follow-up 
     # 617 participants 
  
# Limit datasets to included participants
gp_clinical.df <- gp_clinical.df %>%
  filter(f.eid %in% with_ihd_2$f.eid)
bd <- bd  %>%
  filter(f.eid %in% with_ihd_2$f.eid)

## Get patients with bleeds
# --------------------------
icd10_codes <- c("D62", "I60", "I61", "I62", "N02", "R04", "R31", "R58") # provided
# with_bleeding_3 <- first_occurrence(icd10_codes, bd, gp_clinical.df, dict, first_only_input = "FALSE")
# write_rds(with_bleeding_3, "with_bleeding_3.rds")
with_bleeding_3 <- read_rds("with_bleeding_3.rds")

icd10_codes <- c("D500", "D683", "I690", "I691", "I692", "I850", "J942", "K221", "K223", "K226", "K250", "K252", 
                 "K254", "K256", "K260", "K262", "K264", "K266", "K270", "K272", "K274", "K276", "K280", "K282", 
                 "K284", "K286", "K290", "K625", "K631", "K633", "K920", "K921", "K922", "S062", "S063", "S064", 
                 "S065", "S066", "S068")
# with_bleeding_4 <- first_occurrence(icd10_codes, bd, gp_clinical.df, dict, first_only_input = "FALSE", n_icd_input = 4, include_icd9 = FALSE)
# write_rds(with_bleeding_4, "with_bleeding_4.rds")
with_bleeding_4 <- read_rds("with_bleeding_4.rds")

with_bleeding_3 <- with_bleeding_3[[1]] %>%
  rbind(with_bleeding_3[[2]]) %>%
  rbind(with_bleeding_3[[3]]) %>%
  rbind(with_bleeding_3[[4]]) %>%
  filter(f.eid %in% with_ihd_2$f.eid) %>%
  unique()
length(unique(with_bleeding_3$f.eid)) # 71 bleeding events

with_bleeding_4 <- with_bleeding_4[[1]] %>%
  rbind(with_bleeding_4[[2]]) %>%
  rbind(with_bleeding_4[[3]]) %>%
  rbind(with_bleeding_4[[4]]) %>%
  filter(f.eid %in% with_ihd_2$f.eid) %>%
  unique()
length(unique(with_bleeding_4$f.eid)) # 165 bleeding events

with_bleeding <- with_bleeding_3 %>%
  rbind(with_bleeding_4)
length(unique(with_bleeding$f.eid)) # 213 bleeding events

# Major bleeding
icd10_codes <- c("D500", "K920", "K921", "S064")
# with_bleeding_major <- first_occurrence(icd10_codes, bd, gp_clinical.df, dict, first_only_input = "FALSE", n_icd_input = 4, include_icd9 = FALSE)
# write_rds(with_bleeding_major, "with_bleeding_major.rds")
with_bleeding_major <- read_rds("with_bleeding_major.rds")

with_bleeding_major <- with_bleeding_major[[1]] %>%
  rbind(with_bleeding_major[[2]]) %>%
  rbind(with_bleeding_major[[3]]) %>%
  rbind(with_bleeding_major[[4]]) %>%
  filter(f.eid %in% with_ihd_2$f.eid) %>%
  unique()
length(unique(with_bleeding_major$f.eid)) # 33 major bleeding events

with_ihd_2$bleed <- "no"
with_ihd_2$bleed_major <- "no"
patient_list <- unique(with_bleeding$f.eid)
for (i in seq_along(patient_list)) {
  patient <- patient_list[[i]]
  with_ihd_patient <- with_ihd_2 %>%
    filter(f.eid == patient)
  # all bleeds
  patient_df <- with_bleeding %>%
    filter(f.eid == patient) %>%
    filter(date_diagnosed >= with_ihd_patient$index_date & date_diagnosed <= with_ihd_patient$end_date)
  if (nrow(patient_df) > 0) with_ihd_2$bleed[with_ihd_2$f.eid == patient] <- "yes"
  # major bleeds
  patient_df <- with_bleeding_major %>%
    filter(f.eid == patient) %>%
    filter(date_diagnosed >= with_ihd_patient$index_date & date_diagnosed <= with_ihd_patient$end_date)
  if (nrow(patient_df) > 0) with_ihd_2$bleed_major[with_ihd_2$f.eid == patient] <- "yes"
  
  message(paste(round(i/length(patient_list) * 100, 3), "% complete",sep = ""))
}
table(with_ihd_2$bleed)
# no yes 
# 498 119  
table(with_ihd_2$bleed_major)
# no yes 
# 601  16  


# Survival analysis 
# -----------------
with_ihd_2 <- with_ihd_2 %>%
  mutate(time_days = as.numeric(end_date - index_date),
         status = case_when(bleed == "no" ~ 1, bleed == "yes" ~ 2),
         status2 = case_when(bleed_major == "no" ~ 1, bleed_major == "yes" ~ 2)
         ) 
summary(with_ihd_2$time_days) # Followup days
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6.0   462.0   772.0   834.7  1254.0  2021.0 

# Add age and sex 
patient_age_sex <- select(bd, contains(c("f.eid", "f.53.", "f.31.", "f.21003.")))  %>%
  ukb_reshape_long(dict) %>%
  filter(I == ".0.") %>%
  rename(sex = Sex,
         age = `Age when attended assessment centre`) %>%
  select(f.eid, age, sex) %>%
  mutate(sex = factor(sex))  
with_ihd_2 <- with_ihd_2 %>%
  left_join(patient_age_sex)

# Aspirin status
asa_df <- rbind(drugs_ids[[4]]) %>%
  filter(quantity != "0 refill disks") %>%
  filter(quantity != "-1 INVALID") %>%
  filter(quantity != "") %>%
  select(f.eid, issue_date)
with_ihd_2$asa <- "no"
patient_list <- unique(with_ihd_2$f.eid)
for (i in seq_along(patient_list)) {
  patient <- patient_list[[i]]
  patient_df <- with_ihd_2 %>%
    filter(f.eid == patient)
  asa <- asa_df %>%
    filter(f.eid == patient) %>%
    mutate(issue_date = as.Date(issue_date, format="%d/%m/%Y")) %>%
    filter(issue_date >= patient_df$index_date & issue_date <= patient_df$end_date)
  if (nrow(asa) > 0) with_ihd_2$asa[with_ihd_2$f.eid == patient] <- "yes"
  message(paste(round(i/length(patient_list) * 100, 3), "% complete",sep = ""))
}
table(with_ihd_2$asa) # None on asa, so don't include in analysis

# Univariable analysis (additive mode) 
summary(coxph(Surv(time_days, status) ~ snp_rs3093135 , data = with_ihd_2)) 
summary(coxph(Surv(time_days, status) ~ snp_rs776746 , data = with_ihd_2)) 
summary(coxph(Surv(time_days, status) ~ snp_rs35599367 , data = with_ihd_2)) 

# Multivariable analysis (additive mode)
res.cox <- coxph(Surv(time_days, status) ~ age + sex + snp_rs3093135 + snp_rs776746 + snp_rs35599367, data = with_ihd_2)
summary(res.cox) 

# Univariable analysis (dominant mode) 
with_ihd_2 <- with_ihd_2 %>%
  mutate(snp_rs3093135 = if_else(snp_rs3093135 == 2, 1, snp_rs3093135),
        snp_rs776746 = if_else(snp_rs776746 == 2, 1, snp_rs776746),
        snp_rs35599367 = if_else(snp_rs35599367 == 2, 1, snp_rs35599367))

summary(coxph(Surv(time_days, status) ~ snp_rs3093135 , data = with_ihd_2)) 
summary(coxph(Surv(time_days, status) ~ snp_rs776746 , data = with_ihd_2)) 
summary(coxph(Surv(time_days, status) ~ snp_rs35599367 , data = with_ihd_2)) 

# Multivariable analysis (dominant mode)
res.cox <- coxph(Surv(time_days, status) ~ age + sex + snp_rs3093135 + snp_rs776746 + snp_rs35599367, data = with_ihd_2)
summary(res.cox) 


# Prepare for GWAS
# ----------------
with_ihd_2 <- with_ihd_2 %>%
  select(f.eid, time_days:sex)

# Add the principal components of genetic ancestry
pcs <- as_tibble(fread("C:/Users/asiimwe/Documents/UKBiobank/thiazides/UKBB_principal_components.csv.gz"))[1:11] #load data
colnames(pcs) <- c("f.eid", paste0("C", 1:10))
with_ihd_2 <- left_join(with_ihd_2, pcs)

bd_pheno <- with_ihd_2 %>%
  mutate(FID = f.eid,
         IID = f.eid
    ) %>%
  select(-f.eid) %>%
  relocate(time_days:C10, .after = last_col())
write.table(bd_pheno, file = paste("pheno_tica.txt", sep = ""), sep = " ", row.names = FALSE, quote = FALSE, na = "NA")

bd_keep <- bd_pheno %>% select(FID, IID)
names(bd_keep) <- NULL
write.table(bd_keep, file = paste("keep_tica.txt", sep = ""), sep = " ", row.names = FALSE, quote = FALSE)


# Run in LINUX
# ------------
# Get imputed data
# for i in {1..22}; do
# /pub36/iasiimwe/plink2 --bgen ukb_imp_chr${i}_v3.bgen ref-first --sample ukbb.sample --keep keep_tica.txt --make-bed --out tica_chr_${i} --geno 0.05 --maf 0.01 --hwe 0.000001 --minimac3-r2-filter 0.3
# done

# Perform GWAS in R (used cluster)
# --------------------------------
# Load packages
library(data.table)
library(tidyverse)
library(bigsnpr)
library(survival)

setwd("/pub36/iasiimwe/UKBB/GWAS/TICA")

## Preprocess the bed files (done one)
# for (i in 1:22) snp_readBed(paste0("tica_chr_", i, ".bed"))

# Phenotype file
phenotype <- fread("pheno_tica.txt")

for (i in 1:22) {
  # Attach the genotype object
  obj.bigSNP <- snp_attach(paste0("tica_chr_", i, ".rds"))
  
  # Assign the genotype to a variable 
  genotype <- obj.bigSNP$genotypes
  
  # We assume the fam order is the same across different chromosomes
  fam.order <- as.data.table(obj.bigSNP$fam) %>%
    select(family.ID) %>%
    rename(FID = family.ID)
  
  # Reformat the phenotype file such that pheno is of the same order as the sample ordering in the genotype file
  pheno <- left_join(fam.order, phenotype)
  
  results <- tibble(chr = rep(i, ncol(genotype)),
                    SNP = rep(NA, ncol(genotype)),
                    BP = rep(NA, ncol(genotype)),
                    allele1 = rep(NA, ncol(genotype)),
                    allele2 = rep(NA, ncol(genotype)),
                    N = rep(NA, ncol(genotype)),
                    events = rep(NA, ncol(genotype)),
                    MAF = rep(NA, ncol(genotype)),
                    P_adj = rep(NA, ncol(genotype)),
                    HR_adj = rep(NA, ncol(genotype)),
                    lower_95_adj = rep(NA, ncol(genotype)),
                    upper_95_adj = rep(NA, ncol(genotype)),
                    dom_P_adj = rep(NA, ncol(genotype)),
                    dom_HR_adj = rep(NA, ncol(genotype)),
                    dom_lower_95_adj = rep(NA, ncol(genotype)),
                    dom_upper_95_adj = rep(NA, ncol(genotype)))
  # The logrank test is the score test from a Cox proportional hazards model, so it makes the same assumptions as the Cox model. The LR test, among the three commonly used tests (the other two being Wald and score) is the gold standard.
  
  for(j in 1:ncol(genotype)) {
    pheno$snp <- genotype[, j]
    pheno_2 <- pheno %>% filter(is.na(snp) == FALSE)
    results$SNP[j] <- obj.bigSNP$map$marker.ID[j]
    results$BP[j] <- obj.bigSNP$map$physical.pos[j]
    results$allele1[j] <- obj.bigSNP$map$allele1[j]
    results$allele2[j] <- obj.bigSNP$map$allele2[j]
    results$N[j] <- nrow(pheno_2)
    results$MAF[j] <- (length(pheno_2$snp[pheno_2$snp == 1]) +  2*length(pheno_2$snp[pheno_2$snp == 2]))/(2*nrow(pheno_2))
    
    # Additive mode
    cox_summary <- summary(coxph(Surv(time_days, status) ~ snp + age + sex + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C10 , data = pheno_2))
    results$events[j] <- cox_summary$nevent
    results$P_adj[j] <- cox_summary$coefficients[1, 5]
    results$HR_adj[j] <- cox_summary$conf.int[1, 1]
    results$lower_95_adj[j] <- cox_summary$conf.int[1, 3]
    results$upper_95_adj[j] <- cox_summary$conf.int[1, 4]
    
    # Dominant mode
    pheno_2 <- pheno_2 %>%
      mutate(snp = if_else(snp == 2, 1, snp))
    cox_summary <- summary(coxph(Surv(time_days, status) ~ snp + age + sex + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C10 , data = pheno_2))
    results$dom_P_adj[j] <- cox_summary$coefficients[1, 5]
    results$dom_HR_adj[j] <- cox_summary$conf.int[1, 1]
    results$dom_lower_95_adj[j] <- cox_summary$conf.int[1, 3]
    results$dom_upper_95_adj[j] <- cox_summary$conf.int[1, 4]
    
    message(paste0("Chr ", i, ": ", round(j/ncol(genotype) * 100, 3), "% complete"))
  }
  write.table(results, file = paste0("tica_chr", i, ".txt"), sep = " ", row.names = FALSE, quote = FALSE, na = "NA")
}

# # After GWAS (LINUX)
# # ------------------
# cd /pub36/iasiimwe/UKBB/GWAS/TICA
# /bin/awk 'NR==1' tica_chr1.txt > tica_Manhattan.txt
# for i in {1..22}; do
# /bin/sed '1d' tica_chr${i}.txt >> tica_Manhattan.txt 
# done
# /pub36/iasiimwe/htslib/htslib-1.3.2/bgzip --threads 32 tica_Manhattan.txt

# source /pub36/iasiimwe/miniconda3/bin/activate renv
# In R
# ----
# Libraries and relevant functions
library(data.table) # reading tables faster
library(tidyverse) # data processing
library(qqman) # for generating Manhattan and QQ plots
source("/pub36/iasiimwe/UKBB/GWAS/relevant_functions.R")
results <- as_tibble(fread(paste("tica_Manhattan.txt.gz", sep = "")))
summary(results$N)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   589     612     616     614     619     619

# All bleeds
# Additive
results2 <- results %>%
  rename(CHR = chr,
         P = P_adj) %>%
  filter(is.na(P) == FALSE & MAF >= 0.01)
write.csv(results2[results2$P < 5e-8,],
          file = paste("tica_additive.csv", sep = ""),
          row.names = FALSE)
nrow(results2) # analyzed_snps: 8653366
nrow(results2[results2$P < 5e-8,])  # genome wide significant: 68
# Manhattan plot
ukb_manhattan_plots(results2, "additive", "adjusted")
# Quantile-quantile (QQ) plot
ukb_qq_gwas_plots(results2, "additive", "adjusted")
# Genomic inflation factor (lambda)
ukb_lambda(results2) # 1.061657

# Dominant
results2 <- results %>%
  rename(CHR = chr,
         P = dom_P_adj) %>%
  filter(is.na(P) == FALSE & MAF >= 0.01)
write.csv(results2[results2$P < 5e-8,],
          file = paste("tica_dominant.csv", sep = ""),
          row.names = FALSE)
nrow(results2) # analyzed_snps: 8587357
nrow(results2[results2$P < 5e-8,])  # genome wide significant: 63
# Manhattan plot
ukb_manhattan_plots(results2, "dominant", "adjusted")
# Quantile-quantile (QQ) plot
ukb_qq_gwas_plots(results2, "dominant", "adjusted")
# Genomic inflation factor (lambda)
ukb_lambda(results2) # 1.02165

# LD analysis
# ------------
# Get SNP file (additive, then dominant)
top_additive <- results %>%
  rename(CHR = chr,
         P = P_adj) %>%
  filter(is.na(P) == FALSE & P < 5e-8 & MAF >= 0.01) %>%
  select(CHR, SNP, MAF, P)
write.table(top_additive,
            file = paste0("top_additive.txt"), 
            sep = " ", 
            row.names = FALSE, 
            quote = FALSE)
top_dominant <- results %>%
  rename(CHR = chr,
         P = dom_P_adj) %>%
  filter(is.na(P) == FALSE & P < 5e-8 & MAF >= 0.01) %>%
  select(CHR, SNP, MAF, P)
write.table(top_dominant, file = paste0("top_dominant.txt"), 
              sep = " ", 
              row.names = FALSE, 
              quote = FALSE)

# SNPs to extract
top_additive <- top_additive %>%
  select(SNP)
top_dominant %>%
  select(SNP) %>%
  rbind(top_additive) %>%
  unique() %>%
  write.table(file = paste0("top_snps.txt"), 
              sep = " ", 
              row.names = FALSE, 
              quote = FALSE)

## Merge the files (used plink1.9, in LINUX)
# /plink1.9/plink \
# --bfile tica_chr_1 \
# --merge-list allfiles_tica.txt \
# --make-bed \
# --out tica
## Note: 1. allfiles_tica.txt contains the files to be merged.
##       2. If SNPs producing an error (missnps) exist, exclude them before merging again.

# for i in {1..22}; do
# /pub36/iasiimwe/plink1.9/plink --bfile tica_chr_${i} --exclude tica-merge.missnp --make-bed --out tica_chr_${i}_2
# done

## Now merge again
# /plink1.9/plink \
# --bfile tica_chr_1_2 \
# --merge-list allfiles_tica2.txt \
# --make-bed \
# --out tica 

## Get near-independent SNPs
# genetic_mode="additive"
# # genetic_mode="dominant"
# /plink1.9/plink \
# --bfile tica \
# --clump-p1 5e-8 \
# --clump-p2 5e-8 \
# --clump-r2 0.5 \
# --clump-kb 250 \
# --clump top_${genetic_mode}.txt \
# --clump-snp-field SNP \
# --clump-field P \
# --out ${genetic_mode}

# Outputs (.clumped) opened in Excel and genomic locations obtained from dbSNP (https://www.ncbi.nlm.nih.gov/snp/)
