# 1. Change to long format function
# ---------------------------------
# Uses the pivot_longer() function {tidyr} to "lengthen" UKBiobank datasets.
# Requires a UKBiobank dataset (bd_enc.csv) converted (using ukbconv) - one participant per row and 
# a dictionary to identify the Field IDs (Data_Dictionary_Showcase.csv).
ukb_reshape_long <- function(data, dict) {
  fields <- data %>% 
    colnames() %>% 
    str_extract_all(., "f.[0-9]*") %>% 
    unlist() %>% 
    str_replace_all(., "f.", "") %>%
    unique() %>%
    parse_integer() %>%
    na.omit()
  fields <- append(fields[!fields %in% 53], 53, after = 0) # remove date field, so that it can be added  
  # to the start of the vector - this enables all assessments to be added using left_join().
  for (k in seq_along(fields)) {
    columns <- unlist(str_extract_all(colnames(data), 
                                      str_c("f.", fields[k], ".[0-9]*.[0-9]*", 
                                            sep=""))) # columns to pivot into longer format
    values_name <- dict %>% filter(FieldID == fields[k]) %>% .$Field
    
    if(fields[k] %in% c(41203, 41271)) {
      names_to_change <- data %>% 
        select(contains(c(str_c("f.", fields[k], ".", sep="")))) %>%
        colnames()
      data <- data %>% mutate_at(names_to_change, as.character)
    }
    
    dataQ <- data %>% 
      select(contains(c("f.eid", str_c("f.", fields[k], ".", sep="")))) %>%
      pivot_longer(cols = all_of(columns), 
                   names_to = "I", 
                   values_to = values_name, 
                   values_drop_na = TRUE
      ) %>%
      arrange(f.eid) %>%
      mutate(I = str_replace(.$I, str_c("f.", fields[k], sep=""), "")) %>% 
      mutate(I = unlist(str_extract(.$I, ".[0-9].")))
    if (k == 1) 
      data_long <- dataQ 
    else 
      data_long <- left_join(data_long, dataQ)   #Join by "f.eid" and "I"
  }
  return(data_long)
}

# Obtaining patients on a specific drug
# -------------------------------------
code_lkup_drugs <- function(drug, path, sheet, column_use, column_to_lower) {
  data <- read_excel(path, sheet)
  data <- data[str_detect(tolower(data[[column_to_lower]]), drug), ]
  code_df <- data %>%
    select(all_of(column_use)) %>%
    mutate(code = tolower(.[[column_use]])) %>%
    select(code) %>%
    unique() %>%
    filter(is.na(code) == FALSE) %>%
    mutate(value = "yes") 
  return(code_df)
}
code_search <- function(code_df, column_use, data = gp_scripts.df) {
  code_hits <- data %>%
    mutate(code = tolower(.[[column_use]])) %>%
    left_join(code_df) %>%
    filter(value == "yes") %>%
    select(f.eid, issue_date, drug_name, quantity) %>%
    unique() 
  return(code_hits)
}
number_on_drugs <- function(drugs, gp_scripts.df, write_summary = TRUE) {
  drugs_ids <- vector("list", length(drugs))
  drugs_results <- tibble(brand_names = rep(NA, length(drugs)),
                          names_hits = rep(NA, length(drugs)),
                          bnf_codes = rep(NA, length(drugs)),
                          bnf_hits = rep(NA, length(drugs)),
                          dmd_codes = rep(NA, length(drugs)),
                          dmd_hits = rep(NA, length(drugs)),
                          read_codes = rep(NA, length(drugs)),
                          read_hits = rep(NA, length(drugs)),
                          total_hits = rep(NA, length(drugs)))
  for (i in seq_along(drugs)) {
    ## Getting the codings for drugs
    # BNF
    drug <- code_lkup_drugs(drugs[[i]], "all_lkps_maps_v3.xlsx", "bnf_lkp", "BNF_Product", "BNF_Chemical_Substance") %>%
      pull(code) %>%
      paste0(., collapse = "|")
    drugs_results$brand_names[[i]] <- gsub("\\|", ", ", drug)
    bnf_code <- code_lkup_drugs(drug, "all_lkps_maps_v3.xlsx", "bnf_lkp", "BNF_Presentation_Code", "BNF_Chemical_Substance") 
    drugs_results$bnf_codes[[i]] <- paste0(bnf_code$code, collapse = ", ")
    
    # DMD
    dmd_code <- code_lkup_drugs(drug, "all_lkps_maps_v3.xlsx", "dmd_lkp", "concept_id", "term") 
    drugs_results$dmd_codes[[i]] <- paste0(dmd_code$code, collapse = ", ")
    
    # Read_v2
    read_2_code <- code_lkup_drugs(drug, "all_lkps_maps_v3.xlsx", "read_v2_drugs_lkp", "read_code", "term_description")
    drugs_results$read_codes[[i]] <- paste0(read_2_code$code, collapse = ", ")
    
    # Drug name search
    name_hits <- gp_scripts.df[str_detect(tolower(gp_scripts.df$drug_name), drug), ] %>%
      select(f.eid, issue_date, drug_name, quantity) %>%
      unique() 
    drugs_results$names_hits[[i]] <- length(unique(name_hits$f.eid))
    
    # BNF code search
    bnf_hits <- code_search(bnf_code, "bnf_code") 
    drugs_results$bnf_hits[[i]] <- length(unique(bnf_hits$f.eid)) 
    
    # DMD code search
    dmd_hits <- code_search(dmd_code, "dmd_code")
    drugs_results$dmd_hits[[i]] <- length(unique(dmd_hits$f.eid))
    
    # Read code search
    read_hits <- code_search(read_2_code, "read_2")
    drugs_results$read_hits[[i]] <- length(unique(read_hits$f.eid))
    
    # Total hits
    all_ids <- unique(c(name_hits$f.eid, bnf_hits$f.eid, dmd_hits$f.eid, read_hits$f.eid))
    drugs_results$total_hits[[i]] <- length(all_ids) 
    drugs_ids[[i]] <- name_hits %>%
      rbind(bnf_hits) %>%
      rbind(dmd_hits) %>%
      rbind(read_hits) %>%
      unique()
    
    message(paste(round(i/length(drugs) * 100, 3), "% complete",sep = ""))
  }
  if (write_summary == TRUE) write.csv(drugs_results, paste0("drugs_results_", substr(paste0(drugs, collapse = "_"), 1, 50),".csv"), row.names = FALSE)
  return(drugs_ids) 
}

# Obtaining first occurrence of a condition 
# -----------------------------------------
# Four sources
# 1. Death register - Fields 40000 (date of death), 40001 (primary cause) and 40002 (contributory/secondary) 
# 2. Hospital inpatient data: Field 41234
# 41202	Diagnoses - main ICD10 #	41262	Date of first in-patient diagnosis - main ICD10
# 41203	Diagnoses - main ICD9	# 41263	Date of first in-patient diagnosis - main ICD9
# 41270	Diagnoses - ICD10	# 41280	Date of first in-patient diagnosis - ICD10
# 41271	Diagnoses - ICD9	# 41281	Date of first in-patient diagnosis - ICD9
# 3. Self-report data - Data-field 20002 (condition) and 20008 (Interpolated year when noncancer illness first diagnosed)
# 4. Primary care
code_lkup_helper <- function(provided_codes, sheet = NULL, columns_to_read, path = "all_lkps_maps_v3.xlsx") {
  if (is.null(sheet) == TRUE) {
    data <- fread(path) %>%
      select(all_of(columns_to_read)) %>%
      unique()
  } else {
    data <- read_excel(path, sheet) %>%
      select(all_of(columns_to_read)) %>%
      unique()
  }
  
  colnames(data) <- c("old_code", "new_code")
  data <- data %>%
    filter(is.na(old_code) == FALSE) %>%
    filter(is.na(new_code) == FALSE) %>%
    filter(old_code != "UNDEF") %>%
    filter(new_code != "UNDEF")
  
  new_codes <- provided_codes %>%
    left_join(data) %>%
    select(new_code, value) %>%
    filter(is.na(new_code) == FALSE)
  return(new_codes)
}
code_lkup_conditions <- function(codes_tb, code_type, output, n_icd) { # Read 2 (Can also use "coding1834.tsv", TRUD mapping of Read 2 into 3-character ICD10 (https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=1834)) but it has fewer rows
  if (code_type == "icd9") {
    new_codes_tb <- switch(output,
           icd9 = codes_tb %>% rename(new_code = old_code),
           icd10 = ifelse(n_icd == 3,
                          code_lkup_helper(codes_tb, sheet = NULL, c("coding", "meaning"), path = "coding1836.tsv"),
                          code_lkup_helper(codes_tb, "icd9_icd10", c("ICD9", "ICD10"))),
           read_2 = code_lkup_helper(codes_tb, "read_v2_icd9", c("icd9_code", "read_code")),
           read_3 = code_lkup_helper(codes_tb, "read_ctv3_icd9", c("icd9_code", "read_code")),
           opcs4 = NULL,
           stop("Unknown output!"))
  } else if (code_type == "icd10") {
    new_codes_tb <- switch(output,
           icd9 = ifelse(n_icd == 3,
                         code_lkup_helper(codes_tb, sheet = NULL, c("meaning", "coding"), path = "coding1836.tsv"),
                         code_lkup_helper(codes_tb, "icd9_icd10", c("ICD10", "ICD9"))),
           icd10 = codes_tb %>% rename(new_code = old_code),
           read_2 = ifelse(n_icd == 3,
                           code_lkup_helper(codes_tb, sheet = NULL, c("icd10_3char", "read_2"), path = "GP_Read2_FO_mappings_Jan2022.tsv"),
                           code_lkup_helper(codes_tb, "read_v2_icd10", c("icd10_code", "read_code"))),
           read_3 = ifelse(n_icd == 3,
                           code_lkup_helper(codes_tb, sheet = NULL, c("icd10_3char", "read_3"), path = "GP_Read3_FO_mappings_Jan2022.tsv"),
                           code_lkup_helper(codes_tb, "read_ctv3_icd10", c("icd10_code", "read_code"))),
           opcs4 = NULL,
           stop("Unknown output!"))
  } else if (code_type == "read_2") {
    new_codes_tb <- switch(output,
           icd9 = code_lkup_helper(codes_tb, "read_v2_icd9", c("read_code", "icd9_code")),
           icd10 = ifelse(n_icd == 3,
                          code_lkup_helper(codes_tb, sheet = NULL, c("read_2", "icd10_3char"), path = "GP_Read2_FO_mappings_Jan2022.tsv"),
                          code_lkup_helper(codes_tb, "read_v2_icd10", c("read_code", "icd10_code"))),
           read_2 = codes_tb %>% rename(new_code = old_code),
           read_3 = code_lkup_helper(codes_tb, "read_v2_read_ctv3", c("READV2_CODE", "TERMV3_CODE")),
           opcs4 = code_lkup_helper(codes_tb, "read_v2_opcs4", c("read_code", "opcs_4.2_code")),
           stop("Unknown output!"))
  } else if (code_type == "read_3") {
    new_codes_tb <- switch(output,
           icd9 = code_lkup_helper(codes_tb, "read_ctv3_icd9", c("read_code", "icd9_code")),
           icd10 = ifelse(n_icd == 3,
                          code_lkup_helper(codes_tb, sheet = NULL, c("read_3", "icd10_3char"), path = "GP_Read3_FO_mappings_Jan2022.tsv"),
                          code_lkup_helper(codes_tb, "read_ctv3_icd10", c("read_code", "icd10_code"))),
           read_2 = code_lkup_helper(codes_tb, "read_v2_read_ctv3", c("TERMV3_CODE", "READV2_CODE")),
           read_3 = codes_tb %>% rename(new_code = old_code),
           opcs4 = code_lkup_helper(codes_tb, "read_ctv3_opcs4", c("read_code", "opcs4_code")),
           stop("Unknown output!"))
  } else {
    new_codes_tb <- switch(output,
           icd9 = NULL,
           icd10 = NULL,
           read_2 = code_lkup_helper(codes_tb, "read_v2_opcs4", c("opcs_4.2_code", "read_code")),
           read_3 = code_lkup_helper(codes_tb, "read_ctv3_opcs4", c("opcs4_code", "read_code")),
           opcs4 = codes_tb %>% rename(new_code = old_code),
           stop("Unknown output!"))
  }
  
  if (is_tibble(new_codes_tb) == TRUE) return(new_codes_tb) else return(tibble(new_code = unlist(new_codes_tb), value = "yes"))
}
code_full_description <- function(path, sheet, column_use, code_df, n_icd = 3) {
  data <- read_excel(path, sheet = sheet) %>%    
    select(all_of(column_use))
  colnames(data) <- c("new_code", "description")
  if (sheet == "icd10_lkp" || sheet == "icd9_lkp") data$new_code <- substr(data$new_code, 1, n_icd)
  data <- data %>%
    left_join(code_df) %>%
    filter(value == "yes") %>%
    mutate(code_description = paste(new_code, description, sep = ": ")) %>%
    pull(code_description) %>% 
    paste0(., collapse = ", ")
  return(data)
}
pc_mapping_read <- function(codes, column_use, first_only_mapping, data_gp = gp_clinical.df, data_bd = bd, data_dict = dict) {
  patient_dob <- select(data_bd, contains(c("f.eid", "f.53.", "f.34.")))  %>%
    ukb_reshape_long(data_dict) %>%
    filter(I == ".0.") %>%
    rename(year_dob = `Year of birth`) # Date of birth for removing invalid dates
  data <- data_gp %>%
    mutate(new_code = .[[column_use]]) %>%
    left_join(codes) %>%
    select(f.eid, event_dt, value) %>%
    filter(str_detect(.$value, "yes")) %>%
    filter(is.na(event_dt) == FALSE) %>%
    mutate(event_dt = as.Date(event_dt, format="%d/%m/%Y")) %>%
    left_join(patient_dob) %>%
    filter(year(event_dt) >= year_dob)  # Exclude patients with event year less than participants year of birth 

  if (first_only_mapping == "TRUE") {
    data <- data %>%
      group_by(f.eid) %>%
      summarise(date_first_diagnosed = suppressWarnings(min(event_dt)),
                year_first_diagnosed = as.integer(year(date_first_diagnosed))) %>%
      as_tibble()
  } else {
    data <- data %>%
      select(f.eid, event_dt) %>%
      rename(date_diagnosed = event_dt) %>%
      mutate(year_diagnosed = as.integer(year(date_diagnosed))) %>%
      as_tibble()
  }
  return(data)
}
ons_mapping_icd10 <- function(codes, column_use, new_colnames, first_only_mapping, data_bd = bd, data_dict = dict, n_icd) {
  data <- data_bd %>%
    select(contains(c(c("f.eid", "f.53."), column_use))) %>%
    ukb_reshape_long(dict) 
  colnames(data) <- new_colnames 
  data <- data %>%
    filter(is.na(diagnosis_date) == FALSE) %>% # remove those who haven't died (don't have a cause of death)
    filter(year(diagnosis_date) >= 2006) %>% # Recruited participants cannot have died earlier that year of recruitment
    mutate(new_code = substr(death_cause, 1, n_icd)) %>%
    left_join(codes) %>%
    select(-death_cause, -new_code) %>%
    rename(value1 = value) %>%
    mutate(new_code = substr(death_other_cause, 1, n_icd)) %>%
    left_join(codes) %>%
    mutate(value = paste0(value1, value)) %>%
    select(f.eid, I, date, diagnosis_date, value) %>%
    filter(str_detect(.$value, "yes")) %>%
    filter(is.na(diagnosis_date) == FALSE)

  if (first_only_mapping == "TRUE") {
    data <- data %>%
      group_by(f.eid) %>%
      summarise(date_first_diagnosed = suppressWarnings(min(diagnosis_date)),
                year_first_diagnosed = as.integer(year(date_first_diagnosed)))
  } else {
    data <- data %>%
      select(f.eid, diagnosis_date) %>%
      rename(date_diagnosed = diagnosis_date) %>%
      mutate(year_diagnosed = as.integer(year(date_diagnosed)))
  }
  return(data)
}
pc_mapping_icd <- function(codes, column_use, first_only_mapping, new_colnames = c("f.eid", "I", "date", "new_code", "diagnosis_date"), 
                           data_bd = bd, data_dict = dict, code_type = "icd10", n_icd) {
  data <- data_bd %>%
    select(contains(c(c("f.eid", "f.53."), column_use))) %>%
    ukb_reshape_long(dict) 
  colnames(data) <- new_colnames 
  if (code_type == "icd10") data <- mutate(data, new_code = substr(new_code, 1, n_icd))
  data <- data %>%
    left_join(codes) %>%
    filter(str_detect(.$value, "yes")) %>%
    filter(is.na(diagnosis_date) == FALSE)
    
  if (first_only_mapping == "TRUE") {
    data <- data %>%
      group_by(f.eid) %>%
      summarise(date_first_diagnosed = suppressWarnings(min(diagnosis_date)),
                year_first_diagnosed = as.integer(year(date_first_diagnosed)))
  } else {
    data <- data %>%
      select(f.eid, diagnosis_date) %>%
      rename(date_diagnosed = diagnosis_date) %>%
      mutate(year_diagnosed = as.integer(year(date_diagnosed)))
  }
  
  return(data)
}
first_occurrence <- function(provided_codes, bd_input, gp_clinical.df_input, dict_input, first_only_input = "TRUE", 
                             codes_input_type = "icd10", n_icd_input = 3, write_summary = TRUE, include_icd9 = TRUE) {
  n_rows <- 5
  occurrence_ids <- vector("list", n_rows-1) 
  occurrence_results <- data.frame(icd9_codes_short = rep(NA, n_rows),
                                   icd9_codes_full = rep(NA, n_rows),
                                   icd9_hits = rep(NA, n_rows),
                                   icd10_codes_short = rep(NA, n_rows),
                                   icd10_codes_full = rep(NA, n_rows),
                                   icd10_hits = rep(NA, n_rows),
                                   read_2_codes_short = rep(NA, n_rows),
                                   read_2_codes_full = rep(NA, n_rows),
                                   read_2_hits = rep(NA, n_rows),
                                   read_3_codes_short = rep(NA, n_rows),
                                   read_3_codes_full = rep(NA, n_rows),
                                   read_3_hits = rep(NA, n_rows),
                                   self_reported_codes = rep(NA, n_rows),
                                   self_reported_hits = rep(NA, n_rows),
                                   total_hits = rep(NA, n_rows)
                                   )
  rownames(occurrence_results) <- c("Death_register", "Inpatient_data", "Self_report", "Primary_care", "All")
  
  # Algorithms from https://biobank.ndph.ox.ac.uk/ukb/refer.cgi?id=593 (First occurrences PDF)
  
  ## Get read codes
  # ---------------
  codes_tb <- tibble(old_code = provided_codes, value = "yes")
  
  # ICD-9
  icd9_codes_df <- code_lkup_conditions(codes_tb, codes_input_type, "icd9", n_icd_input)
  # ICD-10
  icd10_codes_df <- code_lkup_conditions(codes_tb, codes_input_type, "icd10", n_icd_input) 
  # Read 2 
  read_2_codes_df <- code_lkup_conditions(codes_tb, codes_input_type, "read_2", n_icd_input) 
  # Read 3 
  read_3_codes_df <- code_lkup_conditions(codes_tb, codes_input_type, "read_3", n_icd_input) 
  
  # Save read codes that will be used in the results file (under the relevant database(s))
  # -------------------------------------------------------------------------------------
  # ICD-9 (inpatient data)
  occurrence_results$icd9_codes_short[[2]] <- paste0(icd9_codes_df$new_code, collapse = ", ")
  
  occurrence_results$icd9_codes_full[[2]] <- code_full_description("all_lkps_maps_v3.xlsx", "icd9_lkp",
                                                                   c("ICD9", "DESCRIPTION_ICD9"),
                                                                   icd9_codes_df, n_icd_input) %>%
    tolower()

  # ICD-10 (death register and inpatient data)
  occurrence_results$icd10_codes_short[1:2] <- paste0(icd10_codes_df$new_code, collapse = ", ")
  occurrence_results$icd10_codes_full[1:2] <- code_full_description("all_lkps_maps_v3.xlsx", "icd10_lkp",
                        c("ALT_CODE", "DESCRIPTION"),
                        icd10_codes_df, n_icd_input) 
  
  # Read 2 (primary care data)
  occurrence_results$read_2_codes_short[[4]] <- paste0(read_2_codes_df$new_code, collapse = " ")
  occurrence_results$read_2_codes_full[[4]] <- code_full_description("all_lkps_maps_v3.xlsx", "read_v2_lkp", 
                                                                     c("read_code", "term_description"), 
                                                                     read_2_codes_df) 
  # Read 3 (primary care data)
  occurrence_results$read_3_codes_short[[4]] <- paste0(read_3_codes_df$new_code, collapse = " ")
  occurrence_results$read_3_codes_full[[4]] <- code_full_description("all_lkps_maps_v3.xlsx", "read_ctv3_lkp", 
                                                                     c("read_code", "term_description"), 
                                                                     read_3_codes_df)
  # Search in the databases
  # -----------------------
  # Primary care - uses read codes
  
  gp_clinical_2.df <- pc_mapping_read(read_2_codes_df, "read_2", first_only_input)
  occurrence_results$read_2_hits[[4]] <- length(unique(gp_clinical_2.df$f.eid))
  
  gp_clinical_3.df <- pc_mapping_read(read_3_codes_df, "read_3", first_only_input)
  occurrence_results$read_3_hits[[4]] <- length(unique(gp_clinical_3.df$f.eid))
  
  gp_clinical_4.df <- unique(rbind(gp_clinical_2.df, gp_clinical_3.df))
  occurrence_ids[[4]] <- gp_clinical_4.df
  occurrence_results$total_hits[[4]] <- length(unique(gp_clinical_4.df$f.eid))
  
  # Death register - uses ICD10
  bd_diagnosis <- ons_mapping_icd10(icd10_codes_df,
                                    c("f.40000.", "f.40001.", "f.40002."),
                                    c("f.eid", "I", "date", "diagnosis_date", "death_cause", "death_other_cause"),
                                    first_only_input, n_icd = n_icd_input)
  occurrence_ids[[1]] <- bd_diagnosis
  occurrence_results$icd10_hits[[1]] <- length(unique(bd_diagnosis$f.eid))
  occurrence_results$total_hits[[1]] <- length(unique(bd_diagnosis$f.eid))
  deaths <- bd_diagnosis
  
  # Inpatient data (uses ICD-9 and ICD-10 codes)
  icd10_diagnosis_1 <- pc_mapping_icd(icd10_codes_df, c("f.41202.", "f.41262."), first_only_input, n_icd = n_icd_input) # Primary/main diagnoses only, i.e. codes in diag_icd10 with level=1 (41202	Diagnoses, 41262	Date)
  icd10_diagnosis_2 <- pc_mapping_icd(icd10_codes_df, c("f.41270.", "f.41280."), first_only_input, n_icd = n_icd_input) # Primary and secondary diagnoses, i.e.  including codes in diag_icd10 with level=1, 2 or 3 (41270	Diagnoses, 41280)	
  icd10_diagnosis_1 <- unique(rbind(icd10_diagnosis_1, icd10_diagnosis_2))
  occurrence_results$icd10_hits[[2]] <- length(unique(icd10_diagnosis_1$f.eid))

  icd9_diagnosis_1 <- pc_mapping_icd(icd9_codes_df, c("f.41203.", "f.41263."), first_only_input, n_icd = n_icd_input, code_type = "icd9") # Primary/main diagnoses only, i.e. codes in diag_icd9 with level=1	 (41203	Diagnoses, 41263	Date)
  icd9_diagnosis_2 <- pc_mapping_icd(icd9_codes_df, c("f.41271.", "f.41281."), first_only_input, n_icd = n_icd_input, code_type = "icd9") # Primary and secondary diagnoses, i.e.  including codes in diag_icd9 with level=1, 2  or 3 (41271	Diagnoses, 41281	Date)
  icd9_diagnosis_1 <- unique(rbind(icd9_diagnosis_1, icd9_diagnosis_2))
  
  if (include_icd9 == FALSE) icd9_diagnosis_1 <- tibble(data.frame())
  
  occurrence_results$icd9_hits[[2]] <- ifelse(include_icd9 == TRUE, length(unique(icd9_diagnosis_1$f.eid)), 0)
  
  bd_diagnosis <- unique(rbind(icd10_diagnosis_1, icd9_diagnosis_1))
  occurrence_ids[[2]] <- bd_diagnosis
  occurrence_results$total_hits[[2]] <- length(unique(bd_diagnosis$f.eid))
  
  # Self report
  self_report_codes <-  icd10_codes_df %>%
    rename(icd10 = new_code) %>%
    left_join(fread("Self_report_FO_mappings_Jan2022.tsv")) %>%
    select(enc6, value) %>%
    filter(is.na(enc6) == FALSE)
  bd_diagnosis <- select(bd, contains(c("f.eid", "f.53.", "f.20002.", "f.20008.")))  %>%
    ukb_reshape_long(dict) %>%
    rename (date = `Date of attending assessment centre`,
            enc6 = `Non-cancer illness code, self-reported`,
            condition_year = `Interpolated Year when non-cancer illness first diagnosed`) %>% 
    mutate(condition_year = as.integer(condition_year)) %>%
    filter(is.na(condition_year) == FALSE) %>% 
    filter(condition_year >= 1930) %>% # Removes special values,  no recruited patient was alive before 1930
    left_join(self_report_codes) %>%
    filter(is.na(value) == FALSE) 
  
  if (first_only_input == "TRUE") {
    bd_diagnosis <- bd_diagnosis %>%
      group_by(f.eid) %>%
      summarise(date_first_diagnosed = as.Date(NA),
                year_first_diagnosed = min(condition_year, na.rm = TRUE))
  } else {
    bd_diagnosis <- bd_diagnosis %>%
      select(f.eid, condition_year) %>%
      rename(year_diagnosed = condition_year) %>%
      mutate(date_diagnosed = as.Date(NA)) %>%
      select(f.eid, date_diagnosed, year_diagnosed)
  }

  occurrence_ids[[3]] <- bd_diagnosis
  Codings <- fread("Codings.csv.gz") %>% #load Codings
    subset(Coding == 6) %>%
    select(Value, Meaning) %>%
    rename(enc6 = Value) 
  occurrence_results$self_reported_codes[[3]] <- self_report_codes %>%
    mutate(enc6 = as.character(enc6)) %>%
    left_join(Codings) %>%
    rowwise() %>%
    mutate(codes = paste(enc6, Meaning, sep = ": ")) %>%
    pull(codes) %>%
    paste0(., collapse = ", ")
  occurrence_results$self_reported_hits[[3]] <- length(unique(bd_diagnosis$f.eid))
  occurrence_results$total_hits[[3]] <- length(unique(bd_diagnosis$f.eid))
  
  # All
  occurrence_results$icd9_hits[[5]] <- length(unique(icd9_diagnosis_1$f.eid))
  occurrence_results$icd10_hits[[5]] <- icd10_diagnosis_1 %>%
    rbind(deaths) %>%
    unique() %>%
    nrow()
  
  
  occurrence_results$read_2_hits[[5]] <- length(unique(gp_clinical_2.df$f.eid))
  occurrence_results$read_3_hits[[5]] <- length(unique(gp_clinical_3.df$f.eid))
  occurrence_results$self_reported_hits[[5]] <- length(unique(bd_diagnosis$f.eid))
  occurrence_results$total_hits[[5]] <- length(unique(c(occurrence_ids[[1]]$f.eid, occurrence_ids[[2]]$f.eid, occurrence_ids[[3]]$f.eid, occurrence_ids[[4]]$f.eid)))
  
  if (write_summary == TRUE) write.csv(occurrence_results, paste0("occurrence_results_", substr(paste0(provided_codes, collapse = "_"), 1, 50), ".csv"), row.names = TRUE)
  return(occurrence_ids) 
}

# Obtaining ticagrelor dose quantity 
# ----------------------------------
get_quantity_ticagrelor <- function(quantity_list) {
  quantity_list <- tolower(quantity_list)
  # Deal with ".000"
  if (TRUE %in% str_detect(quantity_list, ".000")) quantity_list[grepl(".000", quantity_list)] <- gsub(".000", "", quantity_list[grepl(".000", quantity_list)])
  # Deal with ".00"
  if (TRUE %in% str_detect(quantity_list, ".00")) quantity_list[grepl(".00", quantity_list)] <- gsub(".00", "", quantity_list[grepl(".00", quantity_list)])
  # Deal with "-"
  quantity_list <- gsub("- pack", "pack", quantity_list)
  if (TRUE %in% str_detect(quantity_list, " - ")) {
    index <- grepl(" - ", quantity_list)
    x <- quantity_list[index]
    y <- c()
    for (i in seq_along(x)) y[i] <- unlist(str_split(x[[i]], " - "))[1]
    quantity_list[index] <- y
  }
  # Deal with "tablets"
  if (TRUE %in% str_detect(quantity_list, " tab")) {
    index <- grepl(" tab", quantity_list) 
    x <- quantity_list[index]
    y <- c()
    for (i in seq_along(x)) y[i] <- unlist(str_split(x[[i]], " tab"))[1]
    quantity_list[index] <- y
  }
  # Deal with "capsules"
  if (TRUE %in% str_detect(quantity_list, " cap")) {
    index <- grepl(" cap", quantity_list) 
    x <- quantity_list[index]
    y <- c()
    for (i in seq_along(x)) y[i] <- unlist(str_split(x[[i]], " cap"))[1]
    quantity_list[index] <- y
  }
  # Deal with packs and x
  quantity_list <- gsub(" pack of ", "*", gsub("packs", "pack", quantity_list))
  quantity_list <- gsub(" x |x", "*", quantity_list)
  # Deal with *
  if (TRUE %in% str_detect(quantity_list, "[0-9]*\\*[0-9]*")) {
    index <- grepl("[0-9]*\\*[0-9]*", quantity_list) 
    x <- quantity_list[index]
    y <- c()
    for (i in seq_along(x)){
      z <- unlist(str_extract_all(str_extract(x[i], "[0-9]*\\*[0-9]*"), "[0-9]*"))
      z <- as.numeric(z[z != ""])
      y[i] <- z[1] * z[2]
    }
    quantity_list[index] <- y
  }
  # Deal with ""
  quantity_list[quantity_list == ""] <- "0"
  # Deal with "-1"
  quantity_list[quantity_list == "-1"] <- "0"
  # Deal with "-1.000"
  quantity_list[quantity_list == "-1.000"] <- "0"
  # Deal with "-1 INVALID"
  quantity_list[quantity_list == "-1 INVALID"] <- "0"

  return(as.numeric(quantity_list))
}