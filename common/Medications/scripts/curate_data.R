library(data.table)
library(lubridate)

# Check if we have the necessary python libraries installed for dx extract_dataset
if (system("python3 -c 'import pandas'", ignore.stderr=TRUE)) {
  system("pip install pandas", wait=TRUE)
}

# Extract dataset meta-data to local cloud instance, if it hasn't already been downloaded - we need this for the data codigns
if (length(list.files(pattern="data_dictionary")) == 0) {
  project_files <- system("dx ls", intern=TRUE)
  dataset_file <- project_files[project_files %like% "^app" & project_files %like% ".dataset$"]
  system(sprintf("dx extract_dataset %s -ddd", dataset_file), wait=TRUE) # Can ignore warning about pandas version
}

# Auto-detect coding dictionary file then load
code_dict_file <- list.files(pattern="dataset.codings")
code_dict <- fread(code_dict_file)

# Pull down raw data that has been extracted with Table Exporter
system("mkdir -p raw_data", wait=TRUE)
system("dx download 'common/Medications/raw_data/data.csv' -o raw_data/medication.csv", wait=TRUE)

# Load in raw data
raw <- fread("raw_data/medication.csv", na.strings = c("NA", ""))

# Add in visit index 0 for sex field for downstream code compatibility
setnames(raw, "p31", "p31_i0")

# Add in array index 0 for fields without repeated measures
no_repeats <- names(raw)[names(raw) %like% "i[0-9]$"]
setnames(raw, no_repeats, paste0(no_repeats, "_a0"))

# Split out instance (visit) and array index (repeat measure) fields so they
# are rows instead of columns
visit_repeats <- setdiff(unique(gsub("^p[0-9]+_", "", names(raw))), "eid")
raw <- rbindlist(fill=TRUE, use.names=TRUE, lapply(visit_repeats, function(vr) {
  # Find columns matching this visit repeat pair (e.g. ending in _i0)
  this_cols <- names(raw)[grepl(pattern=paste0(vr, "$"), names(raw))]
  
  # Filter to these columns
  this_raw <- raw[, .SD, .SDcols=c("eid", this_cols)]
  
  # Drop repeat visit pair label from column name
  setnames(this_raw, this_cols, gsub(paste0("_", vr, "$"), "", this_cols))
  
  # Add columns for visit and repeat index
  this_raw[, visit_index := as.integer(gsub("^i", "", gsub("_.*", "", vr)))]
  this_raw[, repeat_index := suppressWarnings(as.integer(gsub("^a", "", gsub(".*_", "", vr))))]
  
  # Move to start of data table
  this_raw <- this_raw[,.SD,.SDcols=c("eid", "visit_index", "repeat_index", gsub("_.*", "", this_cols))]
  
  # Drop instance and array index combinations with all missing data
  # eid, visit_index, and array_index always non-missing
  this_raw <- this_raw[apply(this_raw, 1, function(row) { sum(!is.na(row)) > 3L })]
  
  # Return
  this_raw
}))

# Rename fields - note touchsreen questions for standard medications
# had different fields for males and females, hence the column naming 
# and additional subsequent processing
setnames(raw, c("p31", "p2492", "p137", "p20003", "p6671", "p20199", "p6153", "p6177"),
  c("sex", "other_prescription_medication", "num_current_meds", "medication_code",
   "num_antibiotics_last_3mo", "antibiotic_code", "Female", "Male"))

# Convert codes for sex
raw[, sex := ifelse(sex == 1L, "Male", "Female")]

# Propagate sex information to follow-up assessments
raw[raw[visit_index == 0 & repeat_index == 0], on = .(eid), sex := i.sex]

# Extract standard touchscreen survey medication questions
tchscrn <- raw[repeat_index == 0]
tchscrn[, c("num_current_meds", "medication_code", "num_antibiotics_last_3mo", "antibiotic_code") := NULL]

# Combine male and female specific questions into one column
tchscrn[, medications := Male]
tchscrn[is.na(medications), medications := Female]
tchscrn[, c("Male", "Female") := NULL]

# Convert multiple choice selections to wide-format columns
tchscrn[, cholesterol_medication := fcase(
  is.na(medications), NA, # Touchscreen survey question missing
  medications %in% c(-1, -3), NA, # Do not know (-1) or Prefer not to answer (-3)
  medications %like% "1", TRUE, 
  default = FALSE
)]

tchscrn[, blood_pressure_medication := fcase(
  is.na(medications), NA, # Touchscreen survey question missing
  medications %in% c(-1, -3), NA, # Do not know (-1) or Prefer not to answer (-3)
  medications %like% "2", TRUE, 
  default = FALSE
)]

tchscrn[, insulin_medication := fcase(
  is.na(medications), NA, # Touchscreen survey question missing
  medications %in% c(-1, -3), NA, # Do not know (-1) or Prefer not to answer (-3)
  medications %like% "3", TRUE, 
  default = FALSE
)]

tchscrn[, hormone_replacement_therapy := fcase(
  sex == "Male", FALSE, # Question not asked to male participants
  is.na(medications), NA, # Touchscreen survey question missing
  medications %in% c(-1, -3), NA, # Do not know (-1) or Prefer not to answer (-3)
  medications %like% "4", TRUE, 
  default = FALSE
)]

tchscrn[, oral_contraceptive := fcase(
  sex == "Male", FALSE, # Question not asked to male participants
  is.na(medications), NA, # Touchscreen survey question missing
  medications %in% c(-1, -3), NA, # Do not know (-1) or Prefer not to answer (-3)
  medications %like% "5", TRUE, 
  default = FALSE
)]

# Curate information about other prescription medications
tchscrn[, other_prescription_medication := fcase(
  other_prescription_medication == 1, TRUE,
  other_prescription_medication == 0, FALSE,
  default = NA # includes do not know (-1) and Prefer not to answer (-3)
)]

# Filter and organise columns
tchscrn <- tchscrn[,.(eid, visit_index, cholesterol_medication, blood_pressure_medication,
  insulin_medication, hormone_replacement_therapy, oral_contraceptive, other_prescription_medication)]

# Curate information about touchscreen survey medication columns
info <- rbind(use.names=TRUE, fill=TRUE, 
  data.table(var="eid", name="Application-specific Participant ID"),
  data.table(var="visit_index", name="UK Biobank assessment visit: 0 = Baseline Assessment, 1 = First repeat assessment, 2 = First imaging assessment, 3 = Second imaging assessment"),
  data.table(var="cholesterol_medication", name="Cholesterol medication usage reported during touchscreen survey (UKB Field IDs 6177 for males and 6153 for females)"),
  data.table(var="blood_pressure_medication", name="Blood pressure medication usage reported during touchscreen survey (UKB Field IDs 6177 for males and 6153 for females)"),
  data.table(var="insulin_medication", name="Insulin usage reported during touchscreen survey (UKB Field IDs 6177 for males and 6153 for females)"),
  data.table(var="hormone_replacement_therapy", name="Hormone replacement therapy usage reported during touchscreen survey (UKB Field ID: 6153, asked of females only)"),
  data.table(var="oral_contraceptive", name="Oral contraceptive usage reported during touchscreen survey (UKB Field ID: 6153, asked of females only)"),
  data.table(var="other_prescription_medication", name="Touchscreen question asking if taking any other prescription medication (UKB Field ID: 2492)")
)

# Write out 
fwrite(tchscrn, file="Medication/touchscreen_survey_medications.csv")
fwrite(info, file="Medication/touchscreen_survey_medications_column_info.csv")

# Upload to persistent storage
system("dx upload Medication/touchscreen_survey_medications.csv Medication/touchscreen_survey_medications_column_info.csv --destination common/Medications/", wait=TRUE)

# Curate summary information about the verbal interview with the trained nurse
versum <- raw[repeat_index == 0, .(eid, visit_index, num_current_meds, num_antibiotics_last_3mo)]
info <- rbind(use.names=TRUE, fill=TRUE,
  data.table(var="eid", name="Application-specific Participant ID"),
  data.table(var="visit_index", name="UK Biobank assessment visit: 0 = Baseline Assessment, 1 = First repeat assessment, 2 = First imaging assessment, 3 = Second imaging assessment"),
  data.table(var="num_current_meds", name="Number of prescription medications recorded during verbal interview with trained nurse (UKB Field ID: 137)"),
  data.table(var="num_antibiotics_last_3mo", name="Number of antibiotics taken in the last 3 months recorded during verbal interview with trained nurse (UKB Field ID: 6671)")
)

fwrite(versum, file="Medication/verbal_interview_medications_summary.csv")
fwrite(info, file="Medication/verbal_interview_medications_summary_column_info.csv")

system("dx upload Medication/verbal_interview_medications_summary.csv Medication/verbal_interview_medications_summary_column_info.csv --destination common/Medications/", wait=TRUE)

# Extract information about prescription medications curated by trained nurse
meds <- raw[!is.na(medication_code),.(eid, visit_index, medication_code=as.character(medication_code))]
meds[code_dict[coding_name == "data_coding_4"], on = .(medication_code=code), medication_name := i.meaning]

info <- rbind(use.names=TRUE, fill=TRUE,
  data.table(var="eid", name="Application-specific Participant ID"),
  data.table(var="visit_index", name="UK Biobank assessment visit: 0 = Baseline Assessment, 1 = First repeat assessment, 2 = First imaging assessment, 3 = Second imaging assessment"),
  data.table(var="medication_code", name="Prescription medication code recorded by trained nurse during verbal interview (UKB field ID: 20003)"),
  data.table(var="medication_name", name="Corresponding label for medication code (UKB field ID: 20003)")
)

fwrite(meds, file="Medication/verbal_interview_medications.csv")
fwrite(info, file="Medication/verbal_interview_medications_column_info.csv")

system("dx upload Medication/verbal_interview_medications.csv Medication/verbal_interview_medications_column_info.csv --destination common/Medications/", wait=TRUE)

# Extract information about antibiotics taken in the last 3 months
abio <- raw[!is.na(antibiotic_code),.(eid, visit_index, antibiotic_code=as.character(antibiotic_code))]
abio[code_dict[coding_name == "data_coding_744"], on = .(antibiotic_code=code), antibiotic_name := i.meaning]

info <- rbind(use.names=TRUE, fill=TRUE,
  data.table(var="eid", name="Application-specific Participant ID"),
  data.table(var="visit_index", name="UK Biobank assessment visit: 2 = First imaging assessment, 3 = Second imaging assessment"),
  data.table(var="antibiotic_code", name="Code for antibiotic taken in the past 3 months recorded by trained nurse during verbal interview (UKB field ID: 20199)"),
  data.table(var="antibiotic_name", name="Corresponding label for antibiotic code (UKB field ID: 20199)")
)

fwrite(abio, file="Medication/verbal_interview_antibiotics.csv")
fwrite(info, file="Medication/verbal_interview_antibiotics_column_info.csv")

system("dx upload Medication/verbal_interview_antibiotics.csv Medication/verbal_interview_antibiotics_column_info.csv --destination common/Medications/", wait=TRUE)

# Get list of participants/visit combinations which have any medications
tchscrn_any_meds <- tchscrn[cholesterol_medication | blood_pressure_medication | 
                              insulin_medication | hormone_replacement_therapy | oral_contraceptive | 
                              other_prescription_medication, .(eid, visit_index)]
counts_any_meds <- versum[num_current_meds > 0, .(eid, visit_index)]
verbal_any_meds <- unique(meds[,.(eid, visit_index)])
any_meds <- unique(rbind(tchscrn_any_meds, counts_any_meds, verbal_any_meds))

# Get list of participants/visit combinations which definitively had no medications
tchscrn_no_meds <- tchscrn[!cholesterol_medication & !blood_pressure_medication & 
                             !insulin_medication & !hormone_replacement_therapy & !oral_contraceptive & 
                             !other_prescription_medication, .(eid, visit_index)]
counts_no_meds <- versum[num_current_meds == 0, .(eid, visit_index)]
no_meds <- fintersect(tchscrn_no_meds, counts_no_meds)
no_meds <- no_meds[!meds, on = .(eid, visit_index)]

# Sanity check - tables should be non-overlapping
stopifnot(nrow(fintersect(any_meds, no_meds)) == 0)

# Build full table of medication status
med_status <- unique(rbind(
  tchscrn[,.(eid, visit_index)],
  versum[,.(eid, visit_index)],
  meds[,.(eid, visit_index)],
  abio[,.(eid, visit_index)]
))
med_status[any_meds, on = .(eid, visit_index), any_current_medications := TRUE]
med_status[no_meds, on = .(eid, visit_index), any_current_medications := FALSE]

# Do the same for antibiotic usage in past 3 months
med_status[versum[num_antibiotics_last_3mo == 0], on = .(eid, visit_index), antibiotics_in_past_3_months := FALSE]
med_status[versum[num_antibiotics_last_3mo > 0], on = .(eid, visit_index), antibiotics_in_past_3_months := TRUE]
med_status[abio, on = .(eid, visit_index), antibiotics_in_past_3_months := TRUE]

info <- rbind(use.names=TRUE, fill=TRUE,
  data.table(var="eid", name="Application-specific Participant ID"),
  data.table(var="visit_index", name="UK Biobank assessment visit: 0 = Baseline Assessment, 1 = First repeat assessment, 2 = First imaging assessment, 3 = Second imaging assessment"),
  data.table(var="any_current_medications", name="TRUE/FALSE based on combined records on current medications (UKB Field IDs: 137, 197, 2492, 6177, 6153, or 20003), or NA where absence of medication unable to be definitively determined"),
  data.table(var="antibiotics_in_past_3_months", name="TRUE/FALSE based on combined records of antibiotic usage in the past 3 months (UKB field IDs: 6671 or 20199) or NA where absence of antibiotics unable to be definitively determined")
)

fwrite(med_status, file="Medication/medication_status.csv")
fwrite(info, file="Medication/medication_status_column_info.csv")

system("dx upload Medication/medication_status.csv Medication/medication_status_column_info.csv --destination common/Medications/", wait=TRUE)

# Send raw data to deletion folder to reduce storage costs - this needs to be 
# done in two steps, since we can't move two folders with the same name to the
# same location, so here we rename with a random number and then move to trash
rn <- as.integer(Sys.time())
system(sprintf("dx mv 'common/Medications/raw_data/' 'common/Medications/raw_data_%s'", rn), wait=TRUE)
system(sprintf("dx mv 'common/Medications/raw_data_%s/' trash/", rn), wait=TRUE) 
