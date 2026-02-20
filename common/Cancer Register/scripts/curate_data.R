library(data.table)

# Make output directory on local cloud workstations
system("mkdir -p cancer_register")

# Check if we have the necessary python libraries installed for dx extract_dataset
if (system("python3 -c 'import pandas'", ignore.stderr=TRUE)) {
  system("pip install pandas", wait=TRUE)
}

# Extract dataset meta-data to local cloud instance, if it hasn't already been downloaded
if (length(list.files(pattern="data_dictionary")) == 0) {
  project_files <- system("dx ls", intern=TRUE)
  dataset_file <- project_files[project_files %like% "^app" & project_files %like% ".dataset$"]
  system(sprintf("dx extract_dataset %s -ddd", dataset_file), wait=TRUE) # Can ignore warning about pandas version
}

# Auto-detect data-dictionary file then load
code_dict_file <- list.files(pattern="coding")
code_dict <- fread(code_dict_file)

# Pull down raw data that has been extracted with Table Exporter
system("mkdir -p raw_data", wait=TRUE)
system("dx download 'common/Cancer Register/raw_data/data.csv' -o raw_data/cancer_register.csv", wait=TRUE)

# Load in raw data and curated field information
header <- fread("raw_data/cancer_register.csv", nrows = 0)
raw <- fread("raw_data/cancer_register.csv", na.strings=c("", "NA"), colClasses=list("character"=names(header)[names(header) %like% "p40013"]))

# Drop rows where the only non-missing data is the participant ID.
raw <- raw[apply(raw, 1, function(row) { sum(!is.na(row)) > 1L })]

# Split out instances for each field so they are rows instead of columns
instances <- setdiff(unique(gsub("^p[0-9]+_", "", names(raw))), "eid")
raw <- rbindlist(fill=TRUE, use.names=TRUE, lapply(instances, function(vr) {
  # Find columns matching this visit repeat pair (e.g. ending in _i0)
  this_cols <- names(raw)[grepl(pattern=paste0(vr, "$"), names(raw))]
  
  # Filter to these columns
  this_raw <- raw[, .SD, .SDcols=c("eid", this_cols)]
  
  # Drop repeat visit pair label from column name
  setnames(this_raw, this_cols, gsub(paste0("_", vr, "$"), "", this_cols))
  
  # Add columns for the record number (instance)
  this_raw[, record_number := as.integer(gsub("^i", "", gsub("_.*", "", vr)))]
  
  # Move to start of data table
  this_raw <- this_raw[,.SD,.SDcols=c("eid", "record_number",  gsub("_.*", "", this_cols))]
  
  # Return
  this_raw
}))

# Drop rows where the only non-missing data is the eid and record_number
raw <- raw[apply(raw, 1, function(row) { sum(!is.na(row)) > 2L })]

# Name fields
setnames(raw, c("p40005", "p40006", "p40008", "p40009", "p40011", "p40012", "p40013", "p40019", "p40021"),
  c("diagnosis_date", "diag_icd10", "age_at_diagnosis", "number_occurrences", "tumor_histology", "tumor_behaviour",
    "diag_icd9", "record_format", "record_origin"))

# Drop number_occurrences -- this column is redundant with the number of records
# per person once you strip out entries missing diagnosis_date
raw[, number_occurrences := NULL]

# Combine ICD codes into a single column
raw[, icd_version := ifelse(is.na(diag_icd10), 9, 10)]
raw[icd_version == 9, diag_icd := diag_icd9]
raw[icd_version == 10, diag_icd := diag_icd10]
raw[, c("diag_icd9", "diag_icd10") := NULL]

# Create ICD-O-3 morphology and behaviour code
raw[!is.na(tumor_histology), icdO3_code := sprintf("M-%s/%s", tumor_histology, tumor_behaviour)]
raw[, c("tumor_histology", "tumor_behaviour") := NULL]

# Replace codings with labels
raw[, record_format := as.character(record_format)]
raw[, record_origin := as.character(record_origin)]
raw[code_dict[coding_name == "data_coding_262"], on = .(record_format=code), record_format := meaning]
raw[code_dict[coding_name == "data_coding_1970"], on = .(record_origin=code), record_origin := meaning]

# Add unique ID for each event
raw[, cr_id := sprintf('%s-%s', eid, record_number)]
stopifnot(raw[,.N,by=cr_id][N > 1, .N] == 0)

# Curate column information
info <- rbind(use.names=FALSE, fill=TRUE,
  data.table(var="cr_id", name="Unique ID for each cancer register diagnosis"),
  data.table("eid", "Application-specific Participant ID"),
  data.table("record_number", "Record number for the cancer diagnosis; usually one per diagnosis"),
  data.table("record_origin", "Broad grouping of registry providers into England & Wales, Scotland, or NCIN (England)", field.id=40021),
  data.table("record_format", "File format of received cancer record", 40019),
  data.table("diagnosis_date", "Date of cancer diagnosis", 40005),
  data.table("age_at_diagnosis", "Calculated by UKB from date of cancer diagnosis and date of birth", 40008),
  data.table("icd_version", "Indicates the ICD version used for the diagnosis code"),
  data.table("diag_icd", "ICD code for the cancer diagnosis (UKB Field IDs 40013 and 40009 for ICD-9 and ICD-10 codes respectively)"),
  data.table("icdO3_code", "Tumor histology and behaviour code (UKB Field IDs 40011 and 40012 combined)") 
)

# Reorder rows and columns
raw <- raw[,.SD,.SDcols=info$var]

# Write out
fwrite(raw, file="cancer_register/cancer_register.csv")
fwrite(info, file="cancer_register/column_information.csv")

# Curate information about ICD-10 codes
icd10 <- code_dict[coding_name == "data_coding_19"]
icd10 <- icd10[, .(code, meaning, display_order, parent_code)]

chapters <- c("Chapter II", "Chapter XV") # Restrict to universe of possible codes
blocks <- icd10[parent_code %in% chapters, code]
two_digit <- icd10[parent_code %in% blocks, code]
three_digit <- icd10[parent_code %in% two_digit, code]
four_digit <- icd10[parent_code %in% three_digit, code]
icd10 <- icd10[code %in% c(chapters, blocks, two_digit, three_digit, four_digit)]

icd10[, appears_in_records := FALSE]
icd10[raw[icd_version == 10], on = .(code=diag_icd), appears_in_records := TRUE]
fwrite(icd10, "cancer_register/icd10_codes.csv")

# Curate column information
info <- rbind(use.names=FALSE,
  data.table(var="code", name="ICD-10 code as it appears in the 'diag_icd' column of cancer_register.csv"),
  data.table("meaning", "Description of ICD-10 code"),
  data.table("display_order", "Numeric ordering of ICD-10 codes for sorting the table"),
  data.table("parent_code", "Parent code of the ICD-10 code"),
  data.table("appears_in_records", "TRUE where this ICD-10 code is found in any cancer register records, FALSE otherwise")
)
fwrite(info, "cancer_register/icd10_codes_column_information.csv")

# Curate information about ICD-9 codes
icd9 <- code_dict[coding_name == "data_coding_87"]
icd9 <- icd9[, .(code, meaning, display_order, parent_code)]

chapters <- c("Chapter II", "Chapter III", "Chapter XI") # Restrict to universe of possible codes
blocks <- icd9[parent_code %in% chapters, code]
three_digit <- icd9[parent_code %in% blocks, code]
four_digit <- icd9[parent_code %in% three_digit, code]
five_digit <- icd9[parent_code %in% four_digit, code]
six_digit <- icd9[parent_code %in% five_digit, code]
icd9 <- icd9[code %in% c(chapters, blocks, three_digit, four_digit, five_digit, six_digit)]

icd9[, appears_in_records := FALSE]
icd9[raw[icd_version == 9], on = .(code=diag_icd), appears_in_records := TRUE]
fwrite(icd9, "cancer_register/icd9_codes.csv")

# Curate column information
info <- rbind(use.names=FALSE,
  data.table(var="code", name="ICD-9 code as it appears in the 'diag_icd' column of cancer_register.csv"),
  data.table("meaning", "Description of ICD-9 code"),
  data.table("display_order", "Numeric ordering of ICD-9 codes for sorting the table"),
  data.table("parent_code", "Parent code of the ICD-9 code"),
  data.table("appears_in_records", "TRUE where this ICD-9 code is found in any cancer register records, FALSE otherwise")
)
fwrite(info, "cancer_register/icd9_codes_column_information.csv")

# Upload to persistent storage
system("dx upload cancer_register/* --destination 'common/Cancer Register/'", wait=TRUE)

# Send raw data to deletion folder to reduce storage costs - this needs to be 
# done in two steps, since we can't move two folders with the same name to the
# same location, so here we rename with a random number and then move to trash
rn <- as.integer(Sys.time())
system(sprintf("dx mv 'common/Cancer Register/raw_data/' 'common/Cancer Register/raw_data_%s'", rn), wait=TRUE)
system(sprintf("dx mv 'common/Cancer Register/raw_data_%s/' trash/", rn), wait=TRUE) 
