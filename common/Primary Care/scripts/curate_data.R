library(data.table)
library(lubridate)

if (suppressMessages(suppressWarnings(!require(openxlsx)))) {
  # installing missing R packages
  install.packages("openxlsx")
  suppressMessages(library("openxlsx"))
}

# WARNING: You are going to need a *beefy* instance to run this script!
# mem1_ssd1_v2_x72

# create output directory on local cloud workstation
system("mkdir -p primary_care")

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

# Auto-detect coding-dictionary file then load
code_dict_file <- list.files(pattern="coding")
code_dict <- fread(code_dict_file)

# Download raw data extracted by Table Exporter
system("mkdir -p raw_data")
system("dx download 'common/Primary Care/raw_data/*' -o raw_data/")

# Download Read code mappings and documentation
system("wget -P primary_care/ https://biobank.ndph.ox.ac.uk/showcase/ukb/auxdata/primarycare_codings.zip")
system("unzip primary_care/primarycare_codings.zip -d primary_care")
system("mv primary_care/all_lkps_maps_v4.xlsx primary_care/read_code_lookup_maps.xlsx")
system("mv primary_care/all_lkps_maps_variable_names_definitions_v4.pdf primary_care/read_code_lookup_maps_documentation.pdf")
system("rm primary_care/primarycare_codings.zip")

#################################
# Curate GP registration records
#################################
gp_registrations <- fread("raw_data/gp_registrations.csv")

gp_registrations[, data_provider := as.character(data_provider)]
gp_registrations[code_dict[coding_name == "data_coding_626"], on = .(data_provider=code), data_provider := meaning]

# Replace special nonsense date codes with NAs
gp_registrations[, reg_date := as.character(reg_date)]
gp_registrations[code_dict[coding_name == "data_coding_819"], on = .(reg_date=code), reg_date := NA] 
setnames(gp_registrations, "reg_date", "registration_date")

gp_registrations[, deduct_date := as.character(deduct_date)]
gp_registrations[code_dict[coding_name == "data_coding_819"], on = .(deduct_date=code), deduct_date := NA]
setnames(gp_registrations, "deduct_date", "deregistration_date")

fwrite(gp_registrations, "primary_care/gp_registrations.csv")

# Curate column information
info <- rbind(use.names=FALSE,
  data.table(var="eid", name="Application-specific Participant ID"),
  data.table("data_provider", "Data provider: Whether the record originates from England (Vision), Scotland, England (TPP) or Wales"),
  data.table("registration_date", "Date participant registered at a specific general practice"),
  data.table("deregistration_date", "Date particpant was de-registered from that general practice, if applicable")
)
fwrite(info, "primary_care/gp_registrations_column_information.csv")

#################################
# Curate clinical records
#################################
clinical_records <- fread("raw_data/gp_clinical.csv", na.strings=c("", "NA"))

clinical_records[, data_provider := as.character(data_provider)]
clinical_records[code_dict[coding_name == "data_coding_626"], on = .(data_provider=code), data_provider := meaning]

# Replace special nonsense date codes with NAs
clinical_records[, event_dt := as.character(event_dt)]
clinical_records[code_dict[coding_name == "data_coding_819"], on = .(event_dt=code), event_dt := NA] 
setnames(clinical_records, "event_dt", "event_date")

# Combine read code columns
clinical_records[, read_version := ifelse(!is.na(read_2), 2, 3)]
clinical_records[read_version == 2, read_code := read_2]
clinical_records[read_version == 3, read_code := read_3]
clinical_records[, c("read_2", "read_3") := NULL]

# Add in mappings to ICD-10 codes
read2_to_icd10 <- as.data.table(read.xlsx("primary_care/read_code_lookup_maps.xlsx", "read_v2_icd10"))
read3_to_icd10 <- as.data.table(read.xlsx("primary_care/read_code_lookup_maps.xlsx", "read_ctv3_icd10"))
read_to_icd10 <- rbind(idcol="read_version", "2"=read2_to_icd10[,.(read_code, icd10_code)], "3"=read3_to_icd10[,.(read_code, icd10_code)])
read_to_icd10[, read_version := as.numeric(read_version)]
clinical_records[read_to_icd10, on = .(read_version, read_code), icd10_code := i.icd10_code]

# Add in mappings to OPCS-4 codes
read2_to_opcs4 <- as.data.table(read.xlsx("primary_care/read_code_lookup_maps.xlsx", "read_v2_opcs4"))
read3_to_opcs4 <- as.data.table(read.xlsx("primary_care/read_code_lookup_maps.xlsx", "read_ctv3_opcs4"))
read_to_opcs4 <- rbind(idcol="read_version", "2"=read2_to_opcs4[,.(read_code, opcs4_code=opcs_4.2_code)], "3"=read3_to_opcs4[,.(read_code, opcs4_code)])
read_to_opcs4[, read_version := as.numeric(read_version)]
clinical_records[read_to_opcs4, on = .(read_version, read_code), opcs4_code := i.opcs4_code]

# Combine value columns into 1 column:
clinical_records <- melt(clinical_records, measure.vars=c("value1", "value2", "value3"), variable.name="value_number")
clinical_records[, value_number := as.numeric(gsub("value", "", value_number))]
clinical_records <- clinical_records[value_number == 1 | !is.na(value)] # Not all records have values associated
clinical_records[is.na(value), value_number := NA]

# Add unique ID for records
clinical_records[, gp_id := sprintf("%s-%s", eid, 1:.N), by=eid]

# Impute missing event dates as the midpoint in the records for each participant
if (!("IDate" %in% class(clinical_records$event_date))) clinical_records[, event_date := as.IDate(event_date)]
clinical_records[, imputed_date := ifelse(is.na(event_date), NA, FALSE)]
record_range <- clinical_records[!is.na(event_date), .(min=min(event_date), max=max(event_date)), by=eid]
record_range[, midpoint := as.IDate(date_decimal((decimal_date(max) - decimal_date(min))/2 + decimal_date(min)))]
missing_dates <- clinical_records[is.na(event_date)]
missing_dates[record_range, on = .(eid), event_date := i.midpoint]
clinical_records[missing_dates[!is.na(event_date)], on = .(gp_id), 
  c("event_date", "imputed_date") := .(i.event_date, TRUE)]
  
# Curate column information
info <- rbind(use.names=FALSE,
  data.table(var="gp_id", name="Unique ID for this clinical record"),
  data.table("eid", "Application-specific Participant ID"),
  data.table("data_provider", "Data provider: Whether the record originates from England (Vision), Scotland, England (TPP) or Wales"),
  data.table("event_date", "Date of the clinical code was recorded or imputed if missing"),
  data.table("imputed_date", "TRUE where 'event_date' was imputed"),
  data.table("read_version", "Read code version used for the corresponding read code; either 2 or 3 depending on the data provider"),
  data.table("read_code", "Read code for the event"),
  data.table("value_number", "Index for any values (up to 3) recorded alongside the read code"),
  data.table("value", "value associated with the read code - if any"),
  data.table("icd10_code", "ICD-10 code mapped to read code by UK Biobank; see read_code_lookup_maps.xlsx"),
  data.table("opcs4_code", "OPCS-4 code mapped to read code by UK Biobank; see read_code_lookup_maps.xlsx")
)

clinical_records <- clinical_records[,.SD,.SDcols=info$var]
fwrite(clinical_records, "primary_care/gp_clinical_records.csv")
fwrite(info, "primary_care/gp_clinical_records_column_information.csv")

#################################
# Curate prescription records
#################################
prescriptions <- fread("raw_data/gp_scripts.csv", na.strings=c("", "NA"))

prescriptions[, data_provider := as.character(data_provider)]
prescriptions[code_dict[coding_name == "data_coding_626"], on = .(data_provider=code), data_provider := meaning]

# Replace special nonsense date codes with NAs
prescriptions[, issue_date := as.character(issue_date)]
prescriptions[code_dict[coding_name == "data_coding_819"], on = .(issue_date=code), issue_date := NA] 

fwrite(prescriptions, "primary_care/gp_prescriptions.csv")

# Curate column information
info <- rbind(use.names=FALSE,
  data.table(var="eid", name="Application-specific Participant ID"),
  data.table("data_provider", "Data provider: Whether the record originates from England (Vision), Scotland, England (TPP) or Wales"),
  data.table("issue_date", "Date prescription was issued"),
  data.table("read_2", "Read code version 2 for the drug; drugs may be coded with read_2, bnf_code, and/or dmd_code"),
  data.table("bnf_code", "British National Formulary code for the drug; drugs may be coded with read_2, bnf_code, and/or dmd_code"),
  data.table("dmd_code", "Dictionary of Medicines and Devices code for the drug; drugs may be coded with read_2, bnf_code, and/or dmd_code"),
  data.table("drug_name", "Name of the prescribed drug"),
  data.table("quantity", "Amount and/or dosage prescribed")
)
fwrite(info, "primary_care/gp_prescriptions_column_information.csv")

#########################################################
# Curate information on Read-2 codes present in the data
#########################################################
read2 <- as.data.table(read.xlsx("primary_care/read_code_lookup_maps.xlsx", sheet="read_v2_lkp"))
read2 <- read2[!is.na(term_code),.(read_code, description=term_description)]
read2[read2_to_icd10, on = .(read_code), icd10 := icd10_code]
read2_to_opcs4 <- as.data.table(read.xlsx("primary_care/read_code_lookup_maps.xlsx", sheet="read_v2_opcs4"))
read2[read2_to_opcs4, on = .(read_code), opcs4 := opcs_4.2_code]
read2[, in_clinical_records := FALSE]
read2[read_code %in% clinical_records[read_version == 2, unique(read_code)], in_clinical_records := TRUE]
read2[(in_clinical_records), has_value := FALSE]
read2[unique(clinical_records[!is.na(value) & read_version == 2, .(read_code)]), on = .(read_code), has_value := TRUE]
fwrite(read2, "primary_care/read2_codes.csv")

info <- rbind(use.names=FALSE,
  data.table(var="read_code", name="Read-2 code"),
  data.table("description", "Description of Read code from read_v2_lkp sheet in read_code_lookup_maps.xlsx"),
  data.table("icd10", "ICD-10 code mapped to Read code by UK Biobank in read_v2_icd10 sheet in read_code_lookup_maps.xlsx"),
  data.table("opcs4", "OPCS-4 code mapped to Read code by UK Biobank in read_v2_opcs4 sheet in read_code_lookup_maps.xlsx"),
  data.table("in_clinical_records", "TRUE/FALSE indicating whether there are any instances of the Read code in gp_clincal_records.csv"),
  data.table("has_value", "TRUE/FALSE indicating whether any records with the Read code in gp_clinical_records.csv have an associated value")
)
fwrite(info, "primary_care/read2_codes_column_information.csv")

#########################################################
# Curate information on READ-3 codes present in the data
#########################################################
read3 <- as.data.table(read.xlsx("primary_care/read_code_lookup_maps.xlsx", sheet="read_ctv3_lkp"))
read3 <- read3[!is.na(status),.(read_code, description=term_description)]
read3[read3_to_icd10, on = .(read_code), icd10 := icd10_code]
read3_to_opcs4 <- as.data.table(read.xlsx("primary_care/read_code_lookup_maps.xlsx", sheet="read_ctv3_opcs4"))
read3[read3_to_opcs4, on = .(read_code), opcs4 := opcs4_code]
read3[, in_clinical_records := FALSE]
read3[read_code %in% clinical_records[read_version == 3, unique(read_code)], in_clinical_records := TRUE]
read3[(in_clinical_records), has_value := FALSE]
read3[unique(clinical_records[!is.na(value) & read_version == 3, .(read_code)]), on = .(read_code), has_value := TRUE]
fwrite(read3, "primary_care/read3_codes.csv")

info <- rbind(use.names=FALSE,
  data.table(var="read_code", name="Read-3 code"),
  data.table("description", "Description of Read code from read_ctv3_lkp sheet in read_code_lookup_maps.xlsx"),
  data.table("icd10", "ICD-10 code mapped to Read code by UK Biobank in read_ctv3_icd10 sheet in read_code_lookup_maps.xlsx"),
  data.table("opcs4", "OPCS-4 code mapped to Read code by UK Biobank in read_ctv3_opcs4 sheet in read_code_lookup_maps.xlsx"),
  data.table("in_clinical_records", "TRUE/FALSE indicating whether there are any instances of the Read code in gp_clincal_records.csv"),
  data.table("has_value", "TRUE/FALSE indicating whether any records with the Read code in gp_clinical_records.csv have an associated value")
)
fwrite(info, "primary_care/read3_codes_column_information.csv")

###############################
# Upload to persistent storage
###############################
system("dx upload primary_care/* --destination 'common/Primary Care/'", wait=TRUE)

###############################
# Cleanup
###############################
# Send raw data to deletion folder to reduce storage costs - this needs to be 
# done in two steps, since we can't move two folders with the same name to the
# same location, so here we rename with a random number and then move to trash
rn <- as.integer(Sys.time())
system(sprintf("dx mv 'common/Primary Care/raw_data/' 'common/Primary Care/raw_data_%s'", rn), wait=TRUE)
system(sprintf("dx mv 'common/Primary Care/raw_data_%s/' trash/", rn), wait=TRUE) 
