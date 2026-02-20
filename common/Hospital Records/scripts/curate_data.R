library(data.table)

# Make output directory on local cloud workstation
system("mkdir -p hospital_records/")

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

# Load in ICD codings map
code_dict_file <- list.files(pattern="codings")
code_dict <- fread(code_dict_file)

# Pull down information that has been extracted with Table Exporter
system("mkdir -p raw_data", wait=TRUE)
system("dx download 'common/Hospital Records/raw_data/*' -o raw_data/", wait=TRUE)

# Process main hospital records table
hes <- fread("raw_data/hesin.csv")

setnames(hes, "ins_index", "record_number")

setnames(hes, "dsource", "hospital_nation")
hes[, hospital_nation := fcase(
  hospital_nation == "HES", "England",
  hospital_nation == "PEDW", "Wales",
  hospital_nation == "SMR", "Scotland"
)]

hes[, source := as.character(source)]
hes[code_dict[coding_name == "data_coding_263"], on = .(source=code), source := meaning]
setnames(hes, "source", "record_source")

setnames(hes, "epistart", "episode_start")
setnames(hes, "epiend", "episode_end")
setnames(hes, "epidur", "episode_duration")
setnames(hes, "bedyear", "hes_year_duration")

hes[, epistat := as.character(epistat)]
hes[code_dict[coding_name == "data_coding_237"], on = .(epistat=code), epistat := meaning]
hes[, epistat := ifelse(epistat == "Finished", FALSE, TRUE)]
setnames(hes, "epistat", "unfinished_episode")

hes[, epitype := as.character(epitype)]
hes[code_dict[coding_name == "data_coding_238"], on = .(epitype=code), epitype := meaning]
setnames(hes, "epitype", "episode_type")

hes[, epiorder := as.character(epiorder)]
hes[code_dict[coding_name == "data_coding_236"], on = .(epiorder=code), epiorder := NA] # all flavours of missing data
hes[, epiorder := as.integer(epiorder)]
setnames(hes, "epiorder", "episode_number")

hes[, spell_seq := spell_seq + 1]
setnames(hes, "spell_seq", "episode_number_chrono")

hes[, first_episode := ifelse(spelbgin > 0, TRUE, FALSE)]
hes[, started_previous_year := fcase(
  spelbgin == 1, TRUE,
  spelbgin == 2, FALSE,
  default = NA
)]
setnames(hes, "spelbgin", "first_episode")

hes[,spelend := ifelse(spelend == "Y", TRUE, FALSE)]
setnames(hes, "spelend", "last_episode")

setnames(hes, "speldur", "spell_duration")

hes[code_dict[coding_name == "data_coding_239"], on = .(pctcode=code), pctcode := meaning]
hes[pctcode == "", pctcode := NA]
setnames(hes, "pctcode", "administrative_body")

hes[code_dict[coding_name == "data_coding_239"], on = .(gpprpct=code), gpprpct := meaning] 
hes[gpprpct == "", gpprpct := NA]
setnames(hes, "gpprpct", "gp_authority")

hes[, category := as.character(category)]
hes[code_dict[coding_name == "data_coding_201"], on = .(category=code), category := meaning]
hes[category %in% c("Not applicable", "Not known: a validation error"), category := NA]
setnames(hes, "category", "legal_status")

setnames(hes, "admidate", "admission_date")
setnames(hes, "elecdate", "decision_date")
setnames(hes, "elecdur", "waiting_time")

hes[, admimeth_uni := as.character(admimeth_uni)]
hes[code_dict[coding_name == "data_coding_264"], on = .(admimeth_uni=code), admimeth_uni := meaning]
hes[admimeth_uni %in% c("Not applicable", "Not known"), admimeth_uni := NA]
setnames(hes, "admimeth_uni", "admission_method")

hes[, admisorc_uni := as.character(admisorc_uni)]
hes[code_dict[coding_name == "data_coding_265"], on = .(admisorc_uni=code), admisorc_uni := meaning]
hes[admisorc_uni %in% c("Not applicable", "Not known"), admisorc_uni := NA]
setnames(hes, "admisorc_uni", "admission_source")

hes[, firstreg := fcase(
  firstreg == 0, TRUE,
  firstreg == 1, FALSE,
  default=NA # 8, 9, or NA
)]
setnames(hes, "firstreg", "first_regular")

hes[, classpat_uni := as.character(classpat_uni)]
hes[code_dict[coding_name == "data_coding_266"], on = .(classpat_uni=code), classpat_uni := meaning]
hes[classpat_uni == "Not applicable", classpat_uni := NA]
setnames(hes, "classpat_uni", "patient_classification")

hes[, intmanag_uni := as.character(intmanag_uni)]
hes[code_dict[coding_name == "data_coding_271"], on = .(intmanag_uni=code), intmanag_uni := meaning]
hes[intmanag_uni %in% c("Not applicable", "Not known"), intmanag_uni := NA]
setnames(hes, "intmanag_uni", "intended_management")

hes[, mainspef_uni := as.character(mainspef_uni)]
hes[code_dict[coding_name == "data_coding_270"], on = .(mainspef_uni=code), mainspef_uni := meaning]
setnames(hes, "mainspef_uni", "consultant_speciality")

hes[, tretspef_uni := as.character(tretspef_uni)]
hes[code_dict[coding_name == "data_coding_269"], on = .(tretspef_uni=code), tretspef_uni := meaning]
setnames(hes, "tretspef_uni", "treatment_speciality")

hes[, operstat := fcase(
  operstat == 1, TRUE,
  operstat == 8, FALSE,
  default = NA
)]
setnames(hes, "operstat", "operation_performed")

setnames(hes, "disdate", "discharge_date")

hes[, dismeth_uni := as.character(dismeth_uni)]
hes[code_dict[coding_name == "data_coding_268"], on = .(dismeth_uni=code), dismeth_uni := meaning]
hes[dismeth_uni %in% c("Not applicable", "Not known"), dismeth_uni := NA]
setnames(hes, "dismeth_uni", "discharge_method")

hes[, disdest_uni := as.character(disdest_uni)]
hes[code_dict[coding_name == "data_coding_267"], on = .(disdest_uni=code), disdest_uni := meaning]
hes[disdest_uni %in% c("Not applicable", "Not applicable: patient died", "Not known"), disdest_uni := NA]
setnames(hes, "disdest_uni", "discharge_destination")

hes[, carersi := fcase(
  carersi == 1, TRUE,
  carersi == 2, FALSE,
  default = NA
)]
setnames(hes, "carersi", "carer_support")

# Curate column information files
info <- rbind(use.names=FALSE,
  data.table(var="dnx_hesin_id", name="Unique identifier for this record"),
  data.table("eid", "Application-specific Participant ID"),
  data.table("record_number", "Index for this record for this person, n.b. number does not necessarily reflect chronological ordering of events"),
  data.table("hospital_nation", "England, Wales, and Scotland use different systems (HES, PEDW, and SMR) with differences in linkage available"),
  data.table("record_source", "Specific source of record; different sources use different record formats"),
  data.table("episode_start", "Episode start date; an episode is a continuous period of care administered under one consultant within one hospital provider"),
  data.table("episode_end", "Episode end date; an episode is a continuous period of care administered under one consultant within one hospital provider"),
  data.table("episode_duration", "Duration of episode (days) if both episode start and end dates are given"),
  data.table("hes_year_duration", "For records in England, duration of episode (days) within the HES data year"),
  data.table("unfinished_episode", "For records in England, unfinished episodes within a HES data year are replaced with the finished episode in the subsequent HES year"),
  data.table("episode_type", "Type of episode; e.g. general admissions, psychiatric, or maternity"),
  data.table("spell_index", "For records from England or Wales, episodes are grouped into spells covering a total continuous stay in a hospital from admission to discharge, n.b. spells are not indexed in chronological order"),
  data.table("episode_number", "The position of this episode within a spell"),
  data.table("episode_number_chrono", "Chronological position of this epsidoe within a spell (spearman correlation with episode_number: 0.91)"),
  data.table("first_episode", "For records from England, indicates whether this episode is the first in the spell"),
  data.table("started_previous_year", "If this is the first episode in the spell, did the spell start in the previous (HES?) year?"),
  data.table("last_episode", "For records from England, indicates whether this episode is the last in the spell"),
  data.table("spell_duration", "Duration of the spell (days)"),
  data.table("first_regular", "Is this the first regular day or night of admission to hospital in a series of stays?"),
  data.table("admission_date", "Date of admission to the hospital, not available for Hospitals in Wales or where the record source was 'UKB Scottish inpatient (SMR01) v4'"),
  data.table("decision_date", "Date of decision to admit to hospital, not available for Hospitals in Wales or where the record source was 'UKB Scottish inpatient (SMR01) v4'"),
  data.table("waiting_time", "Waiting time for hospital admission (days), not available for Hospitals in Wales or where the record source was 'UKB Scottish inpatient (SMR01) v4'"),
  data.table("administrative_body", "Administrative body responsible for the patient during this episode"),
  data.table("gp_authority", "Administrative body where the patient's usual general practice is registered"),
  data.table("legal_status", "Administrative and legal status of patient, e.g. if they have been detained under the Mental Health Act"),
  data.table("admission_method", "Method of admission to hospital (harmonized across all sources by UK Biobank)"),
  data.table("admission_source", "Source of admission to hospital (harmonized across all sources by UK Biobank)"),
  data.table("patient_classification", "Patient classification on admission (harmonized across all sources by UK Biobank)"),
  data.table("intended_management", "Intended management of patient (harmonized across all sources by UK Biobank)"),
  data.table("consultant_speciality", "Main speciality of consultant (harmonized across all sources by UK Biobank)"),
  data.table("treatment_speciality", "Treatment speciality of consultant (harmonized across all sources by UK Biobank)"),
  data.table("operation_performed", "Was an operation performed as part of this hospital episode?"),
  data.table("discharge_date", "Date of discharge from the hospital"),
  data.table("discharge_method", "Method of discharge from hospital (harmonized across all sources by UK Biobank)"),
  data.table("discharge_destination", "Destination on discharge from hospital (harmonized across all sources by UK Biobank)"),
  data.table("carer_support", "Is a carer available for the patient at their usual place of residence?")
)
fwrite(info, "hospital_records/hospital_records_column_information.csv")

# Order hospital records to match info column order
hes <- hes[,.SD, .SDcols=info$var]
fwrite(hes, "hospital_records/hospital_records.csv")

# Load in and process diagnosis codes
diag <- fread("raw_data/hesin_diag.csv", na.strings=c("", "NA"), colClasses = c("diag_icd9"="character"))
setnames(diag, "ins_index", "hospital_record")
setnames(diag, "arr_index", "diagnosis_record")
diag[, level := fcase(
  level == 1, "primary",
  level == 2, "secondary",
  level == 3, "external"
)]
setnames(diag, "level", "cause_type")

diag[, icd_version := ifelse(!is.na(diag_icd9), 9, 10)]
diag[icd_version == 9, diag_icd := diag_icd9]
diag[icd_version == 10, diag_icd := diag_icd10]
diag[icd_version == 9, icd_note := as.character(diag_icd9_nb)] # always NA
diag[icd_version == 10, icd_note := diag_icd10_nb]
diag[, c("diag_icd9", "diag_icd10", "diag_icd9_nb", "diag_icd10_nb") := NULL]
diag[icd_note == 1, icd_note := "* (manifestation)"]
diag[icd_note == 2, icd_note := "† (cause)"]

# And in key information from the main hes table to reduce the need for computationally
# expensive joins in downstream analyses e.g. common/Endpoints/curate_endpoint.R
diag[hes, on = .(dnx_hesin_id), hospital_nation := i.hospital_nation]

# Add in minimum and maximum dates for each hospitalisation associated with a
# given diagnosis (minimum date used for incident events, max for prevalent)
#
# If episode start date is missing (n=57,078 of 4,232,145 events), use admission date instead
# If episode end date is missing (n=57,233 of 4,232,145 events), use discharge data instead if available (in all but 9 cases these are the same as the episode_start date)
# In remaining cases (n=182 of 4,232,145 events), use the episode start date as the end date
# 
diag[hes, on = .(dnx_hesin_id), event_start := episode_start]
no_start <- diag[is.na(event_start), unique(dnx_hesin_id)]
diag[hes[dnx_hesin_id %in% no_start], on = .(dnx_hesin_id), event_start := admission_date] 

diag[hes, on = .(dnx_hesin_id), event_end := episode_end]
no_end <- diag[is.na(event_end), unique(dnx_hesin_id)]
diag[hes[dnx_hesin_id %in% no_start], on = .(dnx_hesin_id), event_end := discharge_date]
diag[is.na(event_end), event_end := event_start]

fwrite(diag, "hospital_records/diagnoses.csv")

# Curate column information files
info <- rbind(use.names=FALSE,
  data.table(var="dnx_hesin_diag_id", name="Unique identifier for this record"),
  data.table("dnx_hesin_id", "Unique identifier for entry in hospital_records.csv this diagnosis code is associated with"),
  data.table("eid", "Application-specific Participant ID"),
  data.table("hospital_record", "Index for the record in hospital_records.csv for this participant"),
  data.table("diagnosis_record", "Index for this record for this diagnosis and hospital record for this participant"),
  data.table("cause_type", "Indicates whether this diagnosis is a primary (main reason), secondary (contributing factor), or external cause (e.g. injury) for hospitalisation"),
  data.table("icd_version", "Indicates the ICD version used for the diagnosis code"),
  data.table("icd_code", "ICD code for the diagnosis"),
  data.table("icd_note", "In some cases ICD-10 codes were accompanied by an addendum †, indicating the coded diagnosis is part of the underlying cause, or * indicating the the coded diagnosis is a manifestation of the disease"),
  data.table("hospital_nation", "England, Wales, and Scotland use different systems (HES, PEDW, and SMR) with differences in linkage available"),
  data.table("event_start", "Start date of the hospital episode the diagnoses took place in, or date of admission to hospital if episode information was missing"),
  data.table("event_end", "End date of the hospital episode the diagnoses took place in, or date of discharge from hospital if episode information was missing")
)
fwrite(info, "hospital_records/diagnoses_column_information.csv")

# Load in and process operation codes
oper <- fread("raw_data/hesin_oper.csv", na.strings=c("", "NA"), colClasses = c("oper3"="character"))
setnames(oper, "ins_index", "hospital_record")
setnames(oper, "arr_index", "operation_record")
oper[, procedure_type := ifelse(level == 1, "main", "secondary")]
setnames(oper, "opdate", "operation_date")
setnames(oper, "preopdur", "pre_operation_stay")
setnames(oper, "posopdur", "post_operation_stay")

oper[, opcs_version := ifelse(!is.na(oper3), 3, 4)]
oper[opcs_version == 3, opcs_code := as.character(oper3)]
oper[opcs_version == 4, opcs_code:= oper4]
oper[opcs_version == 3, opcs_note := as.character(oper3_nb)] # always NA
oper[opcs_version == 4, opcs_note := oper4_nb] # also always NA
oper[, c("oper3", "oper4", "oper3_nb", "oper4_nb") := NULL]

# And in key information from the main hes table to reduce the need for computationally
# expensive joins in downstream analyses e.g. common/Endpoints/curate_endpoint.R
oper[hes, on = .(dnx_hesin_id), hospital_nation := i.hospital_nation]

# Do the same for operations and procedures
#
# If the operation date is missing (n=1,019,419 of 4,071,285 events), use the episode start date instead
# If the operation date is still missing (n=15,226 of 4,071,285 events), use admission date instead
# Note 20 cases remaining with NA operation_date, which are all X998 - No procedure performed
#
# Unlike 'diag' there is no separate end date here; as episode_end and discharge_date
# in 'hes' often include post-op recovery time, which we don't want to use as the end 
# date if calculating retrospective follow-up
#
oper[, imputed_operation_date := operation_date]
no_start <- oper[is.na(imputed_operation_date), unique(dnx_hesin_id)]
oper[hes[dnx_hesin_id %in% no_start], on = .(dnx_hesin_id), imputed_operation_date := episode_start] 
no_start <- oper[is.na(imputed_operation_date), unique(dnx_hesin_id)]
oper[hes[dnx_hesin_id %in% no_start], on = .(dnx_hesin_id), imputed_operation_date := admission_date] 

# Curate column information files
info <- rbind(use.names=FALSE,
  data.table(var="dnx_hesin_oper_id", name="Unique identifier for this record"),
  data.table("dnx_hesin_id", "Unique identifier for entry in hospital_records.csv this operation code is associated with"),
  data.table("eid", "Application-specific Participant ID"),
  data.table("hospital_record", "Index for the record in hospital_records.csv for this participant"),
  data.table("operation_record", "Index for this record for this operation and hospital record for this participant"),
  data.table("procedure_type", "Indicates whether the associated code is for the main procedure or any secondary procedures carried out during this operation"),
  data.table("operation_date", "Date of operation"),
  data.table("imputed_operation_date", "Where operation_date is missing, imputed as episode start date or date of admission to hospital"),
  data.table("opcs_version", "Indicates the OPCS version used for the operation code (currently all OPCS version 4)"),
  data.table("opcs_code", "OPCS code for the operation"),
  data.table("opcs_note", "Any accompanying addendum wit the OPCS code"),
  data.table("pre_operation_stay", "Duration (days) of pre-operative stay in hospital"),
  data.table("post_operation_stay", "Duration (days) of post-operative stay in hospital"),
  data.table("hospital_nation", "England, Wales, and Scotland use different systems (HES, PEDW, and SMR) with differences in linkage available")
)
fwrite(info, "hospital_records/operations_column_information.csv")

oper <- oper[,.SD,.SDcols=info$var]
fwrite(oper, "hospital_records/operations.csv")

# Curate information about ICD-10 codes
icd10 <- code_dict[coding_name == "data_coding_19"]
icd10 <- icd10[, .(code, meaning, display_order, parent_code)]
icd10[, appears_in_records := FALSE]
icd10[diag[icd_version == 10], on = .(code=diag_icd), appears_in_records := TRUE]
fwrite(icd10, "hospital_records/icd10_codes.csv")

# Curate column information
info <- rbind(use.names=FALSE,
  data.table(var="code", name="ICD-10 code as it appears in the 'diag_icd' column of diagnoses.csv"),
  data.table("meaning", "Description of ICD-10 code"),
  data.table("display_order", "Numeric ordering of ICD-10 codes for sorting the table"),
  data.table("parent_code", "Parent code of the ICD-10 code"),
  data.table("appears_in_records", "TRUE where this ICD-10 code is found in any hospital records, FALSE otherwise")
)
fwrite(info, "hospital_records/icd10_codes_column_information.csv")

# Curate information about ICD-9 codes
icd9 <- code_dict[coding_name == "data_coding_87"]
icd9 <- icd9[, .(code, meaning, display_order, parent_code)]
icd9[, appears_in_records := FALSE]
icd9[diag[icd_version == 9], on = .(code=diag_icd), appears_in_records := TRUE]
fwrite(icd9, "hospital_records/icd9_codes.csv")

# Curate column information
info <- rbind(use.names=FALSE,
  data.table(var="code", name="ICD-9 code as it appears in the 'diag_icd' column of diagnoses.csv"),
  data.table("meaning", "Description of ICD-9 code"),
  data.table("display_order", "Numeric ordering of ICD-9 codes for sorting the table"),
  data.table("parent_code", "Parent code of the ICD-9 code"),
  data.table("appears_in_records", "TRUE where this ICD-9 code is found in any hospital records, FALSE otherwise")
)
fwrite(info, "hospital_records/icd9_codes_column_information.csv")

# Curate information about OPCS-4 codes
opcs4 <- code_dict[coding_name == "data_coding_240"]
opcs4 <- opcs4[, .(code, meaning, display_order, parent_code)]
opcs4[, appears_in_records := FALSE]
opcs4[oper[opcs_version == 4], on = .(code=opcs_code), appears_in_records := TRUE]
fwrite(opcs4, "hospital_records/opcs4_codes.csv")

# Curate column information
info <- rbind(use.names=FALSE,
  data.table(var="code", name="OPCS-4 code as it appears in the 'opcs_code' column of operations.csv"),
  data.table("meaning", "Description of OPCS-4 code"),
  data.table("display_order", "Numeric ordering of OPCS-4 codes for sorting the table"),
  data.table("parent_code", "Parent code of the OPCS-4  code"),
  data.table("appears_in_records", "TRUE where this OPCS-4 code is found in any hospital records, FALSE otherwise")
)
fwrite(info, "hospital_records/opcs4_codes_column_information.csv")

# Curate information about OPCS-3 codes
opcs3 <- code_dict[coding_name == "data_coding_259"]
opcs3 <- opcs3[, .(code, meaning, display_order, parent_code)]
opcs3[, appears_in_records := FALSE]
opcs3[oper[opcs_version == 3], on = .(code=opcs_code), appears_in_records := TRUE]
fwrite(opcs3, "hospital_records/opcs3_codes.csv")

# Curate column information
info <- rbind(use.names=FALSE,
  data.table(var="code", name="OPCS-3 code as it appears in the 'opcs_code' column of operations.csv"),
  data.table("meaning", "Description of OPCS-3 code"),
  data.table("display_order", "Numeric ordering of OPCS-3 codes for sorting the table"),
  data.table("parent_code", "Parent code of the OPCS-3  code"),
  data.table("appears_in_records", "TRUE where this OPCS-3 code is found in any hospital records, FALSE otherwise")
)
fwrite(info, "hospital_records/opcs3_codes_column_information.csv")

# Upload to persistent storage
system("dx upload hospital_records/* --destination 'common/Hospital Records/'")

# Send raw data to deletion folder to reduce storage costs - this needs to be 
# done in two steps, since we can't move two folders with the same name to the
# same location, so here we rename with a random number and then move to trash
rn <- as.integer(Sys.time())
system(sprintf("dx mv 'common/Hospital Records/raw_data/' 'common/Hospital Records/raw_data_%s'", rn), wait=TRUE)
system(sprintf("dx mv 'common/Hospital Records/raw_data_%s/' trash/", rn), wait=TRUE) 
