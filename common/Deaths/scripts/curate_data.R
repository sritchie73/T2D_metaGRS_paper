library(data.table)

# create output directory on local cloud workstation
system("mkdir -p deaths")

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
system("dx download 'common/Deaths/raw_data/*' -o raw_data/")

# Curate death records
deaths <- fread("raw_data/death.csv")

setnames(deaths, "ins_index", "certificate_number")
deaths[, certificate_number := certificate_number + 1]

setnames(deaths, "dsource", "death_nation")
deaths[, death_nation := ifelse(death_nation == "E/W", "England/Wales", "Scotland")]

setnames(deaths, "source", "record_format")
deaths[, record_format := as.character(record_format)]
deaths[code_dict[coding_name == "data_coding_261"], on = .(record_format=code), record_format := meaning]

o <- unique(deaths[,.(eid)])
deaths <- deaths[order(-certificate_number)][o, on = .(eid)] # order so most recent certificate appears first

fwrite(deaths, "deaths/deaths.csv")

# Curate column information
info <- rbind(use.names=FALSE,
  data.table(var="dnx_death_id", name="Unique identifier for this death record"),
  data.table("eid", "Application-specific Participant ID"),
  data.table("certificate_number", "Death certificate number (sometimes more than one issued; e.g. due to postmortem)"),
  data.table("death_nation", "Deaths in England and Wales are reported by NHS England, while deaths in Scotland are reported by the NHS Central Register"),
  data.table("record_format", "Data format death record was recieved in"),
  data.table("date_of_death", "Date of death")
)

fwrite(info, "deaths/deaths_column_information.csv")

# Curate causes of death
causes <- fread("raw_data/death_cause.csv")

setnames(causes, "ins_index", "certificate_number")
causes[, certificate_number := certificate_number + 1]

setnames(causes, "arr_index", "cause_number")
causes[, cause_number := cause_number + 1]

setnames(causes, "level", "cause_type")
causes[, cause_type := ifelse(cause_type == 1, "primary", "secondary")]

causes <- causes[,.(dnx_death_id, dnx_death_cause_id, eid, certificate_number,
  cause_icd10, cause_type, cause_number)]

o <- unique(causes[,.(eid)])
causes <- causes[order(cause_number)][order(cause_type)][order(-certificate_number)][o, on = .(eid)]

# add in date of death and nation to minimise computationally expensive joins in
# downstream analyses, e.g. in common/Endpoints/curate_endpoints.R
causes[deaths, on = .(dnx_death_id), c("date_of_death", "death_nation") := .(date_of_death, death_nation)]

fwrite(causes, "deaths/death_causes.csv")

# Curate column information
info <- rbind(use.names=FALSE,
  data.table(var="dnx_death_id", name="Unique identifier for this death record"),
  data.table("dnx_death_cause_id", "Unique identifier for this cause of death for this death record"),
  data.table("eid", "Application-specific Participant ID"),
  data.table("certificate_number", "Death certificate number (sometimes more than one issued; e.g. due to postmortem)"),
  data.table("cause_icd10", "ICD 10 code for the cause of death"),
  data.table("cause_type", "Whether this ICD-10 code is the primary or secondary cause of death"),
  data.table("cause_number", "Order of cause of death as it appeared on the death certificate"),
  data.table("date_of_death", "Date of death"),
  data.table("death_nation", "Deaths in England and Wales are reported by NHS England, while deaths in Scotland are reported by the NHS Central Register")
)

fwrite(info, "deaths/death_causes_column_information.csv")

# Curate information about ICD-10 codes
icd10 <- code_dict[coding_name == "data_coding_19"]
icd10 <- icd10[, .(code, meaning, display_order, parent_code)]
icd10[, appears_in_records := FALSE]
icd10[causes, on = .(code=cause_icd10), appears_in_records := TRUE]
fwrite(icd10, "deaths/icd10_codes.csv")

# Curate column information
info <- rbind(use.names=FALSE,
  data.table(var="code", name="ICD-10 code as it appears in the 'cause_icd10' column of death_causes.csv"),
  data.table("meaning", "Description of ICD-10 code"),
  data.table("display_order", "Numeric ordering of ICD-10 codes for sorting the table"),
  data.table("parent_code", "Parent code of the ICD-10 code"),
  data.table("appears_in_records", "TRUE where this ICD-10 code is found in any death records, FALSE otherwise")
)

fwrite(info, "deaths/icd10_codes_column_information.csv")

# Upload to persistent storage
system("dx upload deaths/* --destination 'common/Deaths/'", wait=TRUE)

# Send raw data to deletion folder to reduce storage costs - this needs to be 
# done in two steps, since we can't move two folders with the same name to the
# same location, so here we rename with a random number and then move to trash
rn <- as.integer(Sys.time())
system(sprintf("dx mv 'common/Deaths/raw_data/' 'common/Deaths/raw_data_%s'", rn), wait=TRUE)
system(sprintf("dx mv 'common/Deaths/raw_data_%s/' trash/", rn), wait=TRUE) 
