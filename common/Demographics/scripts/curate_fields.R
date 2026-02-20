library(data.table)
library(foreach)

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
data_dict_file <- list.files(pattern="data_dictionary")
data_dict <- fread(data_dict_file)

# Manually curate field information we want to extract
info <- rbind(use.names=FALSE,
  data.table(field.id=54, name="Assessment centre", var="assessment_centre"),
  data.table(field.id=53, name="Date of attending assessment centre", var="assessment_date"),
  data.table(field.id=22189, name="Townsend deprivation index at recruitment", var="townsend"),
  data.table(field.id=31, name="Sex", var="sex"),
  data.table(field.id=21003, name="Age at assessment", var="age"),
  data.table(field.id=34, name="Year of birth", var="birth_year"),
  data.table(field.id=52, name="Month of birth", var="birth_month"),
  data.table(field.id=1647, name="Country of birth (UK/elsewhere)", var="birth_country"),
  data.table(field.id=20115, name="Country of birth (outside the UK)", var="birth_country_nonUK"),
  data.table(field.id=21000, name="Ethnicity", var="ethnicity"),
  data.table(field.id=3140, name="Pregnant", var="pregnant")
)

# Save table of information
system("mkdir -p Demographics")
fwrite(info, file="Demographics/field_information.csv")

# Extract list of fields for Table Exporter
fields <- foreach(this_field_id = info$field.id, .combine=rbind) %do% {
  data_dict[name %like% sprintf("p%s_", this_field_id) | name %like% sprintf("^p%s$", this_field_id)]
}

# Add in participant ID
fields <- rbind(data_dict[entity == "participant" & name == "eid"], fields)

# Write out field IDs for each distinct entity (only participant in this case,
# but generalizes later if we want to extract different entity types)
system("mkdir -p field_lists/")
for (entity_type in unique(fields$entity)) {
  fwrite(fields[entity == entity_type, .(name)], quote=FALSE, col.names=FALSE,
    file=sprintf("field_lists/demographics_fields_%s_entity.txt", entity_type))
}

# Upload field list to persistent storage
system("dx upload field_lists/demographics_fields_* --destination 'common/Demographics/scripts/'", wait=TRUE)