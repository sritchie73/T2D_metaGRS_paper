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
  data.table(field.id=31, name="Sex"),
  data.table(field.id=6177, name="Touchscreen medication question for males"),
  data.table(field.id=6153, name="Touchscreen medication question for females"),
  data.table(field.id=2492, name="Touchscreen question: are you taking other prescription medications?"),
  data.table(field.id=137, name="Verbal interview with nurse: how many treatments/medications are you taking?"),
  data.table(field.id=20003, name="Verbal interview with nurse: treatment/medication code"),
  data.table(field.id=6671, name="Verbal interview with nurse: how many antibiotics have you taken in the last 3 months?"),
  data.table(field.id=20199, name="Verbal interview with nurse: antibiotic codes for last 3 months")
)

# Save table of information
system("mkdir -p Medication")
fwrite(info, file="Medication/field_information.csv")

# Extract list of fields for Table Exporter
fields <- foreach(this_field_id = info$field.id, .combine=rbind) %do% {
  data_dict[name %like% sprintf("^p%s(_|$)", this_field_id)]
}

# Add in participant ID
fields <- rbind(data_dict[entity == "participant" & name == "eid"], fields)

# Write out field IDs for each distinct entity (only participant in this case,
# but generalizes later if we want to extract different entity types)
system("mkdir -p field_lists/")
for (entity_type in unique(fields$entity)) {
  fwrite(fields[entity == entity_type, .(name)], quote=FALSE, col.names=FALSE,
    file=sprintf("field_lists/medication_fields_%s_entity.txt", entity_type))
}

# Upload field list to persistent storage
system("dx upload field_lists/medication_fields_* --destination 'common/Medications/scripts/'", wait=TRUE)