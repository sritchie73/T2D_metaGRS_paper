library(data.table)
library(lubridate)

# Pull down raw data that has been extracted with Table Exporter
system("mkdir -p raw_data", wait=TRUE)
system("dx download 'common/Smoking/raw_data/data.csv' -o raw_data/smoking.csv", wait=TRUE)

# Load in raw data and curated field information
extracted <- fread("raw_data/smoking.csv")
info <- fread("Smoking/field_information.csv")

# Split out instance (visit) and array index (repeat measure) fields so they
# are rows instead of columns
visit_repeats <- setdiff(unique(gsub("^p[0-9]+_", "", names(extracted))), "eid")
extracted <- rbindlist(fill=TRUE, use.names=TRUE, lapply(visit_repeats, function(vr) {
  # Find columns matching this visit repeat pair (e.g. ending in _i0)
  this_cols <- names(extracted)[grepl(pattern=paste0(vr, "$"), names(extracted))]
  
  # Filter to these columns
  this_extracted <- extracted[, .SD, .SDcols=c("eid", this_cols)]
  
  # Drop repeat visit pair label from column name
  setnames(this_extracted, this_cols, gsub(paste0("_", vr, "$"), "", this_cols))
  
  # Add column for visit  index
  this_extracted[, visit_index := as.integer(gsub("^i", "", gsub("_.*", "", vr)))]

  # Move to start of data table
  this_extracted <- this_extracted[,.SD,.SDcols=c("eid", "visit_index", gsub("_.*", "", this_cols))]
  
  # Drop instance and array index combinations with all missing data
  # eid, visit_index, and array_index always non-missing
  this_extracted <- this_extracted[apply(this_extracted, 1, function(row) { sum(!is.na(row)) > 2L })]
  
  # Return
  this_extracted
}))

# Convert field ids to variable names
setnames(extracted, paste0("p", as.character(info$field.id)), info$var)

# Convert smoking status field codes to relevant labels
extracted[, smoking_status := fcase(
  smoking_status == 0, "Never",
  smoking_status == 1, "Previous",
  smoking_status == 2, "Current",
  smoking_status == -3, "Prefer not to answer",
  default = NA_character_)]

# Define current smoking
extracted[, current_smoker := fcase(
  smoking_status == "Current", TRUE,
  smoking_status == "Prefer not to answer", NA,
  default = FALSE)]

# Add number of cigarettes smoked per day
extracted[, daily_cigarettes := as.numeric(daily_cigarettes)]
extracted[!(current_smoker), daily_cigarettes := 0]
extracted[daily_cigarettes == -1, daily_cigarettes := NA] # Do not know
extracted[daily_cigarettes == -3, daily_cigarettes := NA] # Prefer not to answer
extracted[daily_cigarettes == -10, daily_cigarettes := 0.5] # Less than one a day

# Provide some sort of sensible ordering to columns
extracted <- extracted[, .(eid, visit_index, smoking_status, current_smoker, daily_cigarettes)]

# Add information on additional fields
info <- rbind(use.names=TRUE, fill=TRUE, info,
  data.table(var="eid", name="Application-specific Participant ID"),
  data.table(var="visit_index", name="UK Biobank assessment visit: 0 = Baseline Assessment, 1 = First repeat assessment, 2 = First imaging assessment, 3 = Second imaging assessment"),
  data.table(var="current_smoker", name="TRUE where 'smoking_status' == \"Current\", NA where participant answered \"Prefer not to answer\", FALSE otherwise.")
)

# Reorder rows and columns
info <- info[names(extracted), on = .(var), nomatch=0]
info <- info[, .(var, name, field.id)]

# Write out
fwrite(extracted, quote=FALSE, file="Smoking/smoking.csv")
fwrite(info, file="Smoking/column_information.csv")

# Upload to persistent storage
system("dx upload Smoking/smoking.csv Smoking/column_information.csv --destination 'common/Smoking/'", wait=TRUE)

# Send extracted data to deletion folder to reduce storage costs - this needs to be 
# done in two steps, since we can't move two folders with the same name to the
# same location, so here we rename with a random number and then move to trash
rn <- as.integer(Sys.time())
system(sprintf("dx mv 'common/Smoking/raw_data/' 'common/Smoking/raw_data_%s'", rn), wait=TRUE)
system(sprintf("dx mv 'common/Smoking/raw_data_%s/' trash/", rn), wait=TRUE) 
