library(data.table)
library(lubridate)

# Pull down raw data that has been extracted with Table Exporter
system("mkdir -p raw_data", wait=TRUE)
system("dx download 'common/Anthropometrics/raw_data/data.csv' -o raw_data/anthropometrics.csv", wait=TRUE)

# Load in raw data and curated field information
raw <- fread("raw_data/anthropometrics.csv")
info <- fread("Anthropometrics/field_information.csv")

# Split out instance (visit) fields so they are rows instead of columns
visit_repeats <- setdiff(unique(gsub("^p[0-9]+_", "", names(raw))), "eid")
raw <- rbindlist(fill=TRUE, use.names=TRUE, lapply(visit_repeats, function(vr) {
  # Find columns matching this visit repeat pair (e.g. ending in _i0)
  this_cols <- names(raw)[grepl(pattern=paste0(vr, "$"), names(raw))]
  
  # Filter to these columns
  this_raw <- raw[, .SD, .SDcols=c("eid", this_cols)]
  
  # Drop repeat visit pair label from column name
  setnames(this_raw, this_cols, gsub(paste0("_", vr, "$"), "", this_cols))
  
  # Add columns for visit index
  this_raw[, visit_index := as.integer(gsub("^i", "", gsub("_.*", "", vr)))]
  
  # Move to start of data table
  this_raw <- this_raw[,.SD,.SDcols=c("eid", "visit_index", gsub("_.*", "", this_cols))]
  
  # Drop instance and array index combinations with all missing data
  # eid, visit_index, and array_index always non-missing
  this_raw <- this_raw[apply(this_raw, 1, function(row) { sum(!is.na(row)) > 2L })]
  
  # Return
  this_raw
}))

# Convert field ids to variable names
setnames(raw, paste0("p", as.character(info$field.id)), info$var)

# Compute waist to hip ratio
raw[, waist_hip_ratio := waist / hip]

# Provide some sort of sensible ordering to columns
raw <- raw[, .(eid, visit_index, height, weight, bmi, waist, hip, waist_hip_ratio)]

# Add information on additional fields
info <- rbind(use.names=TRUE, fill=TRUE, info,
  data.table(var="eid", name="Application-specific Participant ID"),
  data.table(var="visit_index", name="UK Biobank assessment visit: 0 = Baseline Assessment, 1 = First repeat assessment, 2 = First imaging assessment, 3 = Second imaging assessment"),
  data.table(var="waist_hip_ratio", name="Waist to hip ratio, computed as waist / hip")
)

# Reorder rows and columns
info <- info[names(raw), on = .(var), nomatch=0]
info <- info[, .(var, name, field.id)]

# Write out
fwrite(raw, quote=FALSE, file="Anthropometrics/anthropometrics.csv")
fwrite(info, file="Anthropometrics/column_information.csv")

# Upload to persistent storage
system("dx upload Anthropometrics/anthropometrics.csv Anthropometrics/column_information.csv --destination 'common/Anthropometrics/'", wait=TRUE)

# Send raw data to deletion folder to reduce storage costs - this needs to be 
# done in two steps, since we can't move two folders with the same name to the
# same location, so here we rename with a random number and then move to trash
rn <- as.integer(Sys.time())
system(sprintf("dx mv 'common/Anthropometrics/raw_data/' 'common/Anthropometrics/raw_data_%s'", rn), wait=TRUE)
system(sprintf("dx mv 'common/Anthropometrics/raw_data_%s/' trash/", rn), wait=TRUE) 
