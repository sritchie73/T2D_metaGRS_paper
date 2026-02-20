library(data.table)
library(foreach)

# Pull down raw data that has been raw with Table Exporter
system("mkdir -p raw_data", wait=TRUE)
system("dx download 'common/Family History/raw_data/data.csv' -o raw_data/family_history.csv", wait=TRUE)

# Load in raw data and curated field information
raw <- fread("raw_data/family_history.csv", na.strings=c("", "NA"))
info <- fread("family_history/field_information.csv")

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

# Melt to long and split out multiple selected answers
raw <- melt(raw, id.vars=c("eid", "visit_index"))
raw <- raw[, .(value=strsplit(value, "\\|")[[1]]), by=.(eid, visit_index, variable)]

# Assign labels to codes
raw[, value := fcase(
  value == 13, "Prostate cancer",
  value == 12, "Severe depression",
  value == 11, "Parkinson's disease",
  value == 10, "Alzheimer's disease/dementia",
  value == 9, "Diabetes",
  value == 8, "High blood pressure",
  value == 6, "Chronic bronchitis/emphysema",
  value == 5, "Breast cancer",
  value == 4, "Bowel cancer",
  value == 3, "Lung cancer",
  value == 2, "Stroke",
  value == 1, "Heart disease",
  value == -11, "Do not know (group 1)",
  value == -13, "Prefer not to answer (group 1)",
  value == -17, "None of the above (group 1)",
  value == -21, "Do not know (group 2)",
  value == -23, "Prefer not to answer (group 2)",
  value == -27, "None of the above (group 2)"
)]

# Build wide datasets for each disease
datasets <- foreach(this_group = c("father", "mother", "siblings")) %do% {
  this_dat <- unique(raw[,.(eid, visit_index)])
  this_raw <- raw[variable == this_group]
  
  this_dat[, no_disease := FALSE]
  none <- this_raw[value == "None of the above (group 1)", .(eid, visit_index)]
  none <- none[this_raw[value == "None of the above (group 2)", .(eid, visit_index)], on = .(eid, visit_index), nomatch=0]
  this_dat[none, on = .(eid, visit_index), no_disease := TRUE]
  this_dat[this_raw[value %like% "Do not know" | value %like% "Prefer not to answer"], on = .(eid, visit_index), no_disease := NA]
  
  this_dat[, heart_disease := FALSE]
  this_dat[this_raw[value == "Heart disease"], on = .(eid, visit_index), heart_disease := TRUE]
  this_dat[this_raw[value %in% c("Do not know (group 1)", "Prefer not to answer (group 1)")], on = .(eid, visit_index), heart_disease := NA]
  
  this_dat[, stroke := FALSE]
  this_dat[this_raw[value == "Stroke"], on = .(eid, visit_index), stroke := TRUE]
  this_dat[this_raw[value %in% c("Do not know (group 1)", "Prefer not to answer (group 1)")], on = .(eid, visit_index), stroke := NA]
  
  this_dat[, cardiovascular_disease := heart_disease | stroke]
  
  this_dat[, copd := FALSE]
  this_dat[this_raw[value == "Chronic bronchitis/emphysema"], on = .(eid, visit_index), copd := TRUE]
  this_dat[this_raw[value %in% c("Do not know (group 1)", "Prefer not to answer (group 1)")], on = .(eid, visit_index), copd := NA]
  
  this_dat[, hypertension := FALSE]
  this_dat[this_raw[value == "High blood pressure"], on = .(eid, visit_index), hypertension := TRUE]
  this_dat[this_raw[value %in% c("Do not know (group 1)", "Prefer not to answer (group 1)")], on = .(eid, visit_index), hypertension := NA]
  
  this_dat[, diabetes := FALSE]
  this_dat[this_raw[value == "Diabetes"], on = .(eid, visit_index), diabetes := TRUE]
  this_dat[this_raw[value %in% c("Do not know (group 1)", "Prefer not to answer (group 1)")], on = .(eid, visit_index), diabetes := NA]
  
  this_dat[, alzheimers := FALSE]
  this_dat[this_raw[value == "Alzheimer's disease/dementia"], on = .(eid, visit_index), alzheimers := TRUE]
  this_dat[this_raw[value %in% c("Do not know (group 1)", "Prefer not to answer (group 1)")], on = .(eid, visit_index), alzheimers := NA]
  
  this_dat[, parkinsons := FALSE]
  this_dat[this_raw[value == "Parkinson's disease"], on = .(eid, visit_index), parkinsons := TRUE]
  this_dat[this_raw[value %in% c("Do not know (group 2)", "Prefer not to answer (group 2)")], on = .(eid, visit_index), parkinsons := NA]
  
  this_dat[, depression := FALSE]
  this_dat[this_raw[value == "Severe depression"], on = .(eid, visit_index), depression := TRUE]
  this_dat[this_raw[value %in% c("Do not know (group 2)", "Prefer not to answer (group 2)")], on = .(eid, visit_index), depression := NA]
  
  this_dat[, lung_cancer := FALSE]
  this_dat[this_raw[value == "Lung cancer"], on = .(eid, visit_index), lung_cancer := TRUE]
  this_dat[this_raw[value %in% c("Do not know (group 2)", "Prefer not to answer (group 2)")], on = .(eid, visit_index), lung_cancer := NA]
  
  this_dat[, bowel_cancer := FALSE]
  this_dat[this_raw[value == "Bowel cancer"], on = .(eid, visit_index), bowel_cancer := TRUE]
  this_dat[this_raw[value %in% c("Do not know (group 2)", "Prefer not to answer (group 2)")], on = .(eid, visit_index), bowel_cancer := NA]
  
  this_dat[, breast_cancer := FALSE]
  if (this_group != "father") {
    this_dat[this_raw[value == "Breast cancer"], on = .(eid, visit_index), breast_cancer := TRUE]
    this_dat[this_raw[value %in% c("Do not know (group 2)", "Prefer not to answer (group 2)")], on = .(eid, visit_index), breast_cancer := NA]
  }

  this_dat[, prostate_cancer := FALSE]
  if (this_group != "mother") {
    this_dat[this_raw[value == "Prostate cancer"], on = .(eid, visit_index), prostate_cancer := TRUE]
    this_dat[this_raw[value %in% c("Do not know (group 2)", "Prefer not to answer (group 2)")], on = .(eid, visit_index), prostate_cancer := NA]
  }
  
  this_dat[, any_cancer := lung_cancer | bowel_cancer | breast_cancer | prostate_cancer]
  
  return(this_dat)
}
names(datasets) <- c("father", "mother", "siblings")

# Write out
fwrite(datasets$father, quote=FALSE, file="family_history/illness_of_father.csv")
fwrite(datasets$mother, quote=FALSE, file="family_history/illness_of_mother.csv")
fwrite(datasets$siblings, quote=FALSE, file="family_history/illness_of_siblings.csv")

# Combine father and mother tables to create illness of parents table
p1 <- melt(datasets$father, id.vars=c("eid", "visit_index"), value.name="father")
p2 <- melt(datasets$mother, id.vars=c("eid", "visit_index"), value.name="mother")
parents <- p1[p2, on = .(eid, visit_index, variable)]
parents[variable != "no_disease", parents := father | mother]
parents[variable == "no_disease", parents := father & mother]
parents <- dcast(parents, eid + visit_index ~ variable, value.var="parents")
parents <- parents[, .SD, .SDcols=names(datasets$father)]
parents <- parents[unique(raw[,.(eid, visit_index)]), on = .(eid, visit_index)]
fwrite(parents, quote=FALSE, file="family_history/illness_of_parents.csv")

# Combine with siblings table to create illness of first degree relatives table
p1 <- melt(datasets$father, id.vars=c("eid", "visit_index"), value.name="father")
p2 <- melt(datasets$mother, id.vars=c("eid", "visit_index"), value.name="mother")
p3 <- melt(datasets$siblings, id.vars=c("eid", "visit_index"), value.name="siblings")
relatives <- p1[p2, on = .(eid, visit_index, variable)][p3, on = .(eid, visit_index, variable)]
relatives[variable != "no_disease", relatives := father | mother | siblings]
relatives[variable == "no_disease", relatives := father & mother & siblings]
relatives <- dcast(relatives, eid + visit_index ~ variable, value.var="relatives")
relatives <- relatives[, .SD, .SDcols=names(datasets$father)]
relatives <- relatives[unique(raw[,.(eid, visit_index)]), on = .(eid, visit_index)]
fwrite(parents, quote=FALSE, file="family_history/illness_of_first_degree_relatives.csv")

# Upload to persistent storage
system("dx upload family_history/illness_of_* --destination 'common/Family History/'", wait=TRUE)

# Send raw data to deletion folder to reduce storage costs - this needs to be 
# done in two steps, since we can't move two folders with the same name to the
# same location, so here we rename with a random number and then move to trash
rn <- as.integer(Sys.time())
system(sprintf("dx mv 'common/Family History/raw_data/' 'common/Family History/raw_data_%s'", rn), wait=TRUE)
system(sprintf("dx mv 'common/Family History/raw_data_%s/' trash/", rn), wait=TRUE) 
