library(data.table)
library(lubridate)

# Pull down raw data that has been extracted with Table Exporter
system("mkdir -p raw_data", wait=TRUE)
system("dx download 'common/Blood Pressure/raw_data/data.csv' -o raw_data/blood_pressure.csv", wait=TRUE)

# Load in raw data and curated field information
raw <- fread("raw_data/blood_pressure.csv")
info <- fread("blood_pressure/field_information.csv")

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

# Set names
setnames(raw, c("p93", "p94", "p95", "p102", "p4079", "p4080"),
         c("sbp_manual", "dbp_manual", "pulse_rate_manual",
           "pulse_rate_auto", "dbp_automatic", "sbp_automatic"))

# Average repeated measures:
avg <- function(xx) {
  if(all(is.na(xx))) {
    return(NA_real_)
  } else {
    return(mean(na.omit(xx)))
  }
}

# SD of two repeated measures (used by QRISK3)
bpsd <- function(xx) {
  if (all(is.na(xx))) {
    return(NA_real_)
  } else {
    return(sd(na.omit(xx)))
  }
}

raw <- raw[, .(
  sbp_manual = avg(sbp_manual),
  dbp_manual = avg(dbp_manual),
  pulse_rate_manual = avg(pulse_rate_manual),
  sbp_automatic = avg(sbp_automatic),
  dbp_automatic = avg(dbp_automatic),
  pulse_rate_auto = avg(pulse_rate_auto),
  sd_sbp_manual = bpsd(sbp_manual),
  sd_sbp_automatic = bpsd(sbp_automatic)
), by=.(eid, visit_index)]

raw[, row := .I]
raw[is.na(sd_sbp_manual) & is.na(sd_sbp_automatic) & !is.na(sbp_manual) & !is.na(sbp_automatic), sd_sbp_manual := sd(c(sbp_manual, sbp_automatic)), by=row]
raw[is.na(sd_sbp_manual) & is.na(sd_sbp_automatic) & !is.na(sbp_manual) & !is.na(sbp_automatic), sd_sbp_automatic := sd_sbp_manual]

# Average automatic and manual measurements where both are present
raw <- raw[, .(
  sbp = avg(c(sbp_manual, sbp_automatic)),
  dbp = avg(c(dbp_manual, dbp_automatic)),
  pulse_rate = avg(c(pulse_rate_manual, pulse_rate_auto)),
  sd_sbp = avg(c(sd_sbp_manual, sd_sbp_automatic))
), by=.(eid, visit_index)]

# Curate information
info <- rbind(use.names=FALSE, fill=TRUE,
  data.table(var="eid", name="Application-specific Participant ID"),
  data.table(var="visit_index", name="UK Biobank assessment visit: 0 = Baseline Assessment, 1 = First repeat assessment, 2 = First imaging assessment, 3 = Second imaging assessment"),
  data.table(var="sbp", name="Systolic Blood Pressure (mmHg); average of manual (UKB Field ID: 93) and/or automatic measurements (UKB Field ID: 4080)"), 
  data.table(var="dbp", name="Diastolic Blood Pressure (mmHg); average of manual (UKB Field ID: 94) and/or automatic measurements (UKB Field ID: 4079)"),
  data.table(var="pulse_rate", name="Pulse Rate (beats per minute); average of manual (UKB Field ID: 95) and/or automatic measurements (UKB Field ID: 102)"),
  data.table(var="sd_sbp", name="Standard Deviation of SBP across the multiple manual (UKB Field ID: 93) and/or automatic measurements (UKB Field ID: 4080); used by QRISK3")
)

# Write out
fwrite(raw, quote=FALSE, file="blood_pressure/blood_pressure.csv")
fwrite(info, file="blood_pressure/column_information.csv")

# Upload to persistent storage
system("dx upload blood_pressure/blood_pressure.csv blood_pressure/column_information.csv --destination 'common/Blood Pressure/'", wait=TRUE)

# Send raw data to deletion folder to reduce storage costs - this needs to be 
# done in two steps, since we can't move two folders with the same name to the
# same location, so here we rename with a random number and then move to trash
rn <- as.integer(Sys.time())
system(sprintf("dx mv 'common/Blood Pressure/raw_data/' 'common/Blood Pressure/raw_data_%s'", rn), wait=TRUE)
system(sprintf("dx mv 'common/Blood Pressure/raw_data_%s/' trash/", rn), wait=TRUE) 
