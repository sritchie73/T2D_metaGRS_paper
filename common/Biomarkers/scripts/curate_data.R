library(data.table)

# Pull down raw data that has been extracted with Table Exporter
system("mkdir -p raw_data", wait=TRUE)
system("dx download 'common/Biomarkers/raw_data/data.csv' -o raw_data/biomarkers.csv", wait=TRUE)

# Load pre-curated information sheets
biomarker_info <- fread("biomarkers/biomarker_information.csv")
sample_info <- fread("biomarkers/sample_information.csv")

# Load in raw data
raw <- fread("raw_data/biomarkers.csv")

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

# Extract blood biomarkers
blood <- raw[, .SD, .SDcols=c("eid", "visit_index",
  biomarker_info[sample_type != "Urine" & !is.na(UKB.Field.ID), paste0("p", UKB.Field.ID)])]

setnames(blood,
   biomarker_info[!is.na(UKB.Field.ID) & sample_type != "Urine", paste0("p", UKB.Field.ID)],
   biomarker_info[!is.na(UKB.Field.ID) & sample_type != "Urine", var])

# Extract blood biomarker missingness reason where missing due to being above or below
# detection limits (Reportability field, coding: https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=4917)
miss <- raw[, .SD, .SDcols=c("eid", "visit_index",
  biomarker_info[sample_type != "Urine" & !is.na(Reportability.Field.ID), paste0("p", Reportability.Field.ID)])]

setnames(miss,
  biomarker_info[!is.na(Reportability.Field.ID) & sample_type != "Urine", paste0("p", Reportability.Field.ID)],
  biomarker_info[!is.na(Reportability.Field.ID) & sample_type != "Urine", var])

# Melt to long to determine lower/upper detection limits
blood <- melt(blood, id.vars=c("eid", "visit_index"))
miss <- melt(miss, id.vars=c("eid", "visit_index"), na.rm=TRUE)

blood[, below_detection := FALSE]
blood[miss[value %in% c(2, 4)], on = .(eid, visit_index, variable), below_detection := TRUE]

blood[, above_detection := FALSE]
blood[miss[value %in% c(3, 5)], on = .(eid, visit_index, variable), above_detection := TRUE]

# Derive non-hdl
nonhdl <- merge(blood[variable == "tchol"], blood[variable == "hdl"],
  by=c("eid", "visit_index"), suffixes=c(".tchol", ".hdl"))
nonhdl[, value := value.tchol - value.hdl]
nonhdl[, below_detection := FALSE]
nonhdl[, above_detection := FALSE]

# Upper and lower limits interact in different ways:
#
# Case numbers:
#
#    below_detection.tchol above_detection.tchol below_detection.hdl above_detection.hdl      N
# 1:                 FALSE                 FALSE               FALSE               FALSE 445445
# 2:                  TRUE                 FALSE                TRUE               FALSE      4
# 3:                 FALSE                 FALSE               FALSE                TRUE      1
# 4:                 FALSE                 FALSE                TRUE               FALSE      1
#
# Types:
#  
# 1. Both tchol and hdl are within detection limits, nonhdl is a precise estimate
# 2. Both tchol and hdl are below detection limits, the true value is in range [nonhdl, tchol].
# 3. hdl is above detection, then the true value is in range [0, nonhdl].
# 4. hdl is below detection, then the true value is in range [nonhdl, tchol].
#
nonhdl[below_detection.tchol & below_detection.hdl, above_detection := TRUE]
nonhdl[!(below_detection.tchol) & above_detection.hdl, below_detection := TRUE]
nonhdl[!(below_detection.tchol) & below_detection.hdl, above_detection := TRUE]

# add to main table
nonhdl <- nonhdl[, .(eid, visit_index, variable="nonhdl", value, below_detection, above_detection)]
blood <- rbind(blood, nonhdl)

# Derive ratio of ApoB to ApoA1
apobapoa1 <- merge(blood[variable == "apob"], blood[variable == "apoa1"],
                   by=c("eid", "visit_index"), suffixes=c(".apob", ".apoa1"))
apobapoa1[, value := value.apob / value.apoa1]
apobapoa1[, below_detection := FALSE]
apobapoa1[, above_detection := FALSE]

# Upper and lower limits interact in different ways:
# 
# Case numbers:
#
#    below_detection.apob above_detection.apob below_detection.apoa1 above_detection.apoa1      N
# 1:                FALSE                FALSE                 FALSE                 FALSE 440741
# 2:                FALSE                FALSE                 FALSE                  TRUE   1792
# 3:                FALSE                 TRUE                 FALSE                 FALSE    306
# 4:                 TRUE                FALSE                 FALSE                 FALSE    730
# 5:                FALSE                FALSE                  TRUE                 FALSE     10
# 6:                 TRUE                FALSE                 FALSE                  TRUE     11
# 7:                 TRUE                FALSE                  TRUE                 FALSE     10
# 8:                FALSE                 TRUE                 FALSE                  TRUE      2
# 
# Types:
#
# 1. Both apob and apoa1 are within detection limits, apobapoa1 is a precise ratio
# 2. Only apoa1 is above detection, then the true value is in range [0, apobapoa1]
# 3. Only apob is above detection, then the true value is in range [apobapoa1, inf)
# 4. Only apob is below detection, then the true value is in range [0, apobapoa1]
# 5. Only apoa1 is below detection, then the true value is in range [apobapoa1, inf)
# 6. apob is below detection and apoa1 is above detection, then the true value is in range [0, apobapoa1]
# 7. apob and apoa1 are both below detection, then the ratio is an imprecise estimate, but could truly be anywhere between [0, Inf).
# 8. apob and apoa1 are both above detection, then the ratio is an imprecise estimate, but could truly be anywhere between [0, Inf).
apobapoa1[, below_detection := fcase(
  !is.na(value) & above_detection.apoa1, TRUE,
  !is.na(value) & below_detection.apob, TRUE,
  !is.na(value) & below_detection.apob & above_detection.apoa1, TRUE,
  below_detection.apoa1 & below_detection.apob, FALSE,
  above_detection.apoa1 & above_detection.apob, FALSE,
  default = FALSE)]
apobapoa1[, above_detection := fcase(
  !is.na(value) & below_detection.apoa1, TRUE,
  !is.na(value) & above_detection.apob, TRUE,
  below_detection.apoa1 & below_detection.apob, FALSE,
  above_detection.apoa1 & above_detection.apob, FALSE,
  default = FALSE)]

# add to main table
apobapoa1 <- apobapoa1[, .(eid, visit_index, variable="apobapoa1", value, below_detection, above_detection)]
blood <- rbind(blood, apobapoa1)

# Set values above/below detection to their respective limits (with small offset of 0.0001)
lims <- blood[!is.na(value), .(min = min(value), max = max(value)), by=variable]
blood <- rbind(
  blood[!(below_detection) & !(above_detection)],
  blood[(below_detection)][lims, on = .(variable), nomatch=0, .(
    eid, visit_index, variable, value = i.min - 0.0001,
    below_detection, above_detection
  )],
  blood[(above_detection)][lims, on = .(variable), nomatch=0, .(
    eid, visit_index, variable, value = i.max + 0.0001,
    below_detection, above_detection
  )]
)

# Define and add hba1c_pct
hba1c_pct <- blood[variable == "hba1c"]
hba1c_pct[, value := value/10.929+2.15]
hba1c_pct[, variable := "hba1c_pct"]
blood <- rbind(blood, hba1c_pct)

# Drop missing values
blood <- blood[!is.na(value)]

# Extract urine biomarkers
urine <- raw[, .SD, .SDcols=c("eid", "visit_index",
  biomarker_info[sample_type == "Urine" & !is.na(UKB.Field.ID), paste0("p", UKB.Field.ID)])]

setnames(urine,
  biomarker_info[!is.na(UKB.Field.ID) & sample_type == "Urine", paste0("p", UKB.Field.ID)],
  biomarker_info[!is.na(UKB.Field.ID) & sample_type == "Urine", var])

# Extract urine biomarker missingness reason where missing due to being above or below
# detection limits (Reportability field, coding: https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=4917)
miss <- raw[, .SD, .SDcols=c("eid", "visit_index",
  biomarker_info[sample_type == "Urine" & !is.na(Reportability.Field.ID), paste0("p", Reportability.Field.ID)])]

setnames(miss,
  biomarker_info[!is.na(Reportability.Field.ID) & sample_type == "Urine", paste0("p", Reportability.Field.ID)],
  biomarker_info[!is.na(Reportability.Field.ID) & sample_type == "Urine", var])

# Melt to long to determine lower/upper detection limits
urine[, uriacc := as.numeric(uriacc)] # otherwise integer column throws warning in melt
urine <- melt(urine, id.vars=c("eid", "visit_index"))
miss <- melt(miss, id.vars=c("eid", "visit_index"), na.rm=TRUE)
miss <- miss[value != ""]

urine[, below_detection := FALSE]
urine[miss[value %like% "<"], on = .(eid, visit_index, variable), below_detection := TRUE]

urine[, above_detection := FALSE]
urine[miss[value %like% ">"], on = .(eid, visit_index, variable), above_detection := TRUE]

# Set values above/below detection to their respective limits (with small offset of 0.0001)
# Unlike blood biomarkers, limits are hard coded into the reportability fields
lims <- unique(miss[,.(variable, value)])
lower <- lims[value %like% "<", .(variable, min=as.numeric(gsub("<", "", value)))]
upper <- lims[value %like% ">", .(variable, max=as.numeric(gsub(">", "", value)))]
lims <- merge(lower, upper, by = "variable", all=TRUE)

urine <- rbind(
  urine[!(below_detection) & !(above_detection)],
  urine[(below_detection)][lims, on = .(variable), nomatch=0, .(
    eid, visit_index, variable, value = i.min - 0.0001,
    below_detection, above_detection
  )],
  urine[(above_detection)][lims, on = .(variable), nomatch=0, .(
    eid, visit_index, variable, value = i.max + 0.0001,
    below_detection, above_detection
  )]
)

# Drop missing
urine <- urine[!is.na(value)]

# Build combined biomarker table
biomarkers <- rbind(blood, urine)

# Build and save table of information about measurements that were either above 
# or below detection limits
limits <- biomarkers[(below_detection) | (above_detection),
   .(eid, visit_index, variable, type = ifelse(below_detection, "below detection limit", "above detection limit"))]
fwrite(limits, file="biomarkers/measurements_outside_detection.csv")

# Cast biomarker table back to wide format
biomarkers <- dcast(biomarkers, eid + visit_index ~ variable, value.var="value")
biomarkers <- biomarkers[raw[,.(eid, visit_index)], on = .(eid, visit_index), nomatch=0]

# Add in QC columns flagging whether missing values are due to samples that were 
# not measured (as opposed to missing data for some other reason)
not_measured <- raw[, .(eid, visit_index,
  no_blood_sample = ifelse(!is.na(p20050), TRUE, FALSE),
  no_urine_sample = ifelse(!is.na(p20072), TRUE, FALSE))]
not_measured <- not_measured[(no_blood_sample) | (no_urine_sample)] # 1 million rows without filtering - i.e. incorrectly counts everyone as having repeat assessment
biomarkers <- merge(not_measured, biomarkers, all=TRUE, by=c("eid", "visit_index"))
biomarkers[is.na(no_blood_sample), no_blood_sample := FALSE]
biomarkers[is.na(no_urine_sample), no_urine_sample := FALSE]

# Filter sample info columns
sample_info <- sample_info[!(var %in% c("no_blood_sample_reason", "no_urine_sample_reason"))]

# Write out sample information, filtering to fields kept
fwrite(sample_info, file="biomarkers/sample_information.csv")

# Add in fasting time
biomarkers[raw, on = .(eid, visit_index), fasting_time := p74]

# Reorganise columns
biomarkers <- biomarkers[, .SD, .SDcols=c(sample_info$var, biomarker_info$var)]

# Write out
fwrite(biomarkers, quote=FALSE, file="biomarkers/biomarkers.csv")

# Upload to persistent storage
system("dx upload biomarkers/biomarkers.csv biomarkers/measurements_outside_detection.csv --destination 'common/Biomarkers/'", wait=TRUE)
system("dx upload biomarkers/sample_information.csv biomarkers/biomarker_information.csv --destination 'common/Biomarkers/'", wait=TRUE)

# Send raw data to deletion folder to reduce storage costs - this needs to be 
# done in two steps, since we can't move two folders with the same name to the
# same location, so here we rename with a random number and then move to trash
rn <- as.integer(Sys.time())
system(sprintf("dx mv 'common/Biomarkers/raw_data/' 'common/Biomarkers/raw_data_%s'", rn), wait=TRUE)
system(sprintf("dx mv 'common/Biomarkers/raw_data_%s/' trash/", rn), wait=TRUE) 
