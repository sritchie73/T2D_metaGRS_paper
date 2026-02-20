library(data.table)
library(lubridate)

# Pull down raw data that has been extracted with Table Exporter
system("mkdir -p raw_data", wait=TRUE)
system("dx download 'common/Demographics/raw_data/data.csv' -o raw_data/demographics.csv", wait=TRUE)

# Load in raw data and curated field information
raw <- fread("raw_data/demographics.csv")
info <- fread("Demographics/field_information.csv")

# For fields only available at baseline, add _i0 to the field name for processing
only_baseline <- names(raw)[!(names(raw) %like% "_") & names(raw) %like% "^p"]
setnames(raw, only_baseline, paste0(only_baseline, "_i0"))

# Split out instances from field names so they are rows instead of columns
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

# Some fields (e.g. sex) only available at baseline, copy to repeat visits
raw[raw[visit_index == 0], on = .(eid), sex := i.sex]
raw[raw[visit_index == 0], on = .(eid), birth_year := i.birth_year]
raw[raw[visit_index == 0], on = .(eid), birth_month := i.birth_month]
raw[raw[visit_index == 0], on = .(eid), townsend := i.townsend]
raw[raw[visit_index == 0], on = .(eid), birth_country_nonUK := i.birth_country_nonUK]

# Ethnicity is missing for visit 3, and for 72% of participants at visit 2.
# For these, propogate the most recent reported ethnicity. Note, in some 
# instances a persons self-reported ethnicity changes between visits.
eth_miss_v2 <- raw[visit_index == 2 & is.na(ethnicity)]
eth_miss_v2[raw[order(visit_index)][visit_index < 2][, .SD[.N], by=.(eid)], on = .(eid), ethnicity := i.ethnicity]
raw[eth_miss_v2, on = .(eid, visit_index), ethnicity := i.ethnicity]

eth_miss_v3 <- raw[visit_index == 3 & is.na(ethnicity)]
eth_miss_v3[raw[order(visit_index)][visit_index < 3][, .SD[.N], by=.(eid)], on = .(eid), ethnicity := i.ethnicity]
raw[eth_miss_v3, on = .(eid, visit_index), ethnicity := i.ethnicity]

# Similar story for country of birth
cb_miss_v2 <- raw[visit_index == 2 & is.na(birth_country)]
cb_miss_v2[raw[order(visit_index)][visit_index < 2][, .SD[.N], by=.(eid)], on = .(eid), birth_country := i.birth_country]
raw[cb_miss_v2, on = .(eid, visit_index), birth_country:= i.birth_country]

cb_miss_v3 <- raw[visit_index == 3 & is.na(birth_country)]
cb_miss_v3[raw[order(visit_index)][visit_index < 3][, .SD[.N], by=.(eid)], on = .(eid), birth_country := i.birth_country]
raw[cb_miss_v3, on = .(eid, visit_index), birth_country := i.birth_country]

# Convert codes for sex
raw[, sex := ifelse(sex == 1L, "Male", "Female")]

# Convert codes for pregnancy
raw[, pregnant := fcase(
  sex == "Male", "No",
  pregnant == 0L, "No",
  pregnant == 1L, "Yes",
  pregnant == 2L, "Unsure"
)]

# Convert assessment centre codes to locations
raw[, assessment_centre := fcase(
  assessment_centre == "11012", "Barts",
  assessment_centre == "11021", "Birmingham",
  assessment_centre == "11011", "Bristol",
  assessment_centre == "11008", "Bury",
  assessment_centre == "11003", "Cardiff",
  assessment_centre == "11024", "Cheadle (revisit)",
  assessment_centre == "11020", "Croydon",
  assessment_centre == "11005", "Edinburgh",
  assessment_centre == "11004", "Glasgow",
  assessment_centre == "11018", "Hounslow",
  assessment_centre == "11010", "Leeds",
  assessment_centre == "11016", "Liverpool",
  assessment_centre == "11001", "Manchester",
  assessment_centre == "11017", "Middlesborough",
  assessment_centre == "11009", "Newcastle",
  assessment_centre == "11013", "Nottingham",
  assessment_centre == "11002", "Oxford",
  assessment_centre == "11007", "Reading",
  assessment_centre == "11014", "Sheffield",
  assessment_centre == "10003", "Stockport (pilot)",
  assessment_centre == "11006", "Stoke",
  assessment_centre == "11022", "Swansea",
  assessment_centre == "11023", "Wrexham",
  assessment_centre == "11025", "Cheadle (imaging)",
  assessment_centre == "11026", "Reading (imaging)",
  assessment_centre == "11027", "Newcastle (imaging)",
  assessment_centre == "11028", "Bristol (imaging)"
)]

# Code Nation of assessment centre
raw[, assessment_nation := "England"]
raw[assessment_centre %in% c("Cardiff", "Swansea", "Wrexham"), assessment_nation := "Wales"]
raw[assessment_centre %in% c("Edinburgh", "Glasgow"), assessment_nation := "Scotland"]

# Convert ethnicity codes to labels. Codes and labels are organised
# into a tree, with codes in the 1000s belonging to top level codes
# starting with the same number. E.g. "British" (1001), "Irish" (1002),
# "Any other white background" (1003) are all sub-groups of "White" (1).
# https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=1001
raw[, ethnicity_subgroup := fcase(
  ethnicity == 1, "White",
  ethnicity == 1001, "British",
  ethnicity == 2001, "White and Black Caribbean",
  ethnicity == 3001, "Indian",
  ethnicity == 4001, "Caribbean",
  ethnicity == 2, "Mixed",
  ethnicity == 1002, "Irish",
  ethnicity == 2002, "White and Black African",
  ethnicity == 3002, "Pakistani",
  ethnicity == 4002, "African",
  ethnicity == 3, "Asian or Asian British",
  ethnicity == 1003, "Any other white background",
  ethnicity == 2003, "White and Asian",
  ethnicity == 3003, "Bangladeshi",
  ethnicity == 4003, "Any other Black background",
  ethnicity == 4, "Black or Black British",
  ethnicity == 2004, "Any other mixed background",
  ethnicity == 3004, "Any other Asian background",
  ethnicity == 5, "Chinese",
  ethnicity == 6, "Other ethnic group",
  ethnicity == -1, "Do not know",
  ethnicity == -3, "Prefer not to answer"
)]

raw[, ethnicity_group := fcase(
  ethnicity %in% c(1, 1001, 1002, 1003), "White",
  ethnicity %in% c(2, 2001, 2002, 2003, 2004), "Mixed",
  ethnicity %in% c(3, 3001, 3002, 3003, 3004), "Asian or Asian British",
  ethnicity %in% c(4, 4001, 4002, 4003), "Black or Black British",
  ethnicity == 5, "Chinese",
  ethnicity == 6, "Other ethnic group",
  ethnicity == -1, "Do not know",
  ethnicity == -3, "Prefer not to answer"
)]

# Compute decimal age which incorporates birth month into age
# Note day of birth not provided for privacy reasons, so we guess
# this as the midpoint of each month
raw[, approx_birth_date := as.IDate(sprintf("%s-%s-%s", birth_year, birth_month, 15))]
raw[, age_decimal := time_length(as.Date(assessment_date) - as.Date(approx_birth_date), unit="years")]
raw[, age_decimal := round(age_decimal, digits=2)] # any further resolution is meaningless due to lack of precision

# Load in codings map

# Code country of birth variable
raw[, birth_country := fcase(
  birth_country == 1, "England",
  birth_country == 2, "Wales",
  birth_country == 3, "Scotland",
  birth_country == 4, "Northern Ireland",
  birth_country == 5, "Republic of Ireland",
  birth_country == 6, "Elsewhere",
  birth_country == -1, "Do not know",
  birth_country == -3, "Prefer not to answer"
)]

code_dict_file <- list.files(pattern="codings")
code_dict <- fread(code_dict_file)
raw[, birth_country_nonUK := as.character(birth_country_nonUK)]
raw[code_dict[coding_name == "data_coding_89"], on = .(birth_country_nonUK=code), birth_country_nonUK := meaning]

raw[birth_country == "Elsewhere", birth_country := birth_country_nonUK] # Combine into single column

# Provide some sort of sensible ordering to columns
raw <- raw[, .(eid, visit_index, assessment_date, assessment_centre, 
  assessment_nation, townsend, sex, pregnant, ethnicity_group, ethnicity_subgroup, 
  birth_country, age, age_decimal, birth_year, birth_month, approx_birth_date)]

# Add information on additional fields
info <- rbind(use.names=TRUE, fill=TRUE, info,
  data.table(var="eid", name="Application-specific Participant ID"),
  data.table(var="visit_index", name="UK Biobank assessment visit: 0 = Baseline Assessment, 1 = First repeat assessment, 2 = First imaging assessment, 3 = Second imaging assessment"),
  data.table(var="assessment_nation", name="Nation assessment centre is located in"),
  data.table(var="age_decimal", name="Age compute using birth date and birth month"),
  data.table(var="approx_birth_date", name="Date of birth, inferred as the 15th of the 'birth_month' in the 'birth_year'"),
  data.table(var="ethnicity_group", name="Top-level ethnicity grouping in field #21000"),
  data.table(field.id=21000, var="ethnicity_subgroup", name="Ethnicity")
)

# Reorder rows and columns
info <- info[names(raw), on = .(var), nomatch=0]
info <- info[, .(var, name, field.id)]

# overwrite birth country entry
info[var == "birth_country", name := "Country of birth (combines fields #1647 and #20115)"]
info[var == "birth_country", field.id := NA]

# Write out
fwrite(raw, quote=FALSE, file="Demographics/demographics.csv")
fwrite(info, file="Demographics/column_information.csv")

# Upload to persistent storage
system("dx upload Demographics/demographics.csv Demographics/column_information.csv --destination 'common/Demographics/'", wait=TRUE)

# Send raw data to deletion folder to reduce storage costs - this needs to be 
# done in two steps, since we can't move two folders with the same name to the
# same location, so here we rename with a random number and then move to trash
rn <- as.integer(Sys.time())
system(sprintf("dx mv 'common/Demographics/raw_data/' 'common/Demographics/raw_data_%s'", rn), wait=TRUE)
system(sprintf("dx mv 'common/Demographics/raw_data_%s/' trash/", rn), wait=TRUE) 
