library(data.table)
library(lubridate)

# Make output directory
system("mkdir -p followup")

# Pull down raw data that has been extracted with Table Exporter
system("mkdir -p raw_data", wait=TRUE)
system("dx download 'common/Follow-up/raw_data/data.csv' -o raw_data/followup.csv", wait=TRUE)

# Pull down other data sources we need for this script
system("dx download common/Demographics/demographics.csv")
system("dx download 'common/Hospital Records/hospital_records.csv'")
system("dx download common/Deaths/deaths.csv")
system("dx download 'common/Primary Care/gp_registrations.csv'")
system("dx download 'common/Primary Care/gp_clinical_records.csv'")
system("dx download 'common/Cancer Register/cancer_register'.csv")

#################
# Curate follow-up available in each record source
#################

# Load hospital records and compute follow-up available per nation - they each
# have different cutoff and start dates
hes <- fread("hospital_records.csv", na.strings=c("", "NA"))
hes[, min_date := pmin(episode_start, episode_end, admission_date, na.rm=TRUE)] # decision_date sometimes decades before linkage!
hes[, max_date := pmax(episode_start, episode_end, admission_date, decision_date, discharge_date, na.rm=TRUE)]
hes_follow <- hes[, .(min_follow = min(min_date, na.rm=TRUE), max_follow = max(max_date, na.rm=TRUE)), by=hospital_nation]

# Load death records and compute follow-up available per nation - they each have
# different cutoff dates
deaths <- fread("deaths.csv")
death_follow <- deaths[,.(max_follow=max(date_of_death, na.rm=TRUE)), by=death_nation]

# Load cancer registry records and compute follow-up available per data source
cancer <- fread("cancer_register.csv", na.strings=c("", "NA"))
cancer_follow <- cancer[!is.na(diagnosis_date), .(min_follow=min(diagnosis_date), max_follow=max(diagnosis_date)), by=record_origin]
cancer_follow <- rbind(idcol="cancer_nation",
  "England"=cancer_follow[record_origin == "Originating from England/Wales"],
  "Wales"=cancer_follow[record_origin == "Originating from England/Wales"],
  "Scotland"=cancer_follow[record_origin == "Originating from Scotland"]
)
cancer_follow <- cancer_follow[, .(cancer_nation, min_follow, max_follow)]

# Load in GP registration records and compute follow-up available
# N.b., Minimum follow-up available per person in the registration records
gp_records <- fread("gp_clinical_records.csv", na.strings=c("", "NA"))
gp_follow <- gp_records[!is.na(event_date), .(max_follow=max(event_date)), by=data_provider]

######
# Hard code UKB-recommended censor dates determined by when they judge records
# for a month to be largely complete for a given data source
# https://biobank.ndph.ox.ac.uk/showcase/exinfo.cgi?src=Data_providers_and_dates
######
hes_follow[, min_censor := as.IDate(fcase(
  hospital_nation == "England", "1997-01-01",
  hospital_nation == "Wales", "1991-01-01",
  hospital_nation == "Scotland", "1981-01-01"
))]

hes_follow[, icd9to10_switch := as.IDate(fcase(
  hospital_nation == "England", "1997-01-01",
  hospital_nation == "Wales", "1999-01-01",
  hospital_nation == "Scotland", "1996-01-01"
))]

hes_follow[, opcs3to4_switch := as.IDate(fcase(
  hospital_nation == "England", "1997-01-01",
  hospital_nation == "Wales", "1999-01-01",
  hospital_nation == "Scotland", "1989-01-01"
))]

hes_follow[, max_censor := as.IDate(fcase(
  hospital_nation == "England", "2023-03-31",
  hospital_nation == "Wales", "2022-05-31",
  hospital_nation == "Scotland", "2022-08-31"
))]

death_follow[, max_censor := as.IDate(fcase(
  death_nation == "Scotland", "2024-11-30",
  death_nation == "England/Wales", "2024-08-31"
))]

cancer_follow[, min_censor := as.IDate(fcase(
  cancer_nation == "England", "1971-01-01",
  cancer_nation == "Wales", "1971-01-01",
  cancer_nation == "Scotland", "1957-01-01"
))]

cancer_follow[, icd9to10_switch := as.IDate(fcase(
  cancer_nation == "England", "1995-01-01",
  cancer_nation == "Wales", "1995-01-01",
  cancer_nation == "Scotland", "1997-01-01"
))]

cancer_follow[, max_censor := as.IDate(fcase(
  cancer_nation == "England", "2023-05-31",
  cancer_nation == "Wales", "2016-12-31",
  cancer_nation == "Scotland", "2023-09-30"
))]

gp_follow[, max_censor := as.IDate(fcase(
  data_provider == "England (TPP)", "2016-05-31",
  data_provider == "England (Vision)", "2017-05-31",
  data_provider == "Wales", "2017-08-31",
  data_provider == "Scotland", "2017-03-31"
))]

###########
# Output available follow-up information
###########
o1 <- options()[["datatable.print.class"]]
o2 <- options()[["datatable.print.rownames"]]
options(datatable.print.class=FALSE)
options(datatable.print.rownames=FALSE)
sink("followup/available_followup.txt")
cat("Death records:\n")
cat("---------------\n")
print(death_follow[])
cat("\nHospital records:\n")
cat("------------------\n")
print(hes_follow[])
cat("\n")
cat("\nCancer registry:\n")
cat("-----------------\n")
print(cancer_follow[])
cat("\n")
cat("\nPrimary care records:\n")
cat("----------------------\n")
print(gp_follow[])
cat("\n")
sink()
options(datatable.print.class=o1)
options(datatable.print.rownames=o2)

#####
# Additional data curation required to determine participant's location during 
# followup
#####

# Load in additional information we will need to curate follow-up information
demo <- fread("demographics.csv")
gp_rego <- fread("gp_registrations.csv")

# Curate information about loss of follow-up and immigration to the UK
# Load and information regarding lost to follow-up information
follow <- fread("raw_data/followup.csv", na.strings=c("", "NA"))

lost <- follow[, .(eid, lost_to_followup_reason=p190, lost_to_followup_date=p191)]
lost <- lost[!is.na(lost_to_followup_reason) | !is.na(lost_to_followup_date)]
lost[, lost_to_followup_reason := fcase(
  lost_to_followup_reason == 1, "Death reported to UK Biobank by a relative",
  lost_to_followup_reason == 2, "NHS records indicate they are lost to follow-up",
  lost_to_followup_reason == 3, "NHS records indicate they have left the UK",
  lost_to_followup_reason == 4, "UK Biobank sources report they have left the UK",
  lost_to_followup_reason == 5, "Participant has withdrawn consent for future linkage"
)]

lost_left_uk <- lost[lost_to_followup_reason %in% c("NHS records indicate they have left the UK", "UK Biobank sources report they have left the UK")]
lost_death <- lost[lost_to_followup_reason == "Death reported to UK Biobank by a relative"]
lost_other <- lost[lost_to_followup_reason == "NHS records indicate they are lost to follow-up"]
withdrawn_consent <- lost[lost_to_followup_reason == "Participant has withdrawn consent for future linkage"]

immigrated <- rbind(idcol="visit_index",
  "0"=follow[, .(eid, year_immigrated=p3659_i0)],
  "1"=follow[, .(eid, year_immigrated=p3659_i1)],
  "2"=follow[, .(eid, year_immigrated=p3659_i2)]
)
immigrated[, visit_index := as.integer(visit_index)]
immigrated <- immigrated[!is.na(year_immigrated)]
immigrated[year_immigrated == -1, year_immigrated := NA] # Do not know 
immigrated[year_immigrated == -3, year_immigrated := NA] # Prefer not to answer
immigrated[!is.na(year_immigrated), date_immigrated := as.IDate(sprintf("%s-07-01", year_immigrated))] # set as midpoint of year for now

# Curate cancer registry events so we can work out location changes
cancer_nation <- cancer[record_origin %like% "England/Wales" | record_origin %like% "Scotland"]
cancer_nation[, record_origin := gsub(".* ", "", record_origin)]

# Split out hospital events into ones where we know the person is attending in 
# their home location, unsure, or definitely not
hes_home <- hes[
  (administrative_body == gp_authority) | # GP is located in the same area as hospital
    admission_source %like% "Usual Place of residence" |
    admission_source %like% "Local Authority" |
    admission_source %like% "Nursing home" |
    admission_source %like% "Hospice" |
    admission_source %like% "Hospital at Home" |
    discharge_destination %like% "Usual Place of residence" |
    discharge_destination %like% "Local Authority" |
    discharge_destination %like% "Nursing home" |
    discharge_destination %like% "Hospice" |
    discharge_destination %like% "Hospital at Home" 
]
hes_away <- hes[!(dnx_hesin_id %in% hes_home$dnx_hesin_id) & (
  administrative_body != gp_authority | # Usual GP is located in a different administrative area
  admission_source %like% "Holiday accomodation" |
  discharge_destination %like% "Holiday accomodation"
)]
hes_other <- hes[!(dnx_hesin_id %in% hes_away$dnx_hesin_id) & !(dnx_hesin_id %in% hes_home$dnx_hesin_id)]

################################
# Construct a sequence of events for each person so we can work out what their
# minimimum and maximum follow-up is in each record source
################################
events <- demo[, .(eid, date=assessment_date, type="UKB assessment", location=assessment_nation)]
events <- rbind(events, demo[,.(eid, date=approx_birth_date, type="Inferred DOB", location=NA)]) # Duplicates/disagreements may need resolving later
events <- rbind(events, immigrated[!is.na(date_immigrated), .(eid, date=date_immigrated, type="Immigration to UK", location=NA)]) # Duplicates/disagreements may need resolving later
events <- rbind(events, lost_left_uk[!is.na(lost_to_followup_date), .(eid, date=lost_to_followup_date, type="Lost to followup (left UK)", location=NA)]) # Duplicates/disagreements may need resolving later
events <- rbind(events, lost_death[!is.na(lost_to_followup_date), .(eid, date=lost_to_followup_date, type="Death (reported by relative)", location=NA)]) # Duplicates/disagreements may need resolving later
events <- rbind(events, lost_other[!is.na(lost_to_followup_date), .(eid, date=lost_to_followup_date, type="Lost to followup (other reason)", location=NA)]) # Duplicates/disagreements may need resolving later
events <- rbind(events, gp_rego[!is.na(registration_date), .(eid, date=registration_date, type="GP registration", location=data_provider)])
events <- rbind(events, gp_rego[!is.na(deregistration_date), .(eid, date=deregistration_date, type="GP deregistration", location=data_provider)])
events <- rbind(events, gp_records[!is.na(data_provider) & !is.na(event_date), .(eid, date=event_date, type="GP record", location=data_provider)])
events <- rbind(events, cancer_nation[!is.na(diagnosis_date), .(eid, date=diagnosis_date, type="Cancer registry", location=record_origin)])
events <- rbind(events, cancer[!is.na(diagnosis_date) & !(record_origin %like% "Originating"), .(eid, date=diagnosis_date, type="Cancer registry", location=NA)]) # Missing location data
events <- rbind(events, deaths[, .(eid, date=date_of_death, type="Death registry", location=death_nation)])
events <- rbind(events, hes_home[!is.na(min_date), .(eid, date=min_date, type="Hospital record", location=hospital_nation)])
events <- rbind(events, hes_away[!is.na(min_date), .(eid, date=min_date, type="Hospital record (unusual location)", location=hospital_nation)])
events <- rbind(events, hes_other[!is.na(min_date), .(eid, date=min_date, type="Hospital record (location uncertain)", location=hospital_nation)])
events <- events[order(date)][order(eid)]
events <- unique(events) # strip duplicates by date

#########
# Adjudicate conflicts based on sequence of events
#########

# Extract date of birth, and resolve conflicts (if any) if inferred date of birth
# differs across UK Biobank assessments
dob <- events[type == "Inferred DOB"] 
stopifnot( unique(dob)[,.N,by=eid][N > 1, .N] == 0 ) # Thankfully no disagreements across assessment visits

# Adjudicate immigration date - these have been inferred as the midpoint of the
# year of immigration reported, but these self reported years of immigration
# sometimes differ wildly across UKB assessments, and there may be other records
# or events during those same years that will allow us to do better than taking
# the midpoint of the year (which in two cases, would result in an immigration
# date later than UKB assessment!)
immigrated <- unique(events[type == "Immigration to UK"])
setnames(immigrated, "date", "immigration_date")
immigrated[, immigration_date_copy := immigration_date] # preserved column in result on non-equi join

# Immigration date has been inferred as the midpoint of the self-reported year
# immigrated; what HES/GP/Cancer records exist (if any) prior to this?
before_immigration <- immigrated[
  events[type != "Inferred DOB" & type != "Immigration to UK"], 
  on = .(eid, immigration_date>date), nomatch=0,
  # non-equi joins replace the original column with contents of the matched column
  .(eid, immigration_date=immigration_date_copy, event_date=immigration_date, type=i.type, location=i.location) 
]
immigrated[, immigration_date_copy := NULL]

# How many?
immigrated[, earlier_records := 0]
immigrated[before_immigration[,.N,by=.(eid, immigration_date)], on = .(eid, immigration_date), earlier_records := i.N]

# For each set of possible immigration dates, take the latest date (if multiple)
# among the set with the smallest number of records prior to the inferred 
# immigration date
immigrated <- immigrated[order(-immigration_date)][order(earlier_records)][order(eid)]
immigrated <- immigrated[,.SD[1], by=eid]
before_immigration <- before_immigration[immigrated[,.(eid, immigration_date)], on = .(eid, immigration_date), nomatch=0]

# Re-infer immigration date as the midpoint between the start of the year and
# the first event in that same year (or, if the same year as date of birth, 
# the midpoint between the inferred date of birth and the end of the year/next 
# record in same year)
immigrated[, year_immigrated := year(immigration_date)]
events[, year := year(date)]
same_year <- events[immigrated, on = .(eid, year = year_immigrated), nomatch=0]
same_year <- same_year[type != "Immigration to UK"]

reinf <- unique(same_year[,.(eid, year)])
reinf[, earliest := as.IDate(sprintf("%s-01-01", year))]
reinf[, latest := as.IDate(sprintf("%s-12-31", year))]
reinf[same_year[type == "Inferred DOB"], on = .(eid), earliest := i.date]
reinf[same_year[type != "Inferred DOB"][order(date)][,.SD[1],by=eid], on = .(eid), latest := i.date]
reinf[, midpoint := as.IDate((as.numeric(latest) - as.numeric(earliest))/2 + as.numeric(earliest))]
immigrated[reinf, on = .(eid), immigration_date := i.midpoint]

# Strip out any events that occur prior to the year of self-reported immigration
before_immigration <- events[immigrated, on = .(eid, year < year_immigrated), nomatch=0,
                             .(eid, immigration_date, event_date=date, type, location)]
before_immigration <- before_immigration[type != "Inferred DOB"] # inferred from age, not necessarily in UK
before_immigration <- before_immigration[type != "Immigration to UK"] # Conflicting reports across UKB assessments adjudicated above
events <- events[!before_immigration, on = .(eid, date=event_date)]
events[, year := NULL] # no longer needed

# Strip out any events that occur after loss of follow-up
lost <- events[type %like% "Lost" | type == "Death (reported by relative)"]
setnames(lost, "date", "date_lost")
after_lost <- events[lost, on = .(eid, date > date_lost), nomatch=0, .(eid, date=x.date, type, location)]
events <- events[!after_lost, on = .(eid, date)]

##########
# Make follow-up table
##########

# Begin constructing follow-up table for each participant
out <- unique(demo[,.(eid)])

# Drop people who have withdrawn consent for EHR linkage
out <- out[!withdrawn_consent, on = .(eid)]

# Add information about presence/absence of each participant in various records
out[, has_died := FALSE]
out[events[type %like% "Death"], on = .(eid), has_died := TRUE]

out[, any_hospitalisations := FALSE]
out[events[type %like% "Hospital"], on = .(eid), any_hospitalisations := TRUE]

out[, in_cancer_register := FALSE]
out[events[type %like% "Cancer"], on = .(eid), in_cancer_register := TRUE]

out[, linked_gp_records := FALSE]
out[events[type == "GP record"], on = .(eid), linked_gp_records := TRUE]

out[demo[visit_index == 0], on = .(eid), baseline_location := assessment_nation]

# Get the set of events we will use to determine location of residence at the 
# time of each possible censor date
#
# First we need to filter the events table to only those where the location is 
# known and we are pretty sure the event record was issued in the same location 
# the person usually resided at that time
events_filtered <- events[!is.na(location) & 
  type != "Hospital record (location uncertain)" & 
  type != "Hospital record (unusual location)" &
  location != "England/Wales" # Cancer registry records with insufficient information for curating follow-up
]

#########
# Determine minimum follow-up time per person in the linked hospital records
#########

# For each person, create a counter-factual for each possible nation to find the
# closest record in time to the minimum censor date
cf <- rbind(idcol="hospital_nation",
  "England"=cbind(events_filtered, "min_censor"=hes_follow[hospital_nation == "England", min_censor]),
  "Wales"=cbind(events_filtered, "min_censor"=hes_follow[hospital_nation == "Wales", min_censor]),
  "Scotland"=cbind(events_filtered, "min_censor"=hes_follow[hospital_nation == "Scotland", min_censor])
)
cf[location == "England (Vision)", location := "England"]
cf[location == "England (TPP)", location := "England"]
cf[, dist := abs(date - min_censor)]
cf <- cf[order(dist)][order(hospital_nation)][order(eid)]

# Among each of the three possible minimum censor dates, select the one where
# the adjacent record is the closest in time for each participant, then amongst
# the three possible minimum censor dates, choose the one with the earliest 
# matching record
closest <- cf[hospital_nation == location, .SD[which.min(dist)], by=.(eid, hospital_nation)]
earliest <- closest[, .SD[which.min(date)], by=eid]

# Update the follow-up table
out[earliest, on = .(eid), c("min_hospital_censor", "min_hospital_location", "min_hospital_reason") := 
  .(i.min_censor, i.location, i.type)]

# Roll forward the follow-up for anyone born after the hospital records started
out[dob, on = .(eid, min_hospital_censor < date), 
  c("min_hospital_censor", "min_hospital_location", "min_hospital_reason") := .
  (i.date, "Birth", "Birth")]

# Roll forward the follow-up for anyone who immigrated to the UK after the start
# of the hospital records
out[immigrated, on = .(eid, min_hospital_censor < immigration_date), 
  c("min_hospital_censor", "min_hospital_location", "min_hospital_reason") := .
  (immigration_date, "Immigrated", "Immigrated")]

#########
# Determine minimum follow-up time per person for cancer registry followup
#########

# For each person, create a counter-factual for each possible nation to find the
# closest record in time to the minimum censor date
cf <- rbind(idcol="cancer_nation",
  "England"=cbind(events_filtered, "min_censor"=cancer_follow[cancer_nation == "England", min_censor]),
  "Wales"=cbind(events_filtered, "min_censor"=cancer_follow[cancer_nation == "Wales", min_censor]),
  "Scotland"=cbind(events_filtered, "min_censor"=cancer_follow[cancer_nation == "Scotland", min_censor])
)
cf[location == "England (Vision)", location := "England"]
cf[location == "England (TPP)", location := "England"]
cf[, dist := abs(date - min_censor)]
cf <- cf[order(dist)][order(cancer_nation)][order(eid)]

# Among each of the three possible minimum censor dates, select the one where
# the adjacent record is the closest in time for each participant, then amongst
# the three possible minimum censor dates, choose the one with the earliest 
# matching record
closest <- cf[cancer_nation == location, .SD[which.min(dist)], by=.(eid, cancer_nation)]
earliest <- closest[, .SD[which.min(date)], by=eid]

# Update the follow-up table
out[earliest, on = .(eid), c("min_cancer_censor", "min_cancer_location", "min_cancer_reason") := 
  .(i.min_censor, i.location, i.type)]

# Roll forward the follow-up for anyone born after the hospital records started
out[dob, on = .(eid, min_cancer_censor < date), 
  c("min_cancer_censor", "min_cancer_location", "min_cancer_reason") := .
  (i.date, "Birth", "Birth")]

# Roll forward the follow-up for anyone who immigrated to the UK after the start
# of the hospital records
out[immigrated, on = .(eid, min_cancer_censor < immigration_date), 
  c("min_cancer_censor", "min_cancer_location", "min_cancer_reason") := .
  (immigration_date, "Immigrated", "Immigrated")]

#########
# Determine minimum follow-up time per person in GP records
#########

# Find the earliest GP record available for each person
# N.b. a handful of people have an earlier GP deregistration record, but we will
# ignore this as they have no records from the GP they deregistered from
earliest_gp <- events[type %in% c("GP record", "GP registration"), .SD[1], by=eid]
earliest_gp[location == "England (Vision)", location := "England"]
earliest_gp[location == "England (TPP)", location := "England"]

# Update the follow-up table - note 17 people with GP records without censoring
# information, as all their GP records had no date attached
out[earliest_gp, on = .(eid), c("min_gp_censor", "min_gp_location", "min_gp_reason") := 
  .(i.date, i.location, i.type)]

# Some people reporting immigration to the UK also have linked GP records (years) 
# prior to these dates - assume these are either linkage errors or 
# discontinunities in follow-up (i.e. left the UK, then re-immigrated) 
out[immigrated, on = .(eid, min_gp_censor <= immigration_date), 
    c("min_gp_censor", "min_gp_location", "min_gp_reason") := .
    (immigration_date, "Immigrated", "GP records prior to immigration")]

#########
# Determine minimum censor date for ICD-10 only analysis (or equivalently, 
# maximum censor date for ICD-9 only analysis) in the linked hospital records
#########

# For each person, create a counter-factual for each possible nation to find the
# closest record in time to the switch date
cf <- rbind(idcol="hospital_nation",
  "England"=cbind(events_filtered, "min_censor"=hes_follow[hospital_nation == "England", icd9to10_switch]),
  "Wales"=cbind(events_filtered, "min_censor"=hes_follow[hospital_nation == "Wales", icd9to10_switch]),
  "Scotland"=cbind(events_filtered, "min_censor"=hes_follow[hospital_nation == "Scotland", icd9to10_switch])
)
cf[location == "England (Vision)", location := "England"]
cf[location == "England (TPP)", location := "England"]
cf[, dist := abs(date - min_censor)]
cf <- cf[order(dist)][order(hospital_nation)][order(eid)]

# Among each of the three possible minimum censor dates, select the one where
# the adjacent record is the closest in time for each participant, then amongst
# the three possible minimum censor dates, choose the one with the earliest 
# matching record
closest <- cf[hospital_nation == location, .SD[which.min(dist)], by=.(eid, hospital_nation)]
earliest <- closest[, .SD[which.min(date)], by=eid]

# Update the follow-up table
out[earliest, on = .(eid), c("hospital_icd9to10_switch", "hospital_icd9to10_location", "hospital_icd9to10_reason") := 
      .(i.min_censor, i.location, i.type)]

# Roll forward the follow-up for anyone born after the switch date
out[dob, on = .(eid, hospital_icd9to10_switch < date), 
    c("hospital_icd9to10_switch", "hospital_icd9to10_location", "hospital_icd9to10_reason") := .
    (i.date, "Birth", "Birth")]

# Roll forward the follow-up for anyone who immigrated to the UK after the 
# switch date
out[immigrated, on = .(eid, hospital_icd9to10_switch < immigration_date), 
    c("hospital_icd9to10_switch", "hospital_icd9to10_location", "hospital_icd9to10_reason") := .
    (immigration_date, "Immigrated", "Immigrated")]

#########
# Determine minimum censor date for OPCS-4 only analysis (or equivalently, 
# maximum censor date for OPCS-3 only analysis) in the linked hospital records
#########

# For each person, create a counter-factual for each possible nation to find the
# closest record in time to the switch date
cf <- rbind(idcol="hospital_nation",
  "England"=cbind(events_filtered, "min_censor"=hes_follow[hospital_nation == "England", opcs3to4_switch]),
  "Wales"=cbind(events_filtered, "min_censor"=hes_follow[hospital_nation == "Wales", opcs3to4_switch]),
  "Scotland"=cbind(events_filtered, "min_censor"=hes_follow[hospital_nation == "Scotland", opcs3to4_switch])
)
cf[location == "England (Vision)", location := "England"]
cf[location == "England (TPP)", location := "England"]
cf[, dist := abs(date - min_censor)]
cf <- cf[order(dist)][order(hospital_nation)][order(eid)]

# Among each of the three possible minimum censor dates, select the one where
# the adjacent record is the closest in time for each participant, then amongst
# the three possible minimum censor dates, choose the one with the earliest 
# matching record
closest <- cf[hospital_nation == location, .SD[which.min(dist)], by=.(eid, hospital_nation)]
earliest <- closest[, .SD[which.min(date)], by=eid]

# Update the follow-up table
out[earliest, on = .(eid), c("hospital_opcs3to4_switch", "hospital_opcs3to4_location", "hospital_opcs3to4_reason") := 
      .(i.min_censor, i.location, i.type)]

# Roll forward the follow-up for anyone born after the switch date
out[dob, on = .(eid, hospital_opcs3to4_switch < date), 
    c("hospital_opcs3to4_switch", "hospital_opcs3to4_location", "hospital_opcs3to4_reason") := .
    (i.date, "Birth", "Birth")]

# Roll forward the follow-up for anyone who immigrated to the UK after the 
# switch date
out[immigrated, on = .(eid, hospital_opcs3to4_switch < immigration_date), 
    c("hospital_opcs3to4_switch", "hospital_opcs3to4_location", "hospital_opcs3to4_reason") := .
    (immigration_date, "Immigrated", "Immigrated")]

#########
# Determine minimum censor date for ICD-10 only analysis (or equivalently, 
# maximum censor date for ICD-9 only analysis) in the cancer register
#########

# For each person, create a counter-factual for each possible nation to find the
# closest record in time to the minimum censor date
cf <- rbind(idcol="cancer_nation",
  "England"=cbind(events_filtered, "min_censor"=cancer_follow[cancer_nation == "England", icd9to10_switch]),
  "Wales"=cbind(events_filtered, "min_censor"=cancer_follow[cancer_nation == "Wales", icd9to10_switch]),
  "Scotland"=cbind(events_filtered, "min_censor"=cancer_follow[cancer_nation == "Scotland", icd9to10_switch])
)
cf[location == "England (Vision)", location := "England"]
cf[location == "England (TPP)", location := "England"]
cf[, dist := abs(date - min_censor)]
cf <- cf[order(dist)][order(cancer_nation)][order(eid)]

# Among each of the three possible minimum censor dates, select the one where
# the adjacent record is the closest in time for each participant, then amongst
# the three possible minimum censor dates, choose the one with the earliest 
# matching record
closest <- cf[cancer_nation == location, .SD[which.min(dist)], by=.(eid, cancer_nation)]
earliest <- closest[, .SD[which.min(date)], by=eid]

# Update the follow-up table
out[earliest, on = .(eid), c("cancer_icd9to10_switch", "cancer_icd9to10_location", "cancer_icd9to10_reason") := 
      .(i.min_censor, i.location, i.type)]

# Roll forward the follow-up for anyone born after the hospital records started
out[dob, on = .(eid, cancer_icd9to10_switch < date), 
    c("cancer_icd9to10_switch", "cancer_icd9to10_location", "cancer_icd9to10_reason") := .
    (i.date, "Birth", "Birth")]

# Roll forward the follow-up for anyone who immigrated to the UK after the start
# of the hospital records
out[immigrated, on = .(eid, cancer_icd9to10_switch < immigration_date), 
    c("cancer_icd9to10_switch", "cancer_icd9to10_location", "cancer_icd9to10_reason") := .
    (immigration_date, "Immigrated", "Immigrated")]

#########
# Determine maximum follow-up time in the death registry
#########

# For each person, create a counter-factual for each possible nation to find the
# closest record in time to the maximum censor date
cf <- rbind(idcol="death_nation",
  "England/Wales"=cbind(events_filtered, "max_censor"=death_follow[death_nation == "England/Wales", max_censor]),
  "Scotland"=cbind(events_filtered, "max_censor"=death_follow[death_nation == "Scotland", max_censor])
)
cf[location == "England (Vision)", location := "England"]
cf[location == "England (TPP)", location := "England"]
cf[, dist := abs(date - max_censor)]
cf <- cf[order(dist)][order(death_nation)][order(eid)]

# Among each of the three possible minimum censor dates, select the one where
# the adjacent record is the closest in time for each participant, then amongst
# the three possible minimum censor dates, choose the one with the latest
# matching record
closest <- cf[death_nation == location | (death_nation == "England/Wales" & location %in% c("England", "Wales")),
 .SD[which.min(dist)], by=.(eid, death_nation)]
latest <- closest[, .SD[which.max(date)], by=eid]

# Update the follow-up table
out[latest, on = .(eid), c("max_death_censor", "max_death_location", "max_death_reason") := 
  .(i.max_censor, i.location, i.type)]

# Roll back the max follow-up where a death has occurred
deaths <- events[type %like% "Death"]
deaths <- deaths[,.SD[1], by=eid]
out[deaths, on = .(eid), c("max_death_censor", "max_death_location", "max_death_reason") := 
 .(i.date, "Death", i.type)]

# Roll back the max follow-up where the person has been lost to follow-up
lost <- events[type %like% "Lost to followup"]
lost <- lost[,.SD[1], by=eid]
out[lost, on = .(eid, max_death_censor>date), 
  c("max_death_censor", "max_death_location", "max_death_reason") := 
  .(i.date, "Lost", i.type)]

#########
# Determine maximum follow-up time in the hospital follow-up
#########

# For each person, create a counter-factual for each possible nation to find the
# closest record in time to the maximum censor date
cf <- rbind(idcol="hospital_nation",
  "England"=cbind(events_filtered, "max_censor"=hes_follow[hospital_nation == "England", max_censor]),
  "Wales"=cbind(events_filtered, "max_censor"=hes_follow[hospital_nation == "Wales", max_censor]),
  "Scotland"=cbind(events_filtered, "max_censor"=hes_follow[hospital_nation == "Scotland", max_censor])
)
cf[location == "England (Vision)", location := "England"]
cf[location == "England (TPP)", location := "England"]
cf[, dist := abs(date - max_censor)]
cf <- cf[order(dist)][order(hospital_nation)][order(eid)]

# Among each of the three possible minimum censor dates, select the one where
# the adjacent record is the closest in time for each participant, then amongst
# the three possible minimum censor dates, choose the one with the latest
# matching record
closest <- cf[hospital_nation == location, .SD[which.min(dist)], by=.(eid, hospital_nation)]
latest <- closest[, .SD[which.max(date)], by=eid]

# Update the follow-up table
out[latest, on = .(eid), c("max_hospital_censor", "max_hospital_location", "max_hospital_reason") := 
  .(i.max_censor, i.location, i.type)]

# Roll back the max follow-up where a death has occurred before the censor date
out[deaths, on = .(eid, max_hospital_censor>date), 
  c("max_hospital_censor", "max_hospital_location", "max_hospital_reason") := 
  .(i.date, "Death", i.type)]

# Roll back the max follow-up where the person has been lost to follow-up
out[lost, on = .(eid, max_hospital_censor>date), 
  c("max_hospital_censor", "max_hospital_location", "max_hospital_reason") := 
  .(i.date, "Lost", i.type)]

#########
# Determine maximum follow-up time in the cancer registries
#########

# For each person, create a counter-factual for each possible nation to find the
# closest record in time to the maximum censor date
cf <- rbind(idcol="cancer_nation",
  "England"=cbind(events_filtered, "max_censor"=cancer_follow[cancer_nation == "England", max_censor]),
  "Wales"=cbind(events_filtered, "max_censor"=cancer_follow[cancer_nation == "Wales", max_censor]),
  "Scotland"=cbind(events_filtered, "max_censor"=cancer_follow[cancer_nation == "Scotland", max_censor])
)
cf[location == "England (Vision)", location := "England"]
cf[location == "England (TPP)", location := "England"]
cf[, dist := abs(date - max_censor)]
cf <- cf[order(dist)][order(cancer_nation)][order(eid)]

# Among each of the three possible minimum censor dates, select the one where
# the adjacent record is the closest in time for each participant, then amongst
# the three possible minimum censor dates, choose the one with the latest
# matching record
closest <- cf[cancer_nation == location, .SD[which.min(dist)], by=.(eid, cancer_nation)]
latest <- closest[, .SD[which.max(date)], by=eid]

# Update the follow-up table
out[latest, on = .(eid), c("max_cancer_censor", "max_cancer_location", "max_cancer_reason") := 
  .(i.max_censor, i.location, i.type)]

# Roll back the max follow-up where a death has occurred before the censor date
out[deaths, on = .(eid, max_cancer_censor>date), 
  c("max_cancer_censor", "max_cancer_location", "max_cancer_reason") := 
  .(i.date, "Death", i.type)]

# Roll back the max follow-up where the person has been lost to follow-up
out[lost, on = .(eid, max_cancer_censor>date), 
  c("max_cancer_censor", "max_cancer_location", "max_cancer_reason") := 
  .(i.date, "Lost", i.type)]

#########
# Determine maximum follow-up time per person in GP records
#########

# Find the latest GP record available for each person
latest_gp <- events[type %like% "GP", .SD[.N], by=eid]

# Split out events where the last record is a de-registration, these will be
# taken as the max follow-up for that person
derego <- latest_gp[type == "GP deregistration"]
latest_gp <- latest_gp[type != "GP deregistration"]

# If the latest record isn't a GP de-registration, set the maximum follow-up
# based on the record source
latest_gp[gp_follow, on = .(location=data_provider), date := max_censor]

# Recombine
latest_gp <- rbind(derego, latest_gp)

# Recode reason and location
latest_gp[type == "GP deregistration" & location == "England (TPP)", location := "England"]
latest_gp[type == "GP deregistration" & location == "England (Vision)", location := "England"]
latest_gp[location == "England (TPP)", c("location", "type") := .("England", paste(type, "(TPP)"))]
latest_gp[location == "England (Vision)", c("location", "type") := .("England", paste(type, "(Vision))"))]

# Recombine and add to output table
out[latest_gp, on = .(eid), c("max_gp_censor", "max_gp_location", "max_gp_reason") := 
  .(i.date, i.location, i.type)]

# Roll back the max follow-up where a death has occurred before the censor date
out[deaths, on = .(eid, max_gp_censor>date), 
  c("max_gp_censor", "max_gp_location", "max_gp_reason") := 
  .(i.date, "Death", i.type)]

# Roll back the max follow-up where the person has been lost to follow-up
out[lost, on = .(eid, max_gp_censor>date), 
  c("max_gp_censor", "max_gp_location", "max_gp_reason") := 
  .(i.date, "Lost", i.type)]

############################
# Curate column information
############################

# Curate information about each column:
info <- rbind(
  data.table(var="eid", name="Application-specific participant ID for participants consenting to electronic health record linkage"),
  data.table(var="has_died", name="TRUE where the participant has died since baseline assessment"),
  data.table(var="any_hospitalisations", name="TRUE where any hospital records exist for the participant"),
  data.table(var="in_cancer_register", name="TRUE where participant has any records in the cancer registry data"),
  data.table(var="linked_gp_records", name="TRUE where participant has linkage to primary care data"),
  data.table(var="baseline_nation", name="Nation within the UK the participant was located in at baseline assessment"),
  data.table(var="min_hospital_censor", name="Earliest date for which the participant has hospital record linkage based on date of earliest complete hospital record linkage available in each nation, inferred date of birth, or inferred date of immigration to the UK"),
  data.table(var="min_hospital_location", name="Inferred nation of residence used to determine earliest hospital record linkage, or \"Birth\" or \"Immigrated\" where these occurred after the start of hospital record linkage in the participant's earliest inferred nation of residence"),
  data.table(var="min_hospital_reason", name="Record source used to infer the participant's earliest nation of residence or otherwise set the earliest date of hospital record linkage"),
  data.table(var="min_cancer_censor", name="Earliest date for which the participant could possibly have cancer registry data based on date of earliest complete cancer registry linkage available in each nation, inferred date of birth, or inferred date of immigration to the UK"),
  data.table(var="min_cancer_location", name="Inferred nation of residence used to determine earliest possible cancer registry linkage, or \"Birth\" or \"Immigrated\" where these occurred after the start of hospital record linkage in the participants earliest inferred nation of residence"),
  data.table(var="min_cancer_reason", name="Record source used to infer the participant's earliest nation of residence or otherwise set the earliest date of cancer registry linkage"),
  data.table(var="min_gp_censor", name="Earliest date for which primary care record linkage is available for this participant, if any"),
  data.table(var="min_gp_location", name="Nation of residence of the participant at the earliest available primary care record, if any"),
  data.table(var="min_gp_reason", name="Type of first available primary care record available for this participant, if any"),
  data.table(var="hospital_icd9to10_switch", name="Minimum censor date for ICD-10 only analysis/maximum censor date for ICD-9 only analyses in the hospital records for this participant"),
  data.table(var="hospital_icd9to10_location", name="Inferred nation of residence used to determine the ICD-9 to ICD-10 switch over date in the hospital records, \"Birth\", or \"Immigrated\""),
  data.table(var="hospital_icd9to10_reason", name="Record source used to infer the participant's nation of residence when determine the ICD-9 to ICD-10 switch over date in the hospital records"),
  data.table(var="hospital_opcs3to4_switch", name="Minimum censor date for OPCS-4 only analysis/maximum censor date for OPCS-3 only analyses in the hospital records for this participant"),
  data.table(var="hospital_opcs3to4_location", name="Inferred nation of residence used to determine the OPCS-3 to OPCS-4 switch over date in the hospital records, \"Birth\", or \"Immigrated\""),
  data.table(var="hospital_opcs3to4_reason", name="Record source used to infer the participant's nation of residence when determine the OPCS-3 to OPCS-4 switch over date in the hospital records"),
  data.table(var="cancer_icd9to10_switch", name="Minimum censor date for ICD-10 only analysis/maximum censor date for ICD-9 only analyses in the cancer register for this participant"),
  data.table(var="cancer_icd9to10_location", name="Inferred nation of residence used to determine the ICD-9 to ICD-10 switch over date in the cancer register, \"Birth\", or \"Immigrated\""),
  data.table(var="cancer_icd9to10_reason", name="Record source used to infer the participant's nation of residence when determine the ICD-9 to ICD-10 switch over date in the cancer register"),
  data.table(var="max_death_censor", name="Latest date for which the participant has death registry record linkage based on date of latest complete death record linkage available in each nation, date of death, or loss of follow-up"),
  data.table(var="max_death_location", name="Inferred nation of residence used to determine latest death registry record linkage, or \"Death\" or \"Lost\" where the participant died or was lost to follow-up prior to the end of death record linkage in the participant's latest inferred nation of residence"),
  data.table(var="max_death_reason", name="Record source used to infer the participant's latest nation of residence or otherwise set the latest date of death record linkage"),
  data.table(var="max_hospital_censor", name="Latest date for which the participant has hospital record linkage based on date of latest complete hospital record linkage available in each nation, date of death, or loss of follow-up"),
  data.table(var="max_hospital_location", name="Inferred nation of residence used to determine latest hospital record linkage, or \"Death\" or \"Lost\" where the participant died or was lost to follow-up prior to the end of hospital record linkage in the participant's latest inferred nation of residence"),
  data.table(var="max_hospital_reason", name="Record source used to infer the participant's latest nation of residence or otherwise set the latest date of hospital record linkage"),
  data.table(var="max_cancer_censor", name="Latest date for which the participant could possibly have cancer registry data based on date of latest complete cancer registry linkage available in each nation, date of death, or loss of follow-up"),
  data.table(var="max_cancer_location", name="Inferred nation of residence used to determine latest possible cancer registry linkage, or \"Death\" or \"Lost\" where the participant died or was lost to follow-up prior to the end of cancer registry linkage in the participants latest inferred nation of residence"),
  data.table(var="max_cancer_reason", name="Record source used to infer the participant's latest nation of residence or otherwise set the latest date of cancer registry linkage"),
  data.table(var="max_gp_censor", name="Latest date for which primary care record linkage is available for this participant, if any"),
  data.table(var="max_gp_location", name="Nation of residence of the participant at the latest available primary care record, if any"),
  data.table(var="max_gp_reason", name="Type of last available primary care record available for this participant, if any")
)

##############
# write out
##############
fwrite(out, "followup/follow_up.csv")
fwrite(info, "followup/column_information.csv")

# Upload to persistent storage
system("dx upload followup/* --destination 'common/Follow-up/'")

# Send raw data to deletion folder to reduce storage costs - this needs to be
# done in two steps, since we can't move two folders with the same name to the
# same location, so here we rename with a random number and then move to trash
rn <- as.integer(Sys.time())
system(sprintf("dx mv 'common/Follow-up/raw_data/' 'common/Follow-up/raw_data_%s'", rn), wait=TRUE)
system(sprintf("dx mv 'common/Follow-up/raw_data_%s/' trash/", rn), wait=TRUE)
