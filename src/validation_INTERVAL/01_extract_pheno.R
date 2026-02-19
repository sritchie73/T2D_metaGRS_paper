library(data.table)
source("src/functions/calendar_time.R") # for computing follow-up time in years (handling leap years)

# Load p1074 ID mapping file
idmap_p1074 <- fread("data/INTERVAL/P1074/INTERVAL_OmicsMap_20240202.csv")
idmap_p1074 <- idmap_p1074[,.(identifier, IID=Affymetrix_gwasQC_bl)]

# Drop people without genetic data (or QC'd out)
idmap_p1074 <- idmap_p1074[!is.na(IID)]

# Load p1074 phenotype data
pheno <- fread("data/INTERVAL/P1074/INTERVALdata_02FEB2024.csv")
pheno <- pheno[,.(identifier, attendance_date=as.IDate(attendanceDate, format="%d%b%Y"), age=agePulse, sex=sexPulse)]
pheno <- pheno[idmap_p1074, on = .(identifier), nomatch = 0]

# Load in HES data
hes <- fread("data/INTERVAL/P1074/Caliber_02FEB2024.csv", na.strings=c("", "NA"))
hes[, cenDate := as.IDate(cenDate, format="%d%b%Y")]

# Melt to long to get a sequence of events for each person
events <- copy(hes)
events[, cenDate := NULL]
events <- melt(events, id.vars="identifier", value.name="event_date", na.rm=TRUE)
events[, event_date := as.IDate(event_date, format="%d%b%Y")]

# Add in assessment date
events <- rbind(events, pheno[,.(identifier, variable="baseline", event_date=attendance_date)])

# Order events
events <- events[order(event_date)][order(identifier)]

# Flag diabetes events (see caliber_event_names.csv file)
events[, diabetes_event := ifelse(variable %in% c("cal_ps_62", "cal_p_62"), TRUE, FALSE)]

# Create indicator that creates a flag for every event indicating whether the person has had
# a diabetes event already
events[, diabetic := as.logical(cumsum(diabetes_event)), by=.(identifier)]

# Flag those with diabetes prior to baseline:
pheno[events[variable == "baseline"], on = .(identifier), prevalent_diabetes := i.diabetic] 

# Flag people with incident diabetes
pheno[, incident_diabetes := FALSE]
pheno[events[(diabetic)], on = .(identifier), incident_diabetes := TRUE]
pheno[(prevalent_diabetes), incident_diabetes := NA]

# Set follow-up as the mid-point between the first diabetes event and the previous diabetes-free event
inci_diab <- events[identifier %in% pheno[(incident_diabetes), identifier]]
first_event <- inci_diab[(diabetic),.SD[1],by=identifier]
last_diab_free <- inci_diab[!(diabetic), .SD[.N], by=identifier]
onset <- pheno[(incident_diabetes),.(identifier, attendance_date)]
onset[last_diab_free, on=.(identifier), last_diabetes_free_event := i.event_date]
onset[first_event, on = .(identifier), first_diabetes_event := i.event_date]
onset[, years_to_onset := years_between(attendance_date, midpoint(last_diabetes_free_event, first_diabetes_event)), by=.(identifier)]
pheno[onset, on = .(identifier), incident_censor_years := i.years_to_onset]

# For people without diabetes, set the maximum follow-up time as the max follow
pheno[hes, on = .(identifier), max_follow := i.cenDate]
pheno[!(prevalent_diabetes) & !(incident_diabetes), incident_censor_years := years_between(attendance_date, max_follow)]

# Write out phenotype data
fwrite(pheno, sep="\t", quote=FALSE, file="data/INTERVAL/curated_pheno.tsv")
