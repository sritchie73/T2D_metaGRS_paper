library(data.table)
library(RNOmni)
library(survival)
source("src/functions/cox_test.R")
source("src/functions/factor.R") # set reference group (first level) without specifying all levels
source("src/functions/calendar_time.R") # for computing follow-up time in years (handling leap years)

# Make output directory
system("mkdir -p output/INTERVAL_tests", wait=TRUE)

# Load p1074 ID mapping file
idmap_p1074 <- fread("data/INTERVAL/P1074/INTERVAL_OmicsMap_20240202.csv")
idmap_p1074 <- idmap_p1074[,.(identifier, IID=Affymetrix_gwasQC_bl)]

# Drop people without genetic data (or QC'd out)
idmap_p1074 <- idmap_p1074[!is.na(IID)]

# Load p1074 phenotype data
pheno <- fread("data/INTERVAL/P1074/INTERVALdata_02FEB2024.csv")
pheno <- pheno[,.(identifier, attendance_date=as.IDate(attendanceDate, format="%d%b%Y"), age=agePulse, sex=sexPulse, weight=wt_bl, height=ht_bl)]
pheno <- pheno[idmap_p1074, on = .(identifier), nomatch = 0]

# Recode sex
pheno[, sex := ifelse(sex == 2, "Female", "Male")]

# Compute BMI and remove outliers
pheno[, bmi := weight/height^2]
pheno[weight == 777, bmi := NA_real_] # bad coding
pheno[height < 1.47, bmi := NA_real_] # clinical cutoff for dwarfism
pheno[height > 2.1, bmi := NA_real_] # clinical cutoff for gigantism
pheno[weight < 50 | weight > 160, bmi := NA_real_] # NHS restrictions for weight

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
onset[, years_to_onset := mean(c(years_between(attendance_date, last_diabetes_free_event), years_between(attendance_date, first_diabetes_event))), by=.(identifier)]
pheno[onset, on = .(identifier), incident_censor_years := i.years_to_onset]

# For people without diabetes, set the maximum follow-up time as the max follow
pheno[hes, on = .(identifier), max_follow := i.cenDate]
pheno[!(prevalent_diabetes) & !(incident_diabetes), incident_censor_years := years_between(attendance_date, max_follow)]

# Load metaGRS and existing T2D PGS
pgs <- fread("data/INTERVAL/T2D_PGS/collated_scores.sscore.gz")
metaGRS <- fread("data/INTERVAL/T2D_metaGRS/T2D_metaGRS.sscore.gz")
pgs[metaGRS, on = .(IID), T2D_metaGRS := i.score_sum]
pgs[, Gad_metaGRS_v1 := NULL]
pgs[, Gad_metaGRS_v2 := NULL]

# Add to phenotype data
pheno <- pheno[pgs, on = .(IID), nomatch=0]

# Get PCs
pcs <- fread("data/INTERVAL/genetics/reference/annot_INT_50PCs_pcs.txt")
setnames(pcs, "ID", "IID")
setnames(pcs, gsub("_", "", names(pcs)))
pcs <- pcs[,.(IID, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
              PC11, PC12, PC13, PC14, PC15, PC16, PC17, PC18, PC19, PC20)]

# Add to phenotype data
pheno <- pheno[pcs, on = .(IID), nomatch=0]

# Filter to incident event data with non-missing bmi
pheno <- pheno[!is.na(incident_diabetes) & !(prevalent_diabetes) & !is.na(bmi)]

# Ye2021_PGS001357 has effect allele and other allele incorrectly reported to PGS Catalog;
# as-is the score has strong inverse correlation with T2D status, so we need to flip the 
# sign.
pheno[, Ye2021_PGS001357 := -Ye2021_PGS001357]

# Apply BMI stratification as described in https://www.nature.com/articles/s41588-024-01782-y

# 1. Get sex-specific BMI cut-points in T2D cases
med_bmi_cases <- pheno[!is.na(bmi) & (incident_diabetes), .(median_BMI=median(bmi)), by=sex]

# 2. Determine sample sizes in each group
strata_groups <- rbind(idcol="bmi_strata",
  "low"=pheno[!is.na(incident_diabetes)][med_bmi_cases, on = .(sex, bmi < median_BMI), .N, by=.(sex, incident_diabetes)],
  "high"=pheno[!is.na(incident_diabetes)][med_bmi_cases, on = .(sex, bmi > median_BMI), .N, by=.(sex, incident_diabetes)]
)

# 3. Determine size to down-sample to
target_n <- strata_groups[,.SD[which.min(N)],by=.(sex, incident_diabetes)]

# 4. Create down-sampled dataset
inci_dat <- pheno[!is.na(bmi) & !is.na(incident_diabetes)]
inci_dat[med_bmi_cases, on = .(sex), bmi_strata := ifelse(bmi < median_BMI, "low", "high")]
inci_dat <- foreach(gIdx = strata_groups[,.I], .combine=rbind) %do% {
  this_strata <- strata_groups[gIdx]
  this_strata[, N := NULL]
  this_inci_dat <- inci_dat[this_strata, on = .(bmi_strata, sex, incident_diabetes)]
  this_target_n <- target_n[this_strata, on = .(sex, incident_diabetes), N]
  this_inci_dat[sample(.N, this_target_n)]
}

# For the HuertaChagoya2023, we need to combine their three scores as they describe at
# https://www.pgscatalog.org/score/PGS003443/
# https://www.pgscatalog.org/score/PGS003444/
# https://www.pgscatalog.org/score/PGS003445/
inci_dat[, HuertaChagoya2023_PGS00344X := scale(
  scale(HuertaChagoya2023_EUR_PGS003443)*0.531117 +
  scale(HuertaChagoya2023_EAS_PGS003444)*0.5690198 +
  scale(HuertaChagoya2023_LAT_PGS003445)*0.1465538
)]
inci_dat[, HuertaChagoya2023_EUR_PGS003443 := NULL]
inci_dat[, HuertaChagoya2023_EAS_PGS003444 := NULL]
inci_dat[, HuertaChagoya2023_LAT_PGS003445 := NULL]

# Get list of PGS for testing
test_pgs <- names(pgs)[-1]
test_pgs <- setdiff(test_pgs, c("HuertaChagoya2023_EUR_PGS003443", "HuertaChagoya2023_EAS_PGS003444", "HuertaChagoya2023_LAT_PGS003445"))
test_pgs <- c(test_pgs, "HuertaChagoya2023_PGS00344X")

# Adjust pgs for 20 PCs
for (this_pgs in test_pgs) {
  setnames(inci_dat, this_pgs, "this_pgs")
  inci_dat[, this_pgs := scale(lm(scale(this_pgs) ~ PC1 + PC2 + PC3 + PC4 +
    PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 +
    PC16 + PC17 + PC18 + PC19 + PC20)$residuals)]
  setnames(inci_dat, "this_pgs", this_pgs)
}

# Process risk factors and covariates prior to modelling
inci_dat[, age := scale(age)]
inci_dat[, sex := factor(sex, reference="Female")]

# Assess PRS prediction in each BMI strata for each PRS
mf <- "Surv(incident_censor_years, incident_diabetes) ~ %s + age + sex"
prs_inci <- foreach(this_strata = c("low", "high"), .combine=rbind) %do% {
  foreach(this_pgs = test_pgs, .combine=rbind) %do% {
		res <- cox.test(sprintf(mf, this_pgs), "incident_diabetes", inci_dat[bmi_strata == this_strata])
		res[, model := this_pgs]
    res[, bmi_strata := this_strata]
    return(res[,])
  }
}

# Rename coefficients for readability
prs_inci[coefficient == "age", coefficient := "Age"]
prs_inci[coefficient == "sexMale", coefficient := "Sex: Male vs. Female"]
