library(data.table)
library(RNOmni)
library(survival)
library(boot)
options(boot.parallel="multicore")
options(boot.ncpus=10) # Takes about 20 minutes on icelake, with 20 cores requested (mem > cpu)
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

test_pgs <- names(pgs)[-which(names(pgs) %in% c("IID", "T2D_metaGRS"))]

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

# Filter to incident event data
pheno <- pheno[!is.na(incident_diabetes) & !(prevalent_diabetes)]

# For the HuertaChagoya2023, we need to combine their three scores as they describe at
# https://www.pgscatalog.org/score/PGS003443/
# https://www.pgscatalog.org/score/PGS003444/
# https://www.pgscatalog.org/score/PGS003445/
pheno[, HuertaChagoya2023_PGS00344X := scale(
  scale(HuertaChagoya2023_EUR_PGS003443)*0.531117 +
  scale(HuertaChagoya2023_EAS_PGS003444)*0.5690198 +
  scale(HuertaChagoya2023_LAT_PGS003445)*0.1465538
)]
pheno[, HuertaChagoya2023_EUR_PGS003443 := NULL]
pheno[, HuertaChagoya2023_EAS_PGS003444 := NULL]
pheno[, HuertaChagoya2023_LAT_PGS003445 := NULL]
test_pgs <- c(setdiff(test_pgs, c("HuertaChagoya2023_EUR_PGS003443", "HuertaChagoya2023_EAS_PGS003444", "HuertaChagoya2023_LAT_PGS003445")), "HuertaChagoya2023_PGS00344X")

# Adjust pgs for 20 PCs
pgs_list <- c(intersect(names(pgs)[-1], names(pheno)), "HuertaChagoya2023_PGS00344X")
for (this_pgs in pgs_list) {
  setnames(pheno, this_pgs, "this_pgs")
  pheno[, this_pgs := scale(lm(scale(this_pgs) ~ PC1 + PC2 + PC3 + PC4 +
    PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 +
    PC16 + PC17 + PC18 + PC19 + PC20)$residuals)]
  setnames(pheno, "this_pgs", this_pgs)
}

# Select columns of interest
pheno <- pheno[, .SD, .SDcols=c("IID", "age", "sex", "incident_censor_years", "incident_diabetes", "T2D_metaGRS", test_pgs)]

# bootstrap delta C-index relative to metaPRS
boot.fun <- function(dt) {
  # function to fit Cox proportional hazards model and get C-index
  compute_Cind <- function(PRS, dat) {
    cx <- coxph(as.formula(sprintf("Surv(incident_censor_years, incident_diabetes) ~ scale(%s) + age + sex", PRS)), data=dat)
    summary(cx)$concordance[1]
  }

  # Compute C-index for metaPRS
  cind.ref <- compute_Cind("T2D_metaGRS", dt)
  
  # Compute C-index for test PRSs
  cind.prs <- sapply(test_pgs, compute_Cind, dt)

  # Return delta C-indices
  as.vector(cind.prs - cind.ref)
}

# Run bootstrap analysis 
surv_cols_idx <- match(c("incident_censor_years", "incident_diabetes"), names(pheno))
boot_res <- censboot(pheno, boot.fun, 1000, index=surv_cols_idx)
saveRDS(boot_res, file="output/INTERVAL_tests/metaPRS_delta_cindices_bootstraps.rds")

# Curate results
delta_cind <- data.table(
  prs = test_pgs, delta_C = boot_res$t0, boot_SE = apply(boot_res$t, 2, sd)
)
delta_cind[, boot_L95 := delta_C - qnorm(1-(0.05/2))*boot_SE]
delta_cind[, boot_U95 := delta_C + qnorm(1-(0.05/2))*boot_SE]
delta_cind[, boot_pval := pmin(1, pnorm(abs(delta_C/boot_SE), lower.tail=FALSE)*2)]
fwrite(delta_cind, sep="\t", quote=FALSE, file="output/INTERVAL_tests/metaPRS_delta_cind_bootstraps.txt")

