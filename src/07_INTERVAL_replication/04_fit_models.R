library(data.table)
library(survival)
library(foreach)
source("src/functions/glm_test.R")
source("src/functions/cox_test.R")
source("src/functions/factor.R") # set reference group (first level) without specifying all levels
source("src/functions/factor_by_size.R") # set reference group (first level) as largest group
source("src/functions/calendar_time.R") # for computing follow-up time in years (handling leap years)

# Make output directory
system("mkdir -p output/INTERVAL_tests", wait=TRUE)

# Load p1074 ID mapping file
idmap_p1074 <- fread("data/INTERVAL/p1074/omicsMap.csv")
idmap_p1074 <- idmap_p1074[,.(p1074_id=identifier, IID=Affymetrix_gwasQC_bl)]

# Drop people without genetic data (or QC'd out)
idmap_p1074 <- idmap_p1074[!is.na(IID)]

# Load p1074 phenotype data
pheno <- fread("data/INTERVAL/p1074/INTERVALdata_17JUL2018.csv")
pheno <- pheno[,.(p1074_id=identifier, attendance_date=as.IDate(attendanceDate, format="%d%b%Y"), age=agePulse, sex=sexPulse)]
pheno <- pheno[idmap_p1074, on = .(p1074_id), nomatch = 0]

# Load HES data ID mapping file
fname <- list.files(path="data/INTERVAL/HES", pattern="INTERVAL_OmicsMap_[0-9]*.csv", full.names=TRUE)
idmap_hes <- fread(fname)
idmap_hes <- idmap_hes[,.(hes_id=identifier, IID=Affymetrix_gwasQC_bl)]

# Merge in identifiers
pheno <- pheno[idmap_hes, on = .(IID), nomatch=0]

# Load in HES data
fname <- list.files(path="data/INTERVAL/HES", pattern="INTERVAL_CaliberEP_.*.csv", full.names=TRUE)
hes <- fread(fname, na.strings=c("", "NA"))

# Need to do some processing to get follow-up time
hes <- melt(hes, id.vars="identifier", value.name="event_date", na.rm=TRUE)
hes[, event_date := as.IDate(event_date, format="%d%b%Y")]
maxfollow <- hes[,max(event_date)]

# Extract diabetes (see caliber_event_names.csv file)
diabetes <- hes[variable == "cal_ps_d_062", .(hes_id=identifier, event_date)]

# Add to pheno
pheno[diabetes, on = .(hes_id), diabetes_date := i.event_date]

# Determine prevalent and incident diabetes status
pheno[, prevalent_diabetes := fcase(
  diabetes_date < attendance_date, TRUE,
  diabetes_date > attendance_date, FALSE,
  is.na(diabetes_date), FALSE
)]

pheno[, incident_diabetes := fcase(
  diabetes_date < attendance_date, NA,
  diabetes_date > attendance_date, TRUE,
  is.na(diabetes_date), FALSE
)]

# Determine follow-up time
pheno[, incident_censor_date := diabetes_date]
pheno[is.na(incident_censor_date), incident_censor_date := maxfollow]
pheno[(prevalent_diabetes), incident_censor_date := NA]
pheno[, incident_censor_years := years_between(attendance_date, incident_censor_date)]

# Load metaGRS and existing T2D PGS
pgs <- fread("data/INTERVAL/T2D_PGS/collated_scores.sscore.gz")
metaGRS <- fread("data/INTERVAL/T2D_metaGRS/T2D_metaGRS.sscore.gz")
pgs[metaGRS, on = .(IID), T2D_metaGRS := i.score_sum]

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

# Adjust pgs for 20 PCs
for (this_pgs in names(pgs)[-1]) {
  setnames(pheno, this_pgs, "this_pgs")
  pheno[, this_pgs := scale(lm(scale(this_pgs) ~ PC1 + PC2 + PC3 + PC4 + 
    PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + 
    PC16 + PC17 + PC18 + PC19 + PC20)$residuals)]
  setnames(pheno, "this_pgs", this_pgs)
}

# Ye2021_PGS001357 has effect allele and other allele incorrectly reported to PGS Catalog;
# as-is the score has strong inverse correlation with T2D status, so we need to flip the 
# sign.
pheno[, Ye2021_PGS001357 := -Ye2021_PGS001357]

# Process risk factors and covariates prior to modelling - we want all
# estimates to be per SD change in variable, and for factors, using the
# largest group or lowest risk as reference 
pheno[, age := scale(age)]
pheno[, sex := ifelse(sex == 2, "Female", "Male")]
pheno[, sex := factor(sex, reference="Female")]

# First examine prevalent T2D
prev <- rbind(idcol="model",
  "age"=glm.test(prevalent_diabetes ~ age, "prevalent_diabetes", pheno),
  "sex"=glm.test(prevalent_diabetes ~ sex, "prevalent_diabetes", pheno),
  "age + sex"=glm.test(prevalent_diabetes ~ age + sex, "prevalent_diabetes", pheno)
)

mf <- "prevalent_diabetes ~ %s + age + sex"
prs_prev <- foreach(this_pgs = names(pgs)[-1], .combine=rbind) %do% {
  res <- glm.test(sprintf(mf, this_pgs), "prevalent_diabetes", pheno)
  res[, model := this_pgs]
} 

prev <- rbind(idcol="model_type", "reference"=prev, "pgs"=prs_prev)

# Rename coefficients for readability
prev[coefficient == "age", coefficient := "Age"]
prev[coefficient == "sexMale", coefficient := "Sex: Male vs. Female"]

# Write out
fwrite(prev, sep="\t", quote=FALSE, file="output/INTERVAL_tests/prevalent_diabetes_associations.txt")

# Now do incident T2D
inci <- rbind(idcol="model",
  "age"=cox.test("Surv(incident_censor_years, incident_diabetes) ~ age", "incident_diabetes", pheno),
  "sex"=cox.test("Surv(incident_censor_years, incident_diabetes) ~ sex", "incident_diabetes", pheno),
  "age + sex"=cox.test("Surv(incident_censor_years, incident_diabetes) ~ age + sex", "incident_diabetes", pheno),
  "strata(sex) + age"=cox.test("Surv(incident_censor_years, incident_diabetes) ~ age + strata(sex)", "incident_diabetes", pheno)
)

mf <- "Surv(incident_censor_years, incident_diabetes) ~ %s + age + sex"
pgs_inci <- foreach(this_pgs = names(pgs)[-1], .combine=rbind) %do% {
  res <- cox.test(sprintf(mf, this_pgs), "incident_diabetes", pheno)
  res[, model := this_pgs]
}

inci <- rbind(idcol="model_type", "reference"=inci, "PGS"=pgs_inci)

# Rename coefficients for readability
inci[coefficient == "age", coefficient := "Age"]
inci[coefficient == "sexMale", coefficient := "Sex: Male vs. Female"]

# Write out
fwrite(inci, sep="\t", quote=FALSE, file="output/UKB_tests/incident_diabetes_associations.txt")
