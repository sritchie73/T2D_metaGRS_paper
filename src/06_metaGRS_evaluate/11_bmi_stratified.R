library(data.table)
library(survival)
library(foreach)
library(RNOmni)
source("src/functions/glm_test.R")
source("src/functions/cox_test.R")
source("src/functions/factor.R") # set reference group (first level) without specifying all levels
source("src/functions/factor_by_size.R") # set reference group (first level) as largest group

# Make output directory
system("mkdir -p output/UKB_tests", wait=TRUE)

# Load phenotype data
pheno <- fread("data/UKB/collated_curated_data.txt", tmpdir="tmp", na.strings=c("", "NA"))
pheno <- pheno[(metaGRS_test_samples)]

# Load metaGRS and existing T2D PGS
pgs <- fread("data/UKB/T2D_PGS/collated_scores.sscore.gz")
metaGRS <- fread("data/UKB/T2D_metaGRS_levels.txt")
pgs <- metaGRS[pgs, on = .(eid=IID)]

# Drop invalid PRS
pgs[, c("Mars2022_AJHG_PGS002771", "Gad_metaGRS_v1", "Gad_metaGRS_v2") := NULL]

# Add to phenotype data
pheno <- pheno[pgs, on = .(eid), nomatch=0]

# Ye2021_PGS001357 has effect allele and other allele incorrectly reported to PGS Catalog;
# as-is the score has strong inverse correlation with T2D status, so we need to flip the 
# sign.
pheno[, Ye2021_PGS001357 := -Ye2021_PGS001357]

# Apply BMI stratification as described in https://www.nature.com/articles/s41588-024-01782-y

# 1. Get sex-specific BMI cut-points in T2D cases
med_bmi_cases <- pheno[!is.na(bmi) & (type_2_diabetes), .(median_BMI=median(bmi)), by=sex]

# 2. Determine sample sizes in each group
strata_groups <- rbind(idcol="bmi_strata",
  "low"=pheno[!is.na(type_2_diabetes)][med_bmi_cases, on = .(sex, bmi < median_BMI), .N, by=.(sex, type_2_diabetes)],
  "high"=pheno[!is.na(type_2_diabetes)][med_bmi_cases, on = .(sex, bmi > median_BMI), .N, by=.(sex, type_2_diabetes)]
)

# 3. Determine size to down-sample to
target_n <- strata_groups[,.SD[which.min(N)],by=.(sex, type_2_diabetes)]

# 4. Create down-sampled dataset
prev_dat <- pheno[!is.na(bmi) & !is.na(type_2_diabetes)]
prev_dat[med_bmi_cases, on = .(sex), bmi_strata := ifelse(bmi < median_BMI, "low", "high")]
prev_dat <- foreach(gIdx = strata_groups[,.I], .combine=rbind) %do% {
  this_strata <- strata_groups[gIdx]
  this_strata[, N := NULL]
  this_prev_dat <- prev_dat[this_strata, on = .(bmi_strata, sex, type_2_diabetes)]
  this_target_n <- target_n[this_strata, on = .(sex, type_2_diabetes), N]
  this_prev_dat[sample(.N, this_target_n)]
}

# Adjust pgs for 20 PCs
for (this_pgs in names(pgs)[-1]) {
  setnames(prev_dat, this_pgs, "this_pgs")
  prev_dat[, this_pgs := scale(lm(scale(this_pgs) ~ PC1 + PC2 + PC3 + PC4 +
    PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 +
    PC16 + PC17 + PC18 + PC19 + PC20)$residuals)]
  setnames(prev_dat, "this_pgs", this_pgs)
}

# Process risk factors and covariates prior to modelling
prev_dat[, age := scale(age)]
prev_dat[, sex := factor(sex, reference="Female")]

# Assess PRS prediction in each BMI strata for each PRS
mf <- "type_2_diabetes ~ %s + age + sex"
prs_prev <- foreach(this_strata = c("low", "high"), .combine=rbind) %do% {
  foreach(this_pgs = names(pgs)[-1], .combine=rbind) %do% {
		res <- glm.test(sprintf(mf, this_pgs), "type_2_diabetes", prev_dat[bmi_strata == this_strata])
		res[, model := this_pgs]
    res[, bmi_strata := this_strata]
    return(res[,])
  }
}

# Rename coefficients for readability
prs_prev[coefficient == "age", coefficient := "Age"]
prs_prev[coefficient == "sexMale", coefficient := "Sex: Male vs. Female"]





# Now repeat the process for incident T2D

# 1. Get sex-specific BMI cut-points in T2D cases
med_bmi_cases <- pheno[!is.na(bmi) & (incident_type_2_diabetes), .(median_BMI=median(bmi)), by=sex]

# 2. Determine sample sizes in each group
strata_groups <- rbind(idcol="bmi_strata",
  "low"=pheno[!is.na(incident_type_2_diabetes)][med_bmi_cases, on = .(sex, bmi < median_BMI), .N, by=.(sex, incident_type_2_diabetes)],
  "high"=pheno[!is.na(incident_type_2_diabetes)][med_bmi_cases, on = .(sex, bmi > median_BMI), .N, by=.(sex, incident_type_2_diabetes)]
)

# 3. Determine size to down-sample to
target_n <- strata_groups[,.SD[which.min(N)],by=.(sex, incident_type_2_diabetes)]

# 4. Create down-sampled dataset
inci_dat <- pheno[!is.na(bmi) & !is.na(incident_type_2_diabetes)]
inci_dat[med_bmi_cases, on = .(sex), bmi_strata := ifelse(bmi < median_BMI, "low", "high")]
inci_dat <- foreach(gIdx = strata_groups[,.I], .combine=rbind) %do% {
  this_strata <- strata_groups[gIdx]
  this_strata[, N := NULL]
  this_inci_dat <- inci_dat[this_strata, on = .(bmi_strata, sex, incident_type_2_diabetes)]
  this_target_n <- target_n[this_strata, on = .(sex, incident_type_2_diabetes), N]
  this_inci_dat[sample(.N, this_target_n)]
}

# Adjust pgs for 20 PCs
for (this_pgs in names(pgs)[-1]) {
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
mf <- "Surv(incident_censor_years, incident_type_2_diabetes) ~ %s + age + sex"
prs_inci <- foreach(this_strata = c("low", "high"), .combine=rbind) %do% {
  foreach(this_pgs = names(pgs)[-1], .combine=rbind) %do% {
		res <- cox.test(sprintf(mf, this_pgs), "incident_type_2_diabetes", inci_dat[bmi_strata == this_strata])
		res[, model := this_pgs]
    res[, bmi_strata := this_strata]
    return(res[,])
  }
}

# Rename coefficients for readability
prs_inci[coefficient == "age", coefficient := "Age"]
prs_inci[coefficient == "sexMale", coefficient := "Sex: Male vs. Female"]







