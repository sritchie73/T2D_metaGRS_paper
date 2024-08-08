library(data.table)
library(foreach)
library(survival)
library(boot)
options(boot.parallel="multicore")
options(boot.ncpus=10) # Takes about 20 minutes on icelake, with 20 cores requested (mem > cpu)
source("src/functions/cox_test.R")
source("src/functions/factor_by_size.R") # set reference group (first level) as largest group

# Make output directory
system("mkdir -p output/UKB_tests", wait=TRUE)

# Load phenotype data 
pheno <- fread("data/UKB/collated_curated_data.txt", tmpdir="tmp", na.strings=c("", "NA"))
pheno <- pheno[(metaGRS_test_samples)]

# Load metaGRS 
metaGRS <- fread("data/UKB/T2D_metaGRS_levels.txt")
pheno <- pheno[metaGRS, on = .(eid), nomatch=0]

# Drop participants without QDiabetes data 
pheno <- pheno[!is.na(QDiabetes2018C)]

# Drop participants with uncertain/unknowable diabetes status
pheno <- pheno[!is.na(incident_type_2_diabetes)]

# add in age tertiles for later QDiabetes stratification analysis
pheno[!is.na(QDiabetes2018C), age_tertile := cut(age, breaks=quantile(age, probs=seq(0, 1, length=4)), include.lowest=TRUE)]

# Replace QDiabetes absrisk with linear predictors
lps <- fread("output/UKB_tests/QDiabetes_plus_metaPRS_linear_predictors.txt")
pheno[lps, on = .(eid), QDiabetes2018A := QDiabA_LP]
pheno[lps, on = .(eid), QDiabetes2018B_non_fasting := QDiabB_LP]
pheno[lps, on = .(eid), QDiabetes2018C := QDiabC_LP]

# Filter to columns we actually need
pheno <- pheno[,.(incident_type_2_diabetes, incident_censor_years, age, age_tertile, sex,
  QDiabetes2018A, QDiabetes2018B_non_fasting, QDiabetes2018C, T2D_metaGRS,
  PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, PC11, PC12, PC13, PC14, PC15, PC16, PC17, PC18, PC19, PC20)]

# Bypass NA check for censboot
pheno[is.na(QDiabetes2018B_non_fasting), QDiabetes2018B_non_fasting := -999]

boot.fun <- function(dt) {
  # We need two sets of data to handle the fact that QDiabetes2018 model B has missing data
  dtB <- dt[QDiabetes2018B_non_fasting > -999]

  # Adjust pgs for 20 PCs
  dt[, T2D_metaGRS := lm(scale(T2D_metaGRS) ~ PC1 + PC2 + PC3 + PC4 +
    PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 +
    PC16 + PC17 + PC18 + PC19 + PC20)$residuals]

  dtB[, T2D_metaGRS := lm(scale(T2D_metaGRS) ~ PC1 + PC2 + PC3 + PC4 +
    PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 +
    PC16 + PC17 + PC18 + PC19 + PC20)$residuals]

  # Get C-index for QDiabetes models 
  qdiabA <- summary(coxph(Surv(incident_censor_years, incident_type_2_diabetes) ~ scale(QDiabetes2018A), data=dt))$concordance[1]
  qdiabB <- summary(coxph(Surv(incident_censor_years, incident_type_2_diabetes) ~ scale(QDiabetes2018B_non_fasting), data=dtB))$concordance[1]
  qdiabC <- summary(coxph(Surv(incident_censor_years, incident_type_2_diabetes) ~ scale(QDiabetes2018C), data=dt))$concordance[1]

  # Get C-index when adding metaPRS
  prsA <- summary(coxph(Surv(incident_censor_years, incident_type_2_diabetes) ~ scale(QDiabetes2018A) + scale(T2D_metaGRS), data=dt))$concordance[1]
  prsB <- summary(coxph(Surv(incident_censor_years, incident_type_2_diabetes) ~ scale(QDiabetes2018B_non_fasting) + scale(T2D_metaGRS), data=dtB))$concordance[1]
  prsC <- summary(coxph(Surv(incident_censor_years, incident_type_2_diabetes) ~ scale(QDiabetes2018C) + scale(T2D_metaGRS), data=dt))$concordance[1]

  # Compute the four delta C-indices 
  c(prsA - qdiabA, prsB - qdiabB, prsC - qdiabC)
}

# Run bootstrap analysis 
surv_cols_idx <- match(c("incident_censor_years", "incident_type_2_diabetes"), names(pheno))
boot_res <- censboot(pheno, boot.fun, 1000, index=surv_cols_idx)
saveRDS(boot_res, file="output/UKB_tests/delta_cindices_bootstraps.rds")

# Curate results
delta_cind <- data.table(
  prs_added_to = c("QDiabetes2018A", "QDiabetes2018B_non_fasting", "QDiabetes2018C"), 
  delta_C = boot_res$t0, boot_SE = apply(boot_res$t, 2, sd)
)
delta_cind[, boot_L95 := delta_C - qnorm(1-(0.05/2))*boot_SE]
delta_cind[, boot_U95 := delta_C + qnorm(1-(0.05/2))*boot_SE]
delta_cind[, boot_pval := pmin(1, pnorm(abs(delta_C/boot_SE), lower.tail=FALSE)*2)]
fwrite(delta_cind, sep="\t", quote=FALSE, file="output/UKB_tests/delta_cind_bootstraps.txt")

# Now do the age-tertile stratified bootstraps
boot_res <- lapply(unique(pheno$age_tertile), function(atidx) {
  censboot(pheno[age_tertile == atidx], boot.fun, 1000, index=surv_cols_idx)
})
names(boot_res) <- unique(pheno$age_tertile)
saveRDS(boot_res, file="output/UKB_tests/age_tertiles_delta_cindices_bootstraps.rds")

delta_cind <- rbindlist(lapply(unique(pheno$age_tertile), function(atidx) {
  data.table(
		prs_added_to = c("QDiabetes2018A", "QDiabetes2018B_non_fasting", "QDiabetes2018C"),
    age_tertile = atidx, delta_C = boot_res[[atidx]][["t0"]], boot_SE = apply(boot_res[[atidx]][["t"]], 2, sd)
	)
}))
delta_cind[, boot_L95 := delta_C - qnorm(1-(0.05/2))*boot_SE]
delta_cind[, boot_U95 := delta_C + qnorm(1-(0.05/2))*boot_SE]
delta_cind[, boot_pval := pmin(1, pnorm(abs(delta_C/boot_SE), lower.tail=FALSE)*2)]
fwrite(delta_cind, sep="\t", quote=FALSE, file="output/UKB_tests/age_tertiles_delta_cind_bootstraps.txt")

