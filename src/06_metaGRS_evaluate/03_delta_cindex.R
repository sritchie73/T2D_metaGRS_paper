library(data.table)
library(foreach)
library(survival)
library(boot)
library(doMC)
registerDoMC(10) # Takes about 20 minutes on icelake, with 20 cores requested (mem > cpu)
source("src/functions/cox_test.R")
source("src/functions/factor_by_size.R") # set reference group (first level) as largest group
source("src/functions/logit.R") # like log transform, but for percentages

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

# Filter to columns we actually need
pheno <- pheno[,.(incident_type_2_diabetes, incident_censor_years, age, sex,
  QDiabetes2013, QDiabetes2018A, QDiabetes2018B_non_fasting, QDiabetes2018C, T2D_metaGRS,
  PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, PC11, PC12, PC13, PC14, PC15, PC16, PC17, PC18, PC19, PC20)]

# Bypass NA check for censboot
pheno[is.na(QDiabetes2018B_non_fasting), QDiabetes2018B_non_fasting := -1]

boot.fun <- function(dt) {
  # We need two sets of data to handle the fact that QDiabetes2018 model B has missing data
  dtB <- dt[QDiabetes2018B_non_fasting > 0]

  # Adjust pgs for 20 PCs
  dt[, T2D_metaGRS := lm(scale(T2D_metaGRS) ~ PC1 + PC2 + PC3 + PC4 +
    PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 +
    PC16 + PC17 + PC18 + PC19 + PC20)$residuals]

  dtB[, T2D_metaGRS := lm(scale(T2D_metaGRS) ~ PC1 + PC2 + PC3 + PC4 +
    PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 +
    PC16 + PC17 + PC18 + PC19 + PC20)$residuals]

  # Get C-index for QDiabetes models 
  qdiab2013 <- summary(coxph(Surv(incident_censor_years, incident_type_2_diabetes) ~ scale(logit(QDiabetes2013/100)), data=dt))$concordance[1]
  qdiabA <- summary(coxph(Surv(incident_censor_years, incident_type_2_diabetes) ~ scale(logit(QDiabetes2018A/100)), data=dt))$concordance[1]
  qdiabB <- summary(coxph(Surv(incident_censor_years, incident_type_2_diabetes) ~ scale(logit(QDiabetes2018B_non_fasting/100)), data=dtB))$concordance[1]
  qdiabC <- summary(coxph(Surv(incident_censor_years, incident_type_2_diabetes) ~ scale(logit(QDiabetes2018C/100)), data=dt))$concordance[1]

  # Get C-index when adding metaPRS
  prs2013 <- summary(coxph(Surv(incident_censor_years, incident_type_2_diabetes) ~ scale(logit(QDiabetes2013/100)) + scale(T2D_metaGRS), data=dt))$concordance[1]
  prsA <- summary(coxph(Surv(incident_censor_years, incident_type_2_diabetes) ~ scale(logit(QDiabetes2018A/100)) + scale(T2D_metaGRS), data=dt))$concordance[1]
  prsB <- summary(coxph(Surv(incident_censor_years, incident_type_2_diabetes) ~ scale(logit(QDiabetes2018B_non_fasting/100)) + scale(T2D_metaGRS), data=dtB))$concordance[1]
  prsC <- summary(coxph(Surv(incident_censor_years, incident_type_2_diabetes) ~ scale(logit(QDiabetes2018C/100)) + scale(T2D_metaGRS), data=dt))$concordance[1]

  # Compute the four delta C-indices 
  c(prs2013 - qdiab2013, prsA - qdiabA, prsB - qdiabB, prsC - qdiabC)
}

# Run bootstrap analysis 
surv_cols_idx <- match(c("incident_censor_years", "incident_type_2_diabetes"), names(pheno))
boot_res <- censboot(pheno, boot.fun, 1000, index=surv_cols_idx)
saveRDS(boot_res, file="output/UKB_tests/delta_cindices_bootstraps.rds")

# Curate results
delta_cind <- data.table(
  prs_added_to = c("QDiabetes2013", "QDiabetes2018A", "QDiabetes2018B_non_fasting", "QDiabetes2018C"), 
  delta_C = boot_res$t0, boot_SE = apply(boot_res$t, 2, sd)
)
delta_cind[, boot_L95 := delta_C - qnorm(1-(0.05/2))*boot_SE]
delta_cind[, boot_U95 := delta_C + qnorm(1-(0.05/2))*boot_SE]
delta_cind[, boot_pval := pmin(1, pnorm(abs(delta_C/boot_SE), lower.tail=FALSE)*2)]
fwrite(delta_cind, sep="\t", quote=FALSE, file="output/UKB_tests/delta_cind_bootstraps.txt")













