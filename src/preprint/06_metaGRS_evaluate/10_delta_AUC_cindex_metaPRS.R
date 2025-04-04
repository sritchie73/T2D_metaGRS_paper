library(data.table)
library(RNOmni)
library(survival)
library(boot)
library(MASS)
library(pROC)
options(boot.parallel="multicore")
options(boot.ncpus=10) # Takes about 20 minutes on icelake, with 20 cores requested (mem > cpu)
source("src/functions/factor.R") # set reference group (first level) without specifying all levels

# Make output directory
system("mkdir -p output/UKB_tests", wait=TRUE)

# Load phenotype data
pheno <- fread("data/UKB/collated_curated_data.txt", tmpdir="tmp", na.strings=c("", "NA"))
pheno <- pheno[(metaGRS_test_samples)]

# Load metaGRS and existing T2D PGS
pgs <- fread("data/UKB/T2D_PGS/collated_scores.sscore.gz")
metaGRS <- fread("data/UKB/T2D_metaGRS_levels.txt")
pgs <- metaGRS[pgs, on = .(eid=IID)]
pgs[, Gad_metaGRS_v1 := NULL]
pgs[, Gad_metaGRS_v2 := NULL]
pgs[, Mars2022_AJHG_PGS002771 := NULL] # includes UKB GWAS 

# Get names of test PRSs
test_pgs <- names(pgs)[-c(1:2)]

# Add to phenotype data
pheno <- pheno[pgs, on = .(eid), nomatch=0]

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

# Process risk factors and covariates prior to modelling 
pheno[, age := scale(age)]
pheno[, sex := factor(sex, reference="Female")]

# Drop participants free of T2D at baseline who have withdrawn consent for hospital record
# linkage and people with uncertain diabetes status
pheno <- pheno[!is.na(type_2_diabetes) & !is.na(earliest_hospital_date)]

# Filter to relevant columns for prevalent disease analysis
prev_pheno <- pheno[,.SD,.SDcols=c("eid", "age", "sex", "type_2_diabetes", names(pgs)[-1])]

# Bootstrap delta AUCs relative to metaPRS
boot.fun <- function(dt, idx) {
  # Function to fit logistic regression and get AUC
  compute_AUC <- function(PRS, dat) {
    g1 <- glm(as.formula(sprintf("type_2_diabetes ~ scale(%s) + age + sex", PRS)), data=dat, family="binomial") 
    suppressMessages(c(auc(dat[["type_2_diabetes"]], predict(g1, dat, type="response"))))
  }

  # Compute AUC for metaPRS
  auc.ref <- compute_AUC("T2D_metaGRS", dt[idx])
  
  # Compute AUCs for test PRSs
  auc.prs <- sapply(test_pgs, compute_AUC, dt[idx])

  # Return delta AUCs
  as.vector(auc.prs - auc.ref)
}

# Run bootstrap analysis 
boot_res <- boot(prev_pheno, boot.fun, 1000)
saveRDS(boot_res, file="output/UKB_tests/metaPRS_delta_auc_bootstraps.rds")

# Curate results
delta_AUC <- data.table(
  prs = test_pgs, delta_AUC = boot_res$t0, boot_SE = apply(boot_res$t, 2, sd)
)
delta_AUC[, boot_L95 := delta_AUC - qnorm(1-(0.05/2))*boot_SE]
delta_AUC[, boot_U95 := delta_AUC + qnorm(1-(0.05/2))*boot_SE]
delta_AUC[, boot_pval := pmin(1, pnorm(abs(delta_AUC/boot_SE), lower.tail=FALSE)*2)]
fwrite(delta_AUC, sep="\t", quote=FALSE, file="output/UKB_tests/metaPRS_delta_AUC_bootstraps.txt")

# Now compute delta C-index for incident T2D
inci_pheno <- pheno[!is.na(incident_type_2_diabetes), .SD, .SDcols=c("eid", "age", "sex", "incident_censor_years", "incident_type_2_diabetes", names(pgs)[-1])]

# bootstrap delta C-index relative to metaPRS
boot.fun <- function(dt) {
  # function to fit Cox proportional hazards model and get C-index
  compute_Cind <- function(PRS, dat) {
    cx <- coxph(as.formula(sprintf("Surv(incident_censor_years, incident_type_2_diabetes) ~ scale(%s) + age + sex", PRS)), data=dat)
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
surv_cols_idx <- match(c("incident_censor_years", "incident_type_2_diabetes"), names(inci_pheno))
boot_res <- censboot(inci_pheno, boot.fun, 1000, index=surv_cols_idx)
saveRDS(boot_res, file="output/UKB_tests/metaPRS_delta_cindices_bootstraps.rds")

# Curate results
delta_cind <- data.table(
  prs = test_pgs, delta_C = boot_res$t0, boot_SE = apply(boot_res$t, 2, sd)
)
delta_cind[, boot_L95 := delta_C - qnorm(1-(0.05/2))*boot_SE]
delta_cind[, boot_U95 := delta_C + qnorm(1-(0.05/2))*boot_SE]
delta_cind[, boot_pval := pmin(1, pnorm(abs(delta_C/boot_SE), lower.tail=FALSE)*2)]
fwrite(delta_cind, sep="\t", quote=FALSE, file="output/UKB_tests/metaPRS_delta_cind_bootstraps.txt")

