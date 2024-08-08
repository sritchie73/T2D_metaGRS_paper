library(data.table)
library(RNOmni)
library(survival)
library(foreach)
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

# Filter to people with genetic data and prune first degree relatives
pheno <- pheno[!(first_degree_relative_to_drop)]

# Remove people in train or test sets
pheno <- pheno[!(ldpred2_samples) & !(metaGRS_train_samples) & !(metaGRS_test_samples)]

# Define major non-European ancestry groups similar to the UKB "White British": based on both
# self-reported ethnicity and genetic ancestry clustering
pheno[, ethnicity := fcase(
  ethnicity_subgroup %in% c("African", "Caribbean", "Black or Black British", "Any other Black background"), "AFR",
  ethnicity_subgroup %in% c("Bangladeshi", "Indian", "Pakistani"), "SAS",
  ethnicity_subgroup %in% c("Chinese"), "EAS",
  default = "NA"
)]

ancestry <- fread("data/Carles_UKB_ancestry/UKB_king_pca_projection__InferredAncestry.txt")
ancestry <- ancestry[Pr_1st > 0.95]
pheno[ancestry, on = .(eid=IID), ancestry := Anc_1st]

pheno <- pheno[ancestry == ethnicity]

# Load metaGRS and existing T2D PGS
pgs <- fread("data/UKB/T2D_PGS/collated_scores.sscore.gz")
pgs2 <- fread("data/UKB/T2D_PGS_non_European/collated_scores.sscore.gz")
metaGRS <- fread("data/UKB/T2D_metaGRS_levels.txt")
pgs <- pgs[pgs2, on = .(IID)]
pgs <- metaGRS[pgs, on = .(eid=IID)]
pgs[, Gad_metaGRS_v1 := NULL]
pgs[, Gad_metaGRS_v2 := NULL]

# Remove Shim et al. 2023 PGS, whose T2D GWAS summary statistics include non-European UKB participants
pgs[, Shim2023_PGS003867 := NULL]

# Get names of test PRSs
test_pgs <- names(pgs)[-c(1:2)]
test_pgs <- c(setdiff(test_pgs, c("HuertaChagoya2023_EUR_PGS003443", "HuertaChagoya2023_EAS_PGS003444", "HuertaChagoya2023_LAT_PGS003445")), "HuertaChagoya2023_PGS00344X")

# Add to phenotype data
pheno <- pheno[pgs, on = .(eid), nomatch=0]

# Ye2021_PGS001357 has effect allele and other allele incorrectly reported to PGS Catalog;
# as-is the score has strong inverse correlation with T2D status, so we need to flip the 
# sign.
pheno[, Ye2021_PGS001357 := -Ye2021_PGS001357]

# Drop participants free of T2D at baseline who have withdrawn consent for hospital record
# linkage and people with uncertain diabetes status
pheno <- pheno[!is.na(type_2_diabetes) & !is.na(earliest_hospital_date)]

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

# Run the bootstrap analysis in each ancestry separately
boot_res <- foreach(this_anc = unique(pheno$ancestry)) %do% {
  this_pheno <- pheno[ancestry == this_anc]

  # Process risk factors and covariates prior to modelling
  this_pheno[, age := scale(age)]
  this_pheno[, sex := factor(sex, reference="Female")]

  # For the HuertaChagoya2023, we need to combine their three scores as they describe at
  # https://www.pgscatalog.org/score/PGS003443/
  # https://www.pgscatalog.org/score/PGS003444/
  # https://www.pgscatalog.org/score/PGS003445/
  this_pheno[, HuertaChagoya2023_PGS00344X := scale(
    scale(HuertaChagoya2023_EUR_PGS003443)*0.531117 +
    scale(HuertaChagoya2023_EAS_PGS003444)*0.5690198 +
    scale(HuertaChagoya2023_LAT_PGS003445)*0.1465538
  )]
  this_pheno[, HuertaChagoya2023_EUR_PGS003443 := NULL]
  this_pheno[, HuertaChagoya2023_EAS_PGS003444 := NULL]
  this_pheno[, HuertaChagoya2023_LAT_PGS003445 := NULL]

  # Adjust pgs for 20 PCs
  for (this_pgs in c("T2D_metaGRS", test_pgs)) {
    setnames(this_pheno, this_pgs, "this_pgs")
    this_pheno[, this_pgs := scale(lm(scale(this_pgs) ~ PC1 + PC2 + PC3 + PC4 +
      PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 +
      PC16 + PC17 + PC18 + PC19 + PC20)$residuals)]
    setnames(this_pheno, "this_pgs", this_pgs)
  }

	# Filter to relevant columns for prevalent disease analysis
	prev_pheno <- this_pheno[,.SD,.SDcols=c("eid", "ancestry", "age", "sex", "type_2_diabetes", "T2D_metaGRS", test_pgs)]

  # Run bootstrap analysis
  boot(prev_pheno, boot.fun, 1000)  
}
names(boot_res) <- unique(pheno$ancestry)
saveRDS(boot_res, file="output/UKB_tests/metaPRS_delta_auc_bootstraps_non_EUR.rds")

# Curate results
delta_AUC <- foreach(this_anc = unique(pheno$ancestry), .combine=rbind) %do% { 
  data.table(
    ancestry = this_anc, prs = test_pgs, delta_AUC = boot_res[[this_anc]]$t0, boot_SE = apply(boot_res[[this_anc]]$t, 2, sd)
  )
}
delta_AUC[, boot_L95 := delta_AUC - qnorm(1-(0.05/2))*boot_SE]
delta_AUC[, boot_U95 := delta_AUC + qnorm(1-(0.05/2))*boot_SE]
delta_AUC[, boot_pval := pmin(1, pnorm(abs(delta_AUC/boot_SE), lower.tail=FALSE)*2)]
fwrite(delta_AUC, sep="\t", quote=FALSE, file="output/UKB_tests/metaPRS_delta_AUC_bootstraps_non_EUR.txt")

# Now compute delta C-index for incident T2D
pheno <- pheno[!is.na(incident_type_2_diabetes)]

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

# Run the bootstrap analysis in each ancestry separately
boot_res <- foreach(this_anc = unique(pheno$ancestry)) %do% {
  this_pheno <- pheno[ancestry == this_anc]

  # Process risk factors and covariates prior to modelling
  this_pheno[, age := scale(age)]
  this_pheno[, sex := factor(sex, reference="Female")]

  # For the HuertaChagoya2023, we need to combine their three scores as they describe at
  # https://www.pgscatalog.org/score/PGS003443/
  # https://www.pgscatalog.org/score/PGS003444/
  # https://www.pgscatalog.org/score/PGS003445/
  this_pheno[, HuertaChagoya2023_PGS00344X := scale(
    scale(HuertaChagoya2023_EUR_PGS003443)*0.531117 +
    scale(HuertaChagoya2023_EAS_PGS003444)*0.5690198 +
    scale(HuertaChagoya2023_LAT_PGS003445)*0.1465538
  )]
  this_pheno[, HuertaChagoya2023_EUR_PGS003443 := NULL]
  this_pheno[, HuertaChagoya2023_EAS_PGS003444 := NULL]
  this_pheno[, HuertaChagoya2023_LAT_PGS003445 := NULL]

  # Adjust pgs for 20 PCs
  for (this_pgs in c("T2D_metaGRS", test_pgs)) {
    setnames(this_pheno, this_pgs, "this_pgs")
    this_pheno[, this_pgs := scale(lm(scale(this_pgs) ~ PC1 + PC2 + PC3 + PC4 +
      PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 +
      PC16 + PC17 + PC18 + PC19 + PC20)$residuals)]
    setnames(this_pheno, "this_pgs", this_pgs)
  }

  # Filter to relevant columns for prevalent disease analysis
  inci_pheno <- this_pheno[,.SD,.SDcols=c("eid", "ancestry", "age", "sex", "incident_censor_years", "incident_type_2_diabetes", "T2D_metaGRS", test_pgs)]

  # Run bootstrap analysis
  surv_cols_idx <- match(c("incident_censor_years", "incident_type_2_diabetes"), names(inci_pheno))
  censboot(inci_pheno, boot.fun, 1000, index=surv_cols_idx)
}
names(boot_res) <- unique(pheno$ancestry)
saveRDS(boot_res, file="output/UKB_tests/metaPRS_delta_cindices_bootstraps_non_EUR.rds")

# Curate results
delta_cind <- foreach(this_anc = unique(pheno$ancestry), .combine=rbind) %do% { 
  data.table(
    ancestry = this_anc, prs = test_pgs, delta_C = boot_res[[this_anc]]$t0, boot_SE = apply(boot_res[[this_anc]]$t, 2, sd)
  )
}
delta_cind[, boot_L95 := delta_C - qnorm(1-(0.05/2))*boot_SE]
delta_cind[, boot_U95 := delta_C + qnorm(1-(0.05/2))*boot_SE]
delta_cind[, boot_pval := pmin(1, pnorm(abs(delta_C/boot_SE), lower.tail=FALSE)*2)]
fwrite(delta_cind, sep="\t", quote=FALSE, file="output/UKB_tests/metaPRS_delta_cind_bootstraps_non_EUR.txt")

