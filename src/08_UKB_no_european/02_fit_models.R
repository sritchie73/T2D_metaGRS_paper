library(data.table)
library(survival)
library(foreach)
library(RNOmni)
library(ggplot2)
library(hexbin) # required for geom_hex to work
library(patchwork)
source("src/functions/glm_test.R")
source("src/functions/cox_test.R")
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

# Tabulate case/control numbers
cc <- pheno[!is.na(type_2_diabetes) & !is.na(earliest_hospital_nation),
  .(N=.N, prev_T2D=sum(type_2_diabetes, na.rm=TRUE), inci_T2D=sum(incident_type_2_diabetes, na.rm=TRUE)),
            by=.(ancestry)][order(-N)]
fwrite(cc, sep="\t", quote=FALSE, file="output/UKB_tests/UKB_transferrability_case_control_numbers.txt")

# Load metaGRS and existing T2D PGS
pgs <- fread("data/UKB/T2D_PGS/collated_scores.sscore.gz")
pgs2 <- fread("data/UKB/T2D_PGS_non_European/collated_scores.sscore.gz")
metaGRS <- fread("data/UKB/T2D_metaGRS_levels.txt")
pgs <- pgs[pgs2, on = .(IID)]
pgs <- metaGRS[pgs, on = .(eid=IID)]

# Remove Shim et al. 2023 PGS, whose T2D GWAS summary statistics include non-European UKB participants
pgs[, Shim2023_PGS003867 := NULL]

# Add to phenotype data
pheno <- pheno[pgs, on = .(eid), nomatch=0]

# Ye2021_PGS001357 has effect allele and other allele incorrectly reported to PGS Catalog;
# as-is the score has strong inverse correlation with T2D status, so we need to flip the 
# sign.
pheno[, Ye2021_PGS001357 := -Ye2021_PGS001357]

# Examine prevalent T2D in each test grouping
prev <- foreach(this_ancestry = unique(pheno$ancestry), .combine=rbind) %do% {
  # Extract samples of interest
  dat <- pheno[ancestry == this_ancestry]

	# Drop participants free of T2D at baseline who have withdrawn consent for hospital record
	# linkage and people with uncertain diabetes status
	dat <- dat[!is.na(type_2_diabetes) & !is.na(earliest_hospital_nation)]

	# Process risk factors and covariates prior to modelling - we want all
	# estimates to be per SD change in variable, and for factors, using the
	# largest group or lowest risk as reference 
	dat[, age := scale(age)]
	dat[, sex := factor(sex, reference="Female")]
  dat[, assessment_centre := factor(assessment_centre, reference="Leeds")] # For consistency across tests

  # For the HuertaChagoya2023, we need to combine their three scores as they describe at
  # https://www.pgscatalog.org/score/PGS003443/
  # https://www.pgscatalog.org/score/PGS003444/
  # https://www.pgscatalog.org/score/PGS003445/
  dat[, HuertaChagoya2023_PGS00344X := scale(
    scale(HuertaChagoya2023_EUR_PGS003443)*0.531117 +
    scale(HuertaChagoya2023_EAS_PGS003444)*0.5690198 +
    scale(HuertaChagoya2023_LAT_PGS003445)*0.1465538
  )]
  dat[, HuertaChagoya2023_EUR_PGS003443 := NULL]
  dat[, HuertaChagoya2023_EAS_PGS003444 := NULL]
  dat[, HuertaChagoya2023_LAT_PGS003445 := NULL]

	# Adjust pgs for 20 PCs
  pgs_list <- c(intersect(names(pgs)[-1], names(dat)), "HuertaChagoya2023_PGS00344X")
	for (this_pgs in pgs_list) {
		setnames(dat, this_pgs, "this_pgs")
		dat[, this_pgs := scale(lm(scale(this_pgs) ~ PC1 + PC2 + PC3 + PC4 + 
			PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + 
			PC16 + PC17 + PC18 + PC19 + PC20)$residuals)]
		setnames(dat, "this_pgs", this_pgs)
	}

	# Compute null models
	null <- rbind(idcol="model",
		"age + sex"=glm.test(type_2_diabetes ~ age + sex, "type_2_diabetes", dat)
	)

	mf <- "type_2_diabetes ~ %s + age + sex"
	prs <- foreach(this_pgs = pgs_list, .combine=rbind) %do% {
		res <- glm.test(sprintf(mf, this_pgs), "type_2_diabetes", dat)
		res[, model := this_pgs]
	} 

  cbind("ancestry"=this_ancestry, rbind(null, prs))
}

# Rename coefficients for readability
prev[coefficient == "age", coefficient := "Age"]
prev[coefficient == "sexMale", coefficient := "Sex: Male vs. Female"]

# Write out
fwrite(prev, sep="\t", quote=FALSE, file="output/UKB_tests/non_european_prevalent_T2D_associations.txt")

# Examine incident T2D in each test grouping
inci <- foreach(this_ancestry = unique(pheno$ancestry), .combine=rbind) %do% {
  # Extract samples of interest
  dat <- pheno[ancestry == this_ancestry]

  # Drop participants for whom incident diabetes cannot be fit
	dat <- dat[!is.na(incident_type_2_diabetes)]

	# Process risk factors and covariates prior to modelling - we want all
	# estimates to be per SD change in variable, and for factors, using the
	# largest group or lowest risk as reference 
	dat[, age := scale(age)]
	dat[, sex := factor(sex, reference="Female")]

  # For the HuertaChagoya2023, we need to combine their three scores as they describe at
  # https://www.pgscatalog.org/score/PGS003443/
  # https://www.pgscatalog.org/score/PGS003444/
  # https://www.pgscatalog.org/score/PGS003445/
  dat[, HuertaChagoya2023_PGS00344X := scale(
    scale(HuertaChagoya2023_EUR_PGS003443)*0.531117 +
    scale(HuertaChagoya2023_EAS_PGS003444)*0.5690198 +
    scale(HuertaChagoya2023_LAT_PGS003445)*0.1465538
  )]
  dat[, HuertaChagoya2023_EUR_PGS003443 := NULL]
  dat[, HuertaChagoya2023_EAS_PGS003444 := NULL]
  dat[, HuertaChagoya2023_LAT_PGS003445 := NULL]

	# Adjust pgs for 20 PCs
  pgs_list <- c(intersect(names(pgs)[-1], names(dat)), "HuertaChagoya2023_PGS00344X")
  for (this_pgs in pgs_list) {
		setnames(dat, this_pgs, "this_pgs")
		dat[, this_pgs := scale(lm(scale(this_pgs) ~ PC1 + PC2 + PC3 + PC4 + 
			PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + 
			PC16 + PC17 + PC18 + PC19 + PC20)$residuals)]
		setnames(dat, "this_pgs", this_pgs)
	}

  # Fit null model
  null <- rbind(idcol="model", 
    "age + sex"=cox.test(Surv(incident_censor_years, incident_type_2_diabetes) ~ age + sex, "incident_type_2_diabetes", dat)
  )

  # Fit for each PRS
  mf <- "Surv(incident_censor_years, incident_type_2_diabetes) ~ %s + age + sex"
	prs <- foreach(this_pgs = pgs_list, .combine=rbind) %do% {
		res <- cox.test(sprintf(mf, this_pgs), "incident_type_2_diabetes", dat)
		res[, model := this_pgs]
	}

  cbind("ancestry"=this_ancestry, rbind(null, prs))
}

# Rename coefficients for readability
inci[coefficient == "age", coefficient := "Age"]
inci[coefficient == "sexMale", coefficient := "Sex: Male vs. Female"]

# Write out
fwrite(inci, sep="\t", quote=FALSE, file="output/UKB_tests/non_european_incident_T2D_associations.txt")

