library(data.table)
library(survival)
library(foreach)
library(RNOmni)
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

# Define broad ancestry groups
pheno[, ethnicity_grouping := fcase(
  ethnicity_subgroup %in% c("African", "Caribbean", "Black or Black British", "Any other Black background"), "African",
  ethnicity_subgroup %in% c("Bangladeshi", "Indian", "Pakistani"), "South Asian",
  ethnicity_subgroup %in% c("Chinese", "Asian or Asian British", "Any other Asian background"), "East Asian",
  ethnicity_subgroup %in% c("White", "British", "Irish", "Any other white background"), "European",
  ethnicity_subgroup %in% c("Mixed", "White and Asian", "White and Black Caribbean", "White and Black African", "Any other mixed background"), "Mixed",
  ethnicity_subgroup == "Other ethnic group", "Other",
  ethnicity_subgroup %in% c("Prefer not to answer", "Do not know"), NA_character_
)]

# Tabulate case/control numbers
cc <- pheno[,.(N=.N, prev_T2D=sum(type_2_diabetes, na.rm=TRUE), inci_T2D=sum(incident_type_2_diabetes, na.rm=TRUE)),
            by=.(ethnicity_subgroup, ethnicity_grouping)][order(-N)][order(ethnicity_grouping)]
fwrite(cc, sep="\t", quote=FALSE, file="output/UKB_tests/UKB_transferrability_case_control_numbers.txt")

# Set up information about groupings to test
test_groups <- rbind(
  data.table(grouping="African", group_column="ethnicity_grouping"),
  data.table(grouping="South Asian", group_column="ethnicity_grouping"),
  data.table(grouping="East Asian", group_column="ethnicity_grouping"),
  data.table(grouping="OtherAsian", group_column="ethnicity_qdiabetes_groups"),
  data.table(grouping="Other", group_column="ethnicity_qdiabetes_groups"),
  data.table(grouping="African", group_column="ethnicity_subgroup"),
  data.table(grouping="Caribbean", group_column="ethnicity_subgroup"),
  data.table(grouping="Indian", group_column="ethnicity_subgroup"),
  data.table(grouping="Pakistani", group_column="ethnicity_subgroup"),
  data.table(grouping="Bangladeshi", group_column="ethnicity_subgroup"),
  data.table(grouping="Chinese", group_column="ethnicity_subgroup"),
  data.table(grouping="Any other Asian background", group_column="ethnicity_subgroup"),
  data.table(grouping="Other ethnic group", group_column="ethnicity_subgroup")
)

# Load metaGRS and existing T2D PGS
pgs <- fread("data/UKB/T2D_PGS/collated_scores.sscore.gz")
pgs2 <- fread("data/UKB/T2D_PGS_non_European/collated_scores.sscore.gz")
metaGRS <- fread("data/UKB/T2D_metaGRS_levels.txt")
pgs <- pgs[pgs2, on = .(IID)]
pgs <- metaGRS[pgs, on = .(eid=IID)]

# Add to phenotype data
pheno <- pheno[pgs, on = .(eid), nomatch=0]

# Ye2021_PGS001357 has effect allele and other allele incorrectly reported to PGS Catalog;
# as-is the score has strong inverse correlation with T2D status, so we need to flip the 
# sign.
pheno[, Ye2021_PGS001357 := -Ye2021_PGS001357]

# Examine prevalent T2D in each test grouping
prev <- foreach(gIdx = test_groups[,.I], .combine=rbind) %do% {
  # Extract samples of interest
  this_group <- test_groups[gIdx]
  if (this_group$group_column == "ethnicity_grouping") {
    dat <- pheno[ethnicity_grouping == this_group$grouping]
  } else if (this_group$group_column == "ethnicity_qdiabetes_groups") {
    dat <- pheno[ethnicity_qdiabetes_groups == this_group$grouping]
  } else if (this_group$group_column == "ethnicity_subgroup") {
    dat <- pheno[ethnicity_subgroup == this_group$grouping]
  }

	# Drop participants free of T2D at baseline who have withdrawn consent for hospital record
	# linkage and people with uncertain diabetes status
	dat <- dat[!is.na(type_2_diabetes) & !is.na(earliest_hospital_nation)]

	# Process risk factors and covariates prior to modelling - we want all
	# estimates to be per SD change in variable, and for factors, using the
	# largest group or lowest risk as reference 
	dat[, age := scale(age)]
	dat[, sex := factor(sex, reference="Female")]
  dat[, assessment_centre := factor(assessment_centre, reference="Leeds")] # For consistency across tests

	# Adjust pgs for 20 PCs
	for (this_pgs in names(pgs)[-1]) {
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
	prs <- foreach(this_pgs = names(pgs)[-1], .combine=rbind) %do% {
		res <- glm.test(sprintf(mf, this_pgs), "type_2_diabetes", dat)
		res[, model := this_pgs]
	} 

  cbind(this_group, rbind(null, prs))
}

# Rename coefficients for readability
prev[coefficient == "age", coefficient := "Age"]
prev[coefficient == "sexMale", coefficient := "Sex: Male vs. Female"]

# Write out
fwrite(prev, sep="\t", quote=FALSE, file="output/UKB_tests/non_european_prevalent_T2D_associations.txt")

# Examine incident T2D in each test grouping
inci <- foreach(gIdx = test_groups[,.I], .combine=rbind) %do% {
  # Extract samples of interest
  this_group <- test_groups[gIdx]
  if (this_group$group_column == "ethnicity_grouping") {
    dat <- pheno[ethnicity_grouping == this_group$grouping]
  } else if (this_group$group_column == "ethnicity_qdiabetes_groups") {
    dat <- pheno[ethnicity_qdiabetes_groups == this_group$grouping]
  } else if (this_group$group_column == "ethnicity_subgroup") {
    dat <- pheno[ethnicity_subgroup == this_group$grouping]
  }

  # Drop participants for whom incident diabetes cannot be fit
	dat <- dat[!is.na(incident_type_2_diabetes)]

	# Process risk factors and covariates prior to modelling - we want all
	# estimates to be per SD change in variable, and for factors, using the
	# largest group or lowest risk as reference 
	dat[, age := scale(age)]
	dat[, sex := factor(sex, reference="Female")]

	# Adjust pgs for 20 PCs
	for (this_pgs in names(pgs)[-1]) {
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
	prs <- foreach(this_pgs = names(pgs)[-1], .combine=rbind) %do% {
		res <- cox.test(sprintf(mf, this_pgs), "incident_type_2_diabetes", dat)
		res[, model := this_pgs]
	}

  cbind(this_group, rbind(null, prs))
}

# Rename coefficients for readability
inci[coefficient == "age", coefficient := "Age"]
inci[coefficient == "sexMale", coefficient := "Sex: Male vs. Female"]

# Write out
fwrite(inci, sep="\t", quote=FALSE, file="output/UKB_tests/non_european_incident_T2D_associations.txt")

