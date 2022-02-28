library(data.table)
library(survival)
library(foreach)
source("src/functions/glm_test.R")
source("src/functions/cox_test.R")
source("src/functions/factor_by_size.R")
source("src/functions/define_case_control_status.R")

# Create output directory
system("mkdir -p output/metaGRS/test", wait=TRUE)

# Load dataset for testing metaGRSs
pheno <- fread("data/UKB/collated_curated_data.txt")
pheno <- pheno[(metaGRS_test_samples)]

# curate set of case control definitions and models to test
# 
# Here we drop the combined prevalent + incident diseases definition used for model
# training - we only care how models perform for case/control definitions that we 
# would actually test (prevalent or incident T2D).
case_control_definitions <- case_control_definitions[type != "prevalent + incident"]

# Tabulate case-control numbers for each T2D case/control definition
cc <- foreach(mIdx = case_control_definitions[,.I], .combine=rbind) %do% {
	this_model <- case_control_definitions[mIdx]
	dat <- set_case_control_definition(pheno, this_model$short_name) # Extract case-control definition
	this_model[, N_samples := dat[,.N]]
	this_model[, N_cases := dat[(T2D), .N]]
	this_model[, N_controls := dat[!(T2D), .N]]
	this_model[, N_missing := pheno[,.N] - dat[,.N]]
	this_model
}
fwrite(cc, sep="\t", quote=FALSE, file="output/metaGRS/test/cohort_case_numbers.txt")

# Load all candidate metaGRS and previous T2D PGS for comparison
load_pgs <- function(fname, candidate_metaGRS=TRUE) {
  wide <- fread(fname)
  setnames(wide, "IID", "eid")
  wide <- wide[eid %in% pheno$eid]

  if (candidate_metaGRS) {
    info <- data.table(name=names(wide)[-1])
    info[, lambda := gsub("_", ".", gsub(".*lambda", "lambda", name))]
    info[, alpha := as.numeric(gsub("_", ".", gsub("(.*alpha_)|(_lambda.*)", "", name)))]
    info[, type := gsub("_.*", "", name)]
    info[, prefilter := fcase(
      name %like% "none", "none",     
      name %like% "sig_auc", "sig_auc",
      name %like% "no_auc_ci_overlap", "no_auc_ci_overlap"
    )]
    info[, case_def := gsub(sprintf("(%s_)|(_%s.*)", type, prefilter), "", name), by=.(type, prefilter)]
    info[type == "combined", type := "prevalent + incident"]
  }

  long <- melt(wide, id.vars="eid", variable.name="name", value.name="level")

  if (candidate_metaGRS) {
    long <- long[info, on = .(name)]
    long[, name := "metaGRS"]
  }
  
  return(long)
}

pgs <- rbind(fill=TRUE,
  load_pgs("data/UKB/T2D_PGS/collated_scores.sscore.gz", candidate_metaGRS=FALSE),
  load_pgs("output/metaGRS/all_candidate_metaGRS_lvls/collated_candidate_metaGRSs_prevalent_T2D/collated_scores.sscore.gz"),
  load_pgs("output/metaGRS/all_candidate_metaGRS_lvls/collated_candidate_metaGRSs_incident_T2D/collated_scores.sscore.gz"),
  load_pgs("output/metaGRS/all_candidate_metaGRS_lvls/collated_candidate_metaGRSs_prevalent_plus_incident_T2D/collated_scores.sscore.gz")
)
pgs[, idx := .GRP, by=.(name, lambda, alpha, type, prefilter, case_def)]
invisible(gc())

# Get performance of each candidate metaGRS
for (mIdx in case_control_definitions[,.I]) {
	# Can be run as array job, in which case we will only run the code if the array task
	# matches the current iteration
	if (Sys.getenv("SLURM_ARRAY_TASK_ID") != "" && Sys.getenv("SLURM_ARRAY_TASK_ID") != mIdx) {
		next
	}

	this_model <- case_control_definitions[mIdx]
	cat(this_model$short_name, "\n")
	dat <- set_case_control_definition(pheno, this_model$short_name) # Extract case-control definition

  # Fit null model (no PGS)
  cat("Computing null models...\n")
	if (this_model$type == "prevalent") {
		mf <- T2D ~ genetic_sex + age +
			factor_by_size(earliest_hospital_nation) + # differences in length of retrospective follow-up by nation
			factor_by_size(assessment_centre) # potentially confounded by differences in average date of baseline assessment
		suppressMessages(g1 <- glm.test(formula=mf, event_col="T2D", data=dat, ci.method="wald"))
		null_model <- data.table(model_type=this_model$type, model_case_def=gsub("(prevalent_)|(incident_)", "", this_model$short_name), pgs="none",
      model_fit_metric="AUC", model_fit_estimate=g1$AUC[1], model_fit_L95=g1$AUC.L95[1], model_fit_U95=g1$AUC.U95[1], model_fit_pval=NA)
	} else if (this_model$type == "incident") {
		mf <- Surv(incident_censor_years, T2D) ~ strata(genetic_sex) + age +
			factor_by_size(earliest_hospital_nation) + # differences in length of retrospective follow-up by nation determines dropped prevalent cases
			factor_by_size(censor_hospital_nation) + # differences in follow-up, particularly due to shorter time in Wales
			factor_by_size(assessment_centre) # Differences in average date of baseline assessment by nation potentially confound both
		c1 <- cox.test(formula=mf, event_col="T2D", data=dat)
		null_model <- data.table(model_type=this_model$type, model_case_def=gsub("(prevalent_)|(incident_)", "", this_model$short_name), pgs="none",
      model_fit_metric="C.index", model_fit_estimate=c1$C.index[1], model_fit_L95=c1$C.L95[1], model_fit_U95=c1$C.U95[1], model_fit_pval=NA)
	}


  # Test each PGS for association with T2D
  cat("Evaluating candidate metaGRS...\n")
	pgs_model <- foreach(pIdx = unique(pgs$idx), .combine=rbind) %do% {
    this_pgs <- pgs[idx == pIdx]
    this_pgs_info <- unique(this_pgs[,.(pgs=name, lambda, alpha, type, prefilter, pgs_case_def=case_def)])
    this_model_info <- null_model[, .(model_type, model_case_def)]

	  # Add PGS to dat
		dat[this_pgs, on = .(eid), PGS := i.level]

		# Adjust for 20 PCs and standardise
		dat[, PGS := scale(lm(scale(PGS) ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 +
			PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 +
			PC17 + PC18 + PC19 + PC20)$residuals)]

		# Test
    if (this_model$type == "prevalent") {
			mf <- T2D ~ genetic_sex + PGS + age +
				factor_by_size(earliest_hospital_nation) + # differences in length of retrospective follow-up by nation
				factor_by_size(assessment_centre) # potentially confounded by differences in average date of baseline assessment
			suppressMessages(g1 <- glm.test(formula=mf, event_col="T2D", data=dat, ci.method="wald"))
      g1 <- g1[coefficient == "PGS", .(
        coef_fit_metric="OR", coef_fit_estimate=OR, coef_fit_L95=OR.L95, coef_fit_U95=OR.U95, coef_fit_pval=P.value,
        model_fit_metric="AUC", model_fit_estimate=AUC, model_fit_L95=AUC.L95, model_fit_U95=AUC.U95, model_fit_pval=NA
      )]
      cbind(this_model_info, this_pgs_info, g1)
		} else if (this_model$type == "incident") {
			mf <- Surv(incident_censor_years, T2D) ~ strata(genetic_sex) + PGS + age +
				factor_by_size(earliest_hospital_nation) + # differences in length of retrospective follow-up by nation determines dropped prevalent cases
				factor_by_size(censor_hospital_nation) + # differences in follow-up, particularly due to shorter time in Wales
				factor_by_size(assessment_centre) # Differences in average date of baseline assessment by nation potentially confound both
			c1 <- cox.test(formula=mf, event_col="T2D", data=dat)
      c1 <- c1[coefficient == "PGS", .(
        coef_fit_metric="HR", coef_fit_estimate=HR, coef_fit_L95=L95, coef_fit_U95=U95, coef_fit_pval=Pvalue,
        model_fit_metric="C.index", model_fit_estimate=C.index, model_fit_L95=C.L95, model_fit_U95=C.U95, model_fit_pval=NA
      )]
      cbind(this_model_info, this_pgs_info, c1)
		}
  }

  # Combine and write out
  perf <- rbind(null_model, pgs_model, fill=TRUE)
  fwrite(perf, sep="\t", quote=FALSE, file=sprintf("output/metaGRS/test/%s_model_performance.txt", this_model$short_name))
}
