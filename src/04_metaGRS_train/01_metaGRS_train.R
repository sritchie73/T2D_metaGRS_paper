library(data.table)
library(foreach)
library(caret)
library(survival)
library(glmnet)
library(ggplot2)
source("src/functions/factor_by_size.R")
source("src/functions/define_case_control_status.R")

# Load training data
pheno <- fread("data/UKB/collated_curated_data.txt")
pheno <- pheno[(metaGRS_train_samples)]

# Load list of best models
best_grss <- fread("output/ldpred2/hyperparam_selection/best_models.txt")
best_grss <- best_grss[metric %in% c("AUC", "C.index")]

# Build metaGRS for each case-control definition and model type, with different pre-filtering methods
# on input gwases
iter <- 0
for (mIdx in case_control_definitions[,.I]) {
  for (prefilter in c("none", "sig_auc", "no_auc_ci_overlap")) {
    # Can be run as array job, in which case we will only run the code if the array task 
    # matches the current iteration
    iter <- iter + 1
    if (Sys.getenv("SLURM_ARRAY_TASK_ID") != "" && Sys.getenv("SLURM_ARRAY_TASK_ID") != iter) {
      next
    }

    # Extract details on current case control definition
    this_model <- case_control_definitions[mIdx]

		# Setup output directory
    out_dir <- sprintf("output/metaGRS/train/%s/prefilter_%s/", this_model$short_name, prefilter)
    system(sprintf("mkdir -p %s", out_dir), wait=TRUE)
  
    # Apply case-control definition and subset data
    this_pheno <- set_case_control_definition(pheno, this_model$short_name)

    # Apply pre-filtering to list of GRSs
    # 
    # P-values in 'best_grss' are hardcoded based on the following criteria:
    #
    #  'pval == 0.001': 95% confidence interval for AUC / C-index in model with GRS does not overlap 
    #                   with 95% confidence interval of AUC / C-index for model without GRS
    #  'pval == 0.025': 95% confidence interval for AUC / C-index in model with GRS does not overlap
    #                   with point estimate for AUC / C-index for model without GRS, but does overlap
    #                   with the  95% confidence interval of AUC / C-index for model without GRS
    #    'pval == 0.5': 95% confidence interval for AUC / C-index in model with GRS overlaps the point
    #                   estimate of the AUC / C-index for model without GRS, i.e. this GRS is unlikely
    #                   to add any information to T2D prediction.
    #
    if (prefilter == "none") {
      this_best_grss <- best_grss[case_def == this_model$short_name]
    } else if (prefilter == "sig_auc") {
      this_best_grss <- best_grss[case_def == this_model$short_name & pval < 0.05] 
    } else if (prefilter == "no_auc_ci_overlap") {
      this_best_grss <- best_grss[case_def == this_model$short_name & pval == 0.001]
    }

    # Load GRSs we want to plug into elasticnet
    pgs <- foreach(pIdx = this_best_grss[,.I], .combine=rbind) %do% {
      this_gwas <- this_best_grss[pIdx, gwas]
      this_pgs <- this_best_grss[pIdx, pgs]

      if (this_pgs == "pgs_auto_final2") {
        this_pgs_lvls <- fread(sprintf("output/ldpred2/test/%s/ldpred2_auto_final2_%s.sscore.gz", this_gwas, this_model$short_name), tmpdir="tmp")
      } else {
        this_pgs_lvls <- fread(sprintf("output/ldpred2/all_hyperparam_grs_lvls/%s/collated_scores.sscore.gz", this_gwas),
                               select=c("IID", gsub("^pgs_", "beta_", this_pgs)), tmpdir="tmp")
      }
      setnames(this_pgs_lvls, c("eid", "value"))
      this_pgs_lvls[,.(eid, gwas=this_gwas, value)]
    }

    # Adjust for 20 PCs
    pgs <- pgs[eid %in% this_pheno$eid]
    pcs <- this_pheno[, .(eid, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, 
      PC10, PC11, PC12, PC13, PC14, PC15, PC16, PC17, PC18, PC19, PC20)]
    pgs <- pgs[pcs, on = .(eid)]
    pgs[, value := scale(lm(scale(value) ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + 
      PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20)$residuals),
      by=gwas]
   
    # Cast to wide format
    pgs <- dcast(pgs, eid ~ gwas, value.var="value")

    # Add to phenotype data
    this_pheno <- this_pheno[pgs, on = .(eid)]

    # Split into 10 folds for cross-validation, balancing splits by prevalent and incident T2D, sex, and covariates relating to
    # differences in follow-up time
    this_pheno[, foldgrp := paste(T2D_information, genetic_sex, assessment_centre, earliest_hospital_nation, censor_hospital_nation)]
    this_pheno[, foldid := createFolds(foldgrp, k=10, list=FALSE)]

    # Write out fold ids
    fwrite(this_pheno[, .(eid, foldgrp, foldid)], sep="\t", quote=FALSE, file=sprintf("%s/cross_validation_fold_ids.txt", out_dir))

    # Extract columns we need - differs depending on whether we're training on prevalent T2D, incident, or both combined
    cols <- c("foldid", "T2D", "age", "genetic_sex", "assessment_centre", "earliest_hospital_nation", setdiff(names(pgs), "eid"))
    if (this_model$type != "prevalent") {
      cols <- c(cols, "censor_hospital_nation")
    }
    if (this_model$type == "incident") {
      cols <- c(cols, "incident_censor_years")
    }
    this_pheno <- this_pheno[, .SD, .SDcols=cols]
  
    # Set up factor variables
		this_pheno[, genetic_sex := factor(genetic_sex)]
		this_pheno[, assessment_centre := factor_by_size(assessment_centre)]
		this_pheno[, earliest_hospital_nation := factor_by_size(earliest_hospital_nation)]
    if (this_model$type != "incident") {
      this_pheno[, T2D := factor(T2D, levels=c(FALSE, TRUE))]
    }
    if (this_model$type != "prevalent") {
      this_pheno[, censor_hospital_nation := factor_by_size(censor_hospital_nation)]
    } 

    # Extract x matrix
    if (this_model$type == "incident") {
			x <- this_pheno[,-c("foldid", "T2D", "incident_censor_years")]
    } else {
			x <- this_pheno[,-c("foldid", "T2D")]
    }
    x <- model.matrix(~ 0 + ., x) 

    # Extract y vector(s) and fold-id
    foldid <- this_pheno$foldid
    case_status <- this_pheno$T2D
    if (this_model$type == "incident") {
      follow_time <- this_pheno$incident_censor_years
    }

    # Setup list of alpha mixing parameters to search across (controls balance of ridge vs. lasso)
    alphas <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1) # 0 = ridge, 1 = lasso
    
    # Run elasticnet (logistic or cox depending on case definition)
    if (this_model$type == "incident") {
      enfits <- foreach(alpha=alphas) %do% {
        cv.glmnet(x, Surv(follow_time, case_status), family="cox", foldid=foldid, alpa=alpha, trace.it=TRUE)
      }
    } else {
			enfits <- foreach(alpha=alphas) %do% {
				cv.glmnet(y = case_status, x=x, foldid=foldid, type.measure="auc", family="binomial", alpha=alpha, trace.it=TRUE)
			}
    } 
    names(enfits) <- alphas
    saveRDS(enfits, file=sprintf("%s/enfit.rds", out_dir))

    # build table of AUCs (or cox elasticnet equivalents)
    fitdt <- foreach(idx = seq_along(alphas), .combine=rbind) %do% {
      data.table(alpha=alphas[idx], lambda=enfits[[idx]][["lambda"]],
                 fit_metric=enfits[[idx]][["name"]], 
                 fit_mean=enfits[[idx]][["cvm"]], fit_sd=enfits[[idx]][["cvsd"]],
                 fit_mean_minus_sd=enfits[[idx]][["cvlo"]], fit_mean_plus_sd=enfits[[idx]][["cvup"]],
                 nonzero=enfits[[idx]][["nzero"]],
                 lambda.min=enfits[[idx]][["lambda.min"]],
                 lambda.1se=enfits[[idx]][["lambda.1se"]])
    }
    fwrite(fitdt, sep="\t", quote=FALSE, file=sprintf("%s/enfit_model_fit.txt", out_dir))
    
    # Build table containing information about the best fits in each case:
    bestfit <- fitdt[lambda == lambda.min | lambda == lambda.1se]
    bestfit[lambda == lambda.min, model := "lambda.min"]
    bestfit[lambda == lambda.1se, model := "lambda.1se"]
    bestfit <- bestfit[, .(alpha, lambda, model, nonzero, fit_metric, fit_mean, fit_sd, fit_mean_minus_sd, fit_mean_plus_sd)]
    fwrite(bestfit, sep="\t", quote=FALSE, file=sprintf("%s/enfit_best_model_fit.txt", out_dir))

    # Plot model fits
    g <- ggplot(fitdt) +
      aes(x=log(lambda), y=fit_mean, ymin=fit_mean_minus_sd, ymax=fit_mean_plus_sd,
          fill=factor(alpha), colour=factor(alpha)) +
      geom_ribbon(colour="#00000000", alpha=0.3, show.legend=FALSE) +
      geom_line() +
      geom_vline(data=bestfit, aes(xintercept=log(lambda), colour=factor(alpha), linetype=model), show.legend=FALSE) +
      scale_fill_manual(name="alpha", values=c("0"="#5e4fa2", "0.1"="#3288bd", "0.25"="#66c2a5",
                        "0.5"="#ffff33", "0.75"="#f46d43", "0.9"="#d53e4f", "1"="#9e0142")) +
      scale_colour_manual(name="alpha", values=c("0"="#5e4fa2", "0.1"="#3288bd", "0.25"="#66c2a5",
                        "0.5"="#ffff33", "0.75"="#f46d43", "0.9"="#d53e4f", "1"="#9e0142")) +
      scale_linetype_manual(guide=FALSE, values=c("lambda.min"="dotted", "lambda.1se"="dashed")) +
      xlab("Log(lambda)") + 
      ylab(sprintf("%s (+/- SD)", fitdt$fit_metric[1])) +
      theme_bw()
    ggsave(g, width=10, height=6, units="in", file=sprintf("%s/enfit.png", out_dir)) 

    # Load variant weights so we can build metaGRS
    varweights <- foreach(pIdx = this_best_grss[,.I], .combine=rbind) %do% {
      this_gwas <- this_best_grss[pIdx, gwas]
      this_pgs <- this_best_grss[pIdx, pgs]

      if (this_pgs == "pgs_auto_final2") {
        this_pgs_varweights <- fread(sprintf("output/ldpred2/test/%s/ldpred2_auto_final2_varweights_%s.txt.gz", this_gwas, this_model$short_name), tmpdir="tmp")
      } else {
        this_pgs_varweights <- fread(sprintf("output/ldpred2/train/%s/ldpred2_pgs_varweights.txt.gz", this_gwas), tmpdir="tmp",
                                     select=c("rsid", "chr", "pos", "effect_allele", "other_allele", gsub("^pgs_", "beta_", this_pgs)))
      }
      setnames(this_pgs_varweights, c("rsid", "chr", "pos", "effect_allele", "other_allele", "beta"))
      cbind(gwas=this_gwas, this_pgs_varweights)
    }

    # Collate the coefficients for the best models
    encoef <- foreach(idx = seq_along(enfits), .combine=rbind) %do% {
      lambda.1se <- as.data.table(as.matrix(coef(enfits[[idx]])), keep.rownames="coefficient")
      setnames(lambda.1se, "1", "beta")

      lambda.min <- as.data.table(as.matrix(coef(enfits[[idx]], s="lambda.min")), keep.rownames="coefficient")
      setnames(lambda.min, "1", "beta")

      dt <- rbind(lambda.min=lambda.min, lambda.1se=lambda.1se, idcol="best.model")
      dt <- cbind(alpha=alphas[idx], dt)
      dt
    }
    fwrite(encoef, sep="\t", quote=FALSE, file=sprintf("%s/enfit_best_coefs.txt", out_dir))

    # Multiply out weights by log OR in elasticnet model with best AUC
    pgs_coef <- encoef[coefficient %in% names(pgs)]
    for (this_alpha in alphas) {
      for (this_model in c("lambda.min", "lambda.1se")) {
        pgs_name <- sprintf("alpha.%s_%s", this_alpha, this_model)
        betas <- pgs_coef[alpha == this_alpha & best.model == this_model]
        varweights[betas, on = .(gwas=coefficient), c(pgs_name) := beta * i.beta] # first beta is from GWAS, second is coefficient from elasticnet
      }
    }
    varweights[, beta := NULL]
    metaGRS_weights <- melt(varweights, id.vars=c("rsid", "chr", "pos", "effect_allele", "other_allele", "gwas"), variable.name="metaGRS", value.name="weight")
    metaGRS_weights <- metaGRS_weights[, .(weight=sum(weight)), by=.(rsid, chr, pos, effect_allele, other_allele, metaGRS)]
    metaGRS_weights <- dcast(metaGRS_weights, rsid + chr + pos + effect_allele + other_allele ~ metaGRS, value.var="weight")
    metaGRS_weights <- metaGRS_weights[order(pos)][order(chr)]
    fwrite(metaGRS_weights, sep="\t", quote=FALSE, compress="gzip", file=sprintf("%s/candidate_metaGRS.txt.gz", out_dir))
  }
}

