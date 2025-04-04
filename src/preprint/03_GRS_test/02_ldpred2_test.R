library(data.table)
library(survival)
library(foreach)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(bigsnpr)
source("src/functions/glm_test.R")
source("src/functions/cox_test.R")
source("src/functions/factor_by_size.R")
source("src/functions/define_case_control_status.R")

##############################################
# Determine which GWAS(s) we're working with
##############################################
if (Sys.getenv("SLURM_ARRAY_TASK_ID") != "") {
  gwass <- list.dirs("data/filtered_sumstats", recursive=FALSE, full.names=FALSE)
  gwass <- gwass[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))]
} else if (!exists("gwass") && length(commandArgs(trailingOnly=TRUE)) == 0) {
  gwass <- list.dirs("data/filtered_sumstats", recursive=FALSE, full.names=FALSE)
} else if (!exists("gwass")) {
  gwass <- commandArgs(trailingOnly=TRUE)
}

# Create output directories
for (gwas in gwass) {
  system(sprintf("mkdir -p output/ldpred2/test/%s", gwas), wait=TRUE)
}

# First, lets sanity check the PGS levels we've compute from the saved betas
# match the PGS levels computed in the LDpred2 training script. Note: some 
# differences expected since LDpred2 usage hard call genotypes whereas now
# we've computed PGS levels from probabilistic dosages.
for (gwas in gwass) {
  # Load PGS levels in training data
  train <- fread(sprintf("output/ldpred2/train/%s/traininset_pgs_levels.txt.gz", gwas), tmpdir="tmp")

  # Load PGS levels computed de-novo
  denovo <- fread(sprintf("output/ldpred2/all_hyperparam_grs_lvls/%s/collated_scores.sscore.gz", gwas), tmpdir="tmp")
  setnames(denovo, gsub("beta", "pgs", names(denovo)))

  # Filter to training samples
  denovo <- denovo[IID %in% train$eid]

  train <- melt(train, id.vars="eid", variable.name="score")
  denovo <- melt(denovo, id.vars="IID", variable.name="score")
  comp <- merge(train, denovo, by.x=c("eid", "score"), by.y=c("IID", "score"), suffixes=c(".train", ".denovo"))

  # Make plot
  nscores <- comp[,length(unique(score))]
  g <- ggplot(comp, aes(x=value.train, y=value.denovo)) +
    geom_point(shape=19, size=0.5, alpha=0.5) +
    geom_abline(intercept=0, slope=1, linetype=2, color="red", size=0.5) + 
    facet_wrap(~ score, scales="free", ncol=ceiling(sqrt(nscores))) + 
    xlab("PGS levels in LDpred2 training") +
    ylab("PGS levels computed de-novo from saved betas") +
    theme_bw() + 
    theme(strip.background=element_blank(), strip.text=element_text(face=2, size=6),
          axis.text=element_text(size=6), axis.title=element_text(size=10))
  ggsave(g, width=ceiling(sqrt(nscores)), height=ceiling(sqrt(nscores)), units="in", dpi=100,
         file=sprintf("output/ldpred2/test/%s/pgs_levels_sanity_check.png", gwas))
}

# Load dataset for testing ldpred2 GRSs
pheno <- fread("data/UKB/collated_curated_data.txt")
pheno <- pheno[(metaGRS_train_samples)]

# curate set of case control definitions and models to test
#
# Here we drop QDiabetes models - these use the same definition as 
# short_name == "incident_basic_adjudicated" but additionally drop samples with
# missing data in the given QDiabetes column - i.e. used later only for testing
# metaGRS
case_control_definitions <- case_control_definitions[!(short_name %like% "QDiabetes")]

# Test each model
for (gwas in gwass) {
  cat(gwas, "\n")
  gc()

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
  fwrite(cc, sep="\t", quote=FALSE, file=sprintf("output/ldpred2/test/%s/cohort_case_numbers.txt", gwas))

  # Fit null model (no PGS) in the test data
  cat("Computing null models...\n")
  null_model <- foreach(mIdx = case_control_definitions[,.I], .combine=rbind) %do% {
    this_model <- case_control_definitions[mIdx]
    cat(this_model$short_name, "\n")
    dat <- set_case_control_definition(pheno, this_model$short_name) # Extract case-control definition
    if (this_model$type == "prevalent") {
      mf <- T2D ~ genetic_sex + age + 
        factor_by_size(earliest_hospital_nation) + # differences in length of retrospective follow-up by nation
        factor_by_size(assessment_centre) # potentially confounded by differences in average date of baseline assessment
      suppressMessages(g1 <- glm.test(formula=mf, event_col="T2D", data=dat, ci.method="wald"))
      data.table(pgs="none", type=this_model$type, case_def=this_model$short_name, metric="AUC", 
        estimate=g1$AUC[1], L95=g1$AUC.L95[1], U95=g1$AUC.U95[1], pval=NA)
    } else if (this_model$type == "incident") {
      mf <- Surv(incident_censor_years, T2D) ~ strata(genetic_sex) + age +
        factor_by_size(earliest_hospital_nation) + # differences in length of retrospective follow-up by nation determines dropped prevalent cases
        factor_by_size(censor_hospital_nation) + # differences in follow-up, particularly due to shorter time in Wales
        factor_by_size(assessment_centre) # Differences in average date of baseline assessment by nation potentially confound both
      c1 <- cox.test(formula=mf, event_col="T2D", data=dat)
      data.table(pgs="none", type=this_model$type, case_def=this_model$short_name, metric="C.index", 
        estimate=c1$C.index[1], L95=c1$C.L95[1], U95=c1$C.U95[1], pval=NA)
    } else if (this_model$type == "prevalent + incident") {
      mf <- T2D ~ genetic_sex + age +
        factor_by_size(earliest_hospital_nation) + # differences in length of retrospective follow-up by nation determines dropped prevalent cases
        factor_by_size(censor_hospital_nation) + # differences in follow-up, particularly due to shorter time in Wales
        factor_by_size(assessment_centre) # Differences in average date of baseline assessment by nation potentially confound both
      suppressMessages(g1 <- glm.test(formula=mf, event_col="T2D", data=dat, ci.method="wald"))
      data.table(pgs="none", type=this_model$type, case_def=this_model$short_name, metric="AUC", 
        estimate=g1$AUC[1], L95=g1$AUC.L95[1], U95=g1$AUC.U95[1], pval=NA)
    }
  }

  # Get performance of each candidate PGS
  cat("Evaluating candidate PGS...\n")
  prs_model <- foreach(mIdx = case_control_definitions[,.I], .combine=rbind) %do% {
    this_model <- case_control_definitions[mIdx]
    cat(this_model$short_name, "\n")
    dat <- set_case_control_definition(pheno, this_model$short_name) # Extract case-control definition

		# Load computed PGSs
		pgs <- fread(sprintf("output/ldpred2/all_hyperparam_grs_lvls/%s/collated_scores.sscore.gz", gwas), tmpdir="tmp")
		setnames(pgs, gsub("beta", "pgs", names(pgs)))

		# Filter to test cohort samples
		pgs <- pgs[IID %in% dat$eid]

		# identify "bad" chains for auto model - must be done prior to correction for PCs
		pgs <- melt(pgs, id.vars="IID", variable.name="score")
		auto_chains <- pgs[score %like% "pgs_auto_[0-9]*$",.(sc=sd(value)),by=score]
		auto_chains[, keep := abs(sc - median(sc)) < 3 * mad(sc)]

    if (auto_chains[(keep), .N] > 0) {
			# Determine new PGS auto final model
			pgs_auto_final2 <- pgs[score %in% auto_chains[(keep), score], .(score="pgs_auto_final2", value=mean(value)), by=IID]
      fwrite(dcast(pgs_auto_final2, IID ~ score, value.var="value"), sep="\t", quote=FALSE, compress="gzip", 
             file=sprintf("output/ldpred2/test/%s/ldpred2_auto_final2_%s.sscore.gz", gwas, this_model$short_name))
			pgs <- rbind(pgs, pgs_auto_final2)

			# Get variant weights for new auto final score
			varweights <- fread(sprintf("output/ldpred2/train/%s/ldpred2_pgs_varweights.txt.gz", gwas), tmpdir="tmp",
													select=c("rsid", "chr", "pos", "effect_allele", "other_allele", auto_chains[(keep), gsub("pgs", "beta", score)]))
			varweights <- melt(varweights, id.vars=c("rsid", "chr", "pos", "effect_allele", "other_allele"))
			varweights <- varweights[, .(beta_auto_final2=mean(value)), by=c("rsid", "chr", "pos", "effect_allele", "other_allele")]
			fwrite(varweights, sep="\t", quote=FALSE, compress="gzip", file=sprintf("output/ldpred2/test/%s/ldpred2_auto_final2_varweights_%s.txt.gz", gwas, this_model$short_name))

      # Compare final auto PGS derived from training set to that derived in test set
      comp <- merge(pgs[score == "pgs_auto_final"], pgs[score == "pgs_auto_final2"], by="IID", suffixes=c(".final", ".final2"))

      g <- ggplot(comp, aes(x=value.final, y=value.final2)) +
        geom_point(shape=19, alpha=0.5) +
        geom_abline(intercept=0, slope=1, linetype=2, color="red") +
        xlab("LDpred2 auto final PGS\nderived from good chains in training dataset") +
        ylab("LDpred2 auto final PGS\nderived from good chains in test dataset") +
        theme_bw()
      ggsave(g, width=7.2, height=7.2, units="in", file=sprintf("output/ldpred2/test/%s/ldpred2_auto_compare_in_%s.png", gwas, this_model$short_name))
    }

    # Test each candidate PGS for association with T2D
    foreach(candidate_pgs = unique(pgs$score), .combine=rbind) %do% { 
      # Add PGS to dat
      dat[pgs[score == candidate_pgs], on = .(eid=IID), PGS := i.value]

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
				g1 <- rbind(idcol="metric",
					"AUC"=g1[1, .(estimate=AUC, L95=AUC.L95, U95=AUC.U95, pval=NA)],
					"OR"=g1[coefficient == "PGS", .(estimate=OR, L95=OR.L95, U95=OR.U95, pval=P.value)]
				)
				cbind(pgs=candidate_pgs, type=this_model$type, case_def=this_model$short_name, g1)
			} else if (this_model$type == "incident") {
				mf <- Surv(incident_censor_years, T2D) ~ strata(genetic_sex) + PGS + age +
					factor_by_size(earliest_hospital_nation) + # differences in length of retrospective follow-up by nation determines dropped prevalent cases
					factor_by_size(censor_hospital_nation) + # differences in follow-up, particularly due to shorter time in Wales
					factor_by_size(assessment_centre) # Differences in average date of baseline assessment by nation potentially confound both
				c1 <- cox.test(formula=mf, event_col="T2D", data=dat)
				c1 <- rbind(idcol="metric",
					"C.index"=c1[1, .(estimate=C.index, L95=C.L95, U95=C.U95, pval=NA)],
					"HR"=c1[coefficient == "PGS", .(estimate=HR, L95=L95, U95=U95, pval=Pvalue)]
				)
				cbind(pgs=candidate_pgs, type=this_model$type, case_def=this_model$short_name, c1)
			} else if (this_model$type == "prevalent + incident") {
				mf <- T2D ~ genetic_sex + PGS + age +
					factor_by_size(earliest_hospital_nation) + # differences in length of retrospective follow-up by nation determines dropped prevalent cases
					factor_by_size(censor_hospital_nation) + # differences in follow-up, particularly due to shorter time in Wales
					factor_by_size(assessment_centre) # Differences in average date of baseline assessment by nation potentially confound both
				suppressMessages(g1 <- glm.test(formula=mf, event_col="T2D", data=dat, ci.method="wald"))
				g1 <- rbind(idcol="metric",
					"AUC"=g1[1, .(estimate=AUC, L95=AUC.L95, U95=AUC.U95, pval=NA)],
					"OR"=g1[coefficient == "PGS", .(estimate=OR, L95=OR.L95, U95=OR.U95, pval=P.value)]
				)
				cbind(pgs=candidate_pgs, type=this_model$type, case_def=this_model$short_name, g1)
			}
    } 
  }

  # Combine PRS and NULL model performance
  perf <- rbind(null_model, prs_model, fill=TRUE)
  fwrite(perf, sep="\t", quote=FALSE, file=sprintf("output/ldpred2/test/%s/all_model_performance.txt", gwas))

  # Generate diagnostic plots
  cat("Generating diagnositic plots...\n")

  # Load LDpred2 grid model parameters
  ldsc <- readRDS(sprintf("output/ldpred2/train/%s/ldsc_results.rds", gwas))
	h2_est <- ldsc[["h2"]]
	h2_seq <- round(h2_est * c(0.3, 0.7, 1, 1.4), 4)
	p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2)
	params <- expand.grid(grid_param_p = p_seq, grid_param_h2 = h2_seq, grid_param_sparse = c(FALSE, TRUE))
  setDT(params)
  params[, paramset := .I]

  # Plot performance against grid parameters
  for (mIdx in case_control_definitions[,.I]) {
		this_model <- case_control_definitions[mIdx]
		this_perf <- perf[pgs %like% "_grid_" & case_def == this_model$short_name]
    this_perf[, paramset := as.integer(gsub("pgs_grid_", "", pgs))]
    this_perf <- merge(this_perf, params, by="paramset", all.x=TRUE)
		this_null <- perf[pgs == "none" & case_def == this_model$short_name]

		g1 <- ggplot(this_perf[metric %in% c("AUC", "C.index")]) + 
			aes(x=grid_param_p, y=estimate, ymin=L95, ymax=U95, color=as.factor(grid_param_h2)) +
			theme_bigstatsr() +
			geom_rect(inherit.aes=FALSE, ymin=this_null$L95, ymax=this_null$U95, xmin=-Inf, xmax=Inf, fill="#bdbdbd", alpha=0.3) +
			geom_hline(yintercept=this_null$estimate, linetype=2) +
			geom_errorbar(width=0, alpha=0.5) +
			geom_point() +
      geom_line() + 
			scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
			facet_wrap(~ grid_param_sparse, labeller = label_both) +
			labs(y = sprintf("%s (+/- 95%% CI) for %s T2D", this_null$metric, this_model$type), color = "h2") +
			theme_bw() +
			theme(legend.position = "top", panel.spacing = unit(1, "lines"))

		g2 <- ggplot(this_perf[metric %in% c("OR", "HR")]) +
			aes(x=grid_param_p, y=estimate, ymin=L95, ymax=U95, color=as.factor(grid_param_h2)) +
			theme_bigstatsr() +
			geom_hline(yintercept=1, linetype=2) +
			geom_errorbar(width=0, alpha=0.5) +
			geom_point() +
      geom_line() +
			scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
			facet_wrap(~ grid_param_sparse, labeller = label_both) +
			labs(y = sprintf("%s (+/- 95%% CI) for %s T2D", this_perf[metric %in% c("OR", "HR"), metric][1], this_model$type), color = "h2") +
			theme_bw() +
			theme(legend.position = "top", panel.spacing = unit(1, "lines"))

		g <- plot_grid(g1, g2, nrow = 2)
		ggsave(g, width=7.2, height=7.2, units="in", file=sprintf("output/ldpred2/test/%s/ldpred2_grid_performance_%s.png", gwas, this_model$short_name))
  }

	# Get lasso parameters
	params <- fread(sprintf("output/ldpred2/train/%s/ldpred2_lassosum2_parameters.txt", gwas))

  # Plot performance against lasso parameters
  for (mIdx in case_control_definitions[,.I]) {
		this_model <- case_control_definitions[mIdx]
		this_perf <- perf[pgs %like% "_lassosum2_" & case_def == this_model$short_name]
		this_perf[, paramset := as.integer(gsub("pgs_lassosum2_", "", pgs))]
    this_perf <- merge(this_perf, params, by="paramset", all.x=TRUE)
		this_null <- perf[pgs == "none" & case_def == this_model$short_name]

		g1 <- ggplot(this_perf[metric %in% c("AUC", "C.index")]) + 
			aes(x=lambda, y=estimate, ymin=L95, ymax=U95, color=as.factor(delta)) +
			theme_bigstatsr() +
			geom_rect(inherit.aes=FALSE, ymin=this_null$L95, ymax=this_null$U95, xmin=-Inf, xmax=Inf, fill="#bdbdbd", alpha=0.3) +
			geom_hline(yintercept=this_null$estimate, linetype=2) +
			geom_errorbar(width=0, alpha=0.5) +
			geom_point() +
      geom_line() +
			scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
			labs(y = sprintf("%s (+/- 95%% CI) for %s T2D", this_null$metric, this_model$type), color = "delta") +
			theme_bw() +
			theme(legend.position = "top", panel.spacing = unit(1, "lines"))

		g2 <- ggplot(this_perf[metric %in% c("OR", "HR")]) +
			aes(x=lambda, y=estimate, ymin=L95, ymax=U95, color=as.factor(delta)) +
			theme_bigstatsr() +
			geom_hline(yintercept=1, linetype=2) +
			geom_errorbar(width=0, alpha=0.5) +
			geom_point() +
      geom_line() +
			scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
			labs(y = sprintf("%s (+/- 95%% CI) for %s T2D", this_perf[metric %in% c("OR", "HR"), metric][1], this_model$type), color = "delta") +
			theme_bw() +
			theme(legend.position = "top", panel.spacing = unit(1, "lines"))

		g <- plot_grid(g1, g2, nrow = 2)
		ggsave(g, width=7.2, height=7.2, units="in", file=sprintf("output/ldpred2/test/%s/ldpred2_lassosum2_performance_%s.png", gwas, this_model$short_name))
  }
}
