library(data.table)
library(foreach)
library(ggplot2)
library(RColorBrewer)

# create output directory
system("mkdir -p output/metaGRS/hyperparam_selection", wait=TRUE)

# Load in candidate metaGRS performance in training cross-validation
case_defs <- list.dirs("output/metaGRS/train/", recursive=FALSE, full.names=FALSE)
cvfit <- foreach(this_case_def = case_defs, .combine=rbind) %do% {
  prefilters <- list.dirs(sprintf("output/metaGRS/train/%s", this_case_def), recursive=FALSE, full.names=FALSE)
  foreach(this_prefilter = prefilters, .combine=rbind) %do% {
    this_perf <- fread(sprintf("output/metaGRS/train/%s/%s/enfit_best_model_fit.txt", this_case_def, this_prefilter))
    setnames(this_perf, gsub("fit_", "cvfit_", names(this_perf)))
    this_perf[, .(
      type = fcase(
        this_case_def %like% "^prevalent", "prevalent",
        this_case_def %like% "^incident", "incident",
        this_case_def %like% "^combined", "prevalent + incident"
      ),
      case_def = gsub("(^prevalent_)|(^incident_)|(^combined_)", "", this_case_def),
      prefilter = gsub("(^prefilter_)", "", this_prefilter),
      alpha, lambda=model, cvfit_metric, cvfit_mean, cvfit_sd, cvfit_mean_minus_sd, cvfit_mean_plus_sd
    )]
  }
}

# Load in candidate metaGRS (and other existing T2D PGS) performance in training set
train_perf_files <- list.files(path="output/metaGRS/test/trainingset", pattern=".*model_performance.txt", full.names=TRUE)
train_perf <- rbindlist(lapply(train_perf_files, fread, tmpdir="tmp"))

# Load in candidate metaGRS (and other existing T2D PGS) performance in test set
test_perf_files <- list.files(path="output/metaGRS/test/testset", pattern=".*model_performance.txt", full.names=TRUE)
test_perf <- rbindlist(lapply(test_perf_files, fread, tmpdir="tmp"))

# Here, we want to find the most principled approach to constructing a metaGRS, however, there are multiple
# justifiable ways to do this:
#
# - There are multiple valid ways of classifying controls (non-T2D cases)
# - There are two equally valid models output by elasticnet ('lambda.min' or 'lambda.1se') that loosely 
#   correspond to "assume no overfitting" or "robust to overfitting at the expense of power".
# - We may wish to prefilter the list of GWAS/GRSs input to elasticnet to reduce overfitting
# - We may wish to train against incident T2D instead of prevalent, or combine prevalent and incident cases.

# First, extract the optimal 'lambda.min' and 'lambda.1se' models for each elasticnet (here we ran
# with seven alpha mixing parameters between 0 and 1 corresponding to mixture of ridge and lasso penalties)
cvfit <- cvfit[, .SD[which.max(cvfit_mean)], by=.(type, case_def, prefilter, lambda)]

# Split out model performance:
train_null_perf <- train_perf[pgs == "none"]
train_metaGRS_perf <- train_perf[pgs == "metaGRS"]
train_metaGRS_perf <- train_metaGRS_perf[
  cvfit[,.(type, case_def, prefilter, lambda, alpha)], 
  on = .(type, pgs_case_def=case_def, prefilter, lambda, alpha)
]

test_null_perf <- test_perf[pgs == "none"]
test_pgs_perf <- test_perf[pgs != "none" & pgs != "metaGRS"]
test_metaGRS_perf <- test_perf[pgs == "metaGRS"]
test_metaGRS_perf <- test_metaGRS_perf[
  cvfit[,.(type, case_def, prefilter, lambda, alpha)], 
  on = .(type, pgs_case_def=case_def, prefilter, lambda, alpha)
]

# First, we ask:
#
#  (1) What method of classifying controls leads to the best prediction across possible case/control definitions
#      and prevalent/incident disease prediction, and
#  (2) What method of classifying cases (prevalent, incident, or both combined) leads to best prediction across 
#      possible case/control definitions and prevalent/incident disease prediction
#
# For this, we just look at the performance of various candidate metaGRS across these different scenarios in
# the training data:

g <- ggplot(train_metaGRS_perf) +
  aes(y=model_fit_estimate, ymin=model_fit_L95, ymax=model_fit_U95,
      x=paste(type, pgs_case_def), color=lambda, shape=prefilter) +
  geom_rect(inherit.aes=FALSE, data=train_null_perf,
            aes(ymin=model_fit_L95, ymax=model_fit_U95), 
            xmin=-Inf, xmax=Inf, fill="#bdbdbd", alpha=0.3) +
  geom_hline(data=train_null_perf, aes(yintercept=model_fit_estimate), linetype=2) +
  geom_errorbar(width=0, alpha=0.5, position=position_dodge(width=0.9)) +
  geom_point(position=position_dodge(width=0.9)) +
  scale_colour_manual(name="elasticnet best model", 
                      values=c("lambda.min"="#762a83", "lambda.1se"="#1b7837")) +
  scale_shape_manual(name="GRS prefiltering",
                     values=c("none"=19, "sig_auc"=18, "no_auc_ci_overlap"=15)) +
  facet_wrap(~ model_type + model_case_def, scales="free_y", ncol=3) +
  xlab("") + 
  ylab("AUC or C-index (incident) (+/- 95% CI) for T2D") +
	theme_bw() +
	theme(
		legend.position="bottom",
		panel.grid.major.x=element_blank(),
		panel.grid.minor.x=element_blank(),
		axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=0.5),
		axis.text.y=element_text(size=8),
		axis.title=element_text(size=10),
		strip.text.x=element_text(size=8),
		strip.text.y=element_text(size=8),
    legend.title=element_text(size=8),
    legend.text=element_text(size=8)
	)
ggsave(g, width=30, height=30, units="in", file="output/metaGRS/hyperparam_selection/case_def_selection_trainingset.png")

g <- ggplot(test_metaGRS_perf) +
  aes(y=model_fit_estimate, ymin=model_fit_L95, ymax=model_fit_U95,
      x=paste(type, pgs_case_def), color=lambda, shape=prefilter) +
  geom_rect(inherit.aes=FALSE, data=test_null_perf,
            aes(ymin=model_fit_L95, ymax=model_fit_U95), 
            xmin=-Inf, xmax=Inf, fill="#bdbdbd", alpha=0.3) +
  geom_hline(data=test_null_perf, aes(yintercept=model_fit_estimate), linetype=2) +
  geom_errorbar(width=0, alpha=0.5, position=position_dodge(width=0.9)) +
  geom_point(position=position_dodge(width=0.9)) +
  scale_colour_manual(name="elasticnet best model", 
                      values=c("lambda.min"="#762a83", "lambda.1se"="#1b7837")) +
  scale_shape_manual(name="GRS prefiltering",
                     values=c("none"=19, "sig_auc"=18, "no_auc_ci_overlap"=15)) +
  facet_wrap(~ model_type + model_case_def, scales="free_y", ncol=3) +
  xlab("") + 
  ylab("AUC or C-index (incident) (+/- 95% CI) for T2D") +
	theme_bw() +
	theme(
		legend.position="bottom",
		panel.grid.major.x=element_blank(),
		panel.grid.minor.x=element_blank(),
		axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=0.5),
		axis.text.y=element_text(size=8),
		axis.title=element_text(size=10),
		strip.text.x=element_text(size=8),
		strip.text.y=element_text(size=8),
    legend.title=element_text(size=8),
    legend.text=element_text(size=8)
	)
ggsave(g, width=30, height=30, units="in", file="output/metaGRS/hyperparam_selection/case_def_selection_testset.png")

# We see training on prevalent + incident T2D combined gives best performance overall, so lets restrict to those metaGRS
cvfit <- cvfit[type == "prevalent + incident"]
train_metaGRS_perf <- train_metaGRS_perf[type == "prevalent + incident"]
test_metaGRS_perf <- test_metaGRS_perf[type == "prevalent + incident"]

# Now lets see if we see a difference in performance based on the way controls are defined. For this we need to zoom in
# by dropping NULL model
g <- ggplot(train_metaGRS_perf) +
  aes(y=model_fit_estimate, ymin=model_fit_L95, ymax=model_fit_U95,
      x=pgs_case_def, color=lambda, shape=prefilter) +
  geom_errorbar(width=0, alpha=0.5, position=position_dodge(width=0.9)) +
  geom_point(position=position_dodge(width=0.9)) +
  scale_colour_manual(name="elasticnet best model", 
                      values=c("lambda.min"="#762a83", "lambda.1se"="#1b7837")) +
  scale_shape_manual(name="GRS prefiltering",
                     values=c("none"=19, "sig_auc"=18, "no_auc_ci_overlap"=15)) +
  facet_wrap(~ model_type + model_case_def, scales="free_y", ncol=3) +
  xlab("") + 
  ylab("AUC or C-index (incident) (+/- 95% CI) for T2D") +
	theme_bw() +
	theme(
		legend.position="bottom",
		panel.grid.major.x=element_blank(),
		panel.grid.minor.x=element_blank(),
		axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=0.5),
		axis.text.y=element_text(size=8),
		axis.title=element_text(size=10),
		strip.text.x=element_text(size=8),
		strip.text.y=element_text(size=8),
    legend.title=element_text(size=8),
    legend.text=element_text(size=8)
	)
ggsave(g, width=30, height=30, units="in", file="output/metaGRS/hyperparam_selection/control_def_selection_trainingset.png")

g <- ggplot(test_metaGRS_perf) +
  aes(y=model_fit_estimate, ymin=model_fit_L95, ymax=model_fit_U95,
      x=pgs_case_def, color=lambda, shape=prefilter) +
  geom_errorbar(width=0, alpha=0.5, position=position_dodge(width=0.9)) +
  geom_point(position=position_dodge(width=0.9)) +
  scale_colour_manual(name="elasticnet best model", 
                      values=c("lambda.min"="#762a83", "lambda.1se"="#1b7837")) +
  scale_shape_manual(name="GRS prefiltering",
                     values=c("none"=19, "sig_auc"=18, "no_auc_ci_overlap"=15)) +
  facet_wrap(~ model_type + model_case_def, scales="free_y", ncol=3) +
  xlab("") + 
  ylab("AUC or C-index (incident) (+/- 95% CI) for T2D") +
	theme_bw() +
	theme(
		legend.position="bottom",
		panel.grid.major.x=element_blank(),
		panel.grid.minor.x=element_blank(),
		axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=0.5),
		axis.text.y=element_text(size=8),
		axis.title=element_text(size=10),
		strip.text.x=element_text(size=8),
		strip.text.y=element_text(size=8),
    legend.title=element_text(size=8),
    legend.text=element_text(size=8)
	)
ggsave(g, width=30, height=30, units="in", file="output/metaGRS/hyperparam_selection/control_def_selection_testset.png")

# Its hard to see any differences in performance based on control definition in trainingset, so we can skip deciding that for now.
# In test set, we see lambda.min models consistently outperform lambda.1se models, indicating the models with best AUC in 
# cross validation are not overfit compared to the simpler models within 1 standard error (consistent with other experiments we've
# run in UKB), so we can take the best 'lambda.min' model in each instance
cvfit <- cvfit[lambda == "lambda.min"]
train_metaGRS_perf <- train_metaGRS_perf[lambda == "lambda.min"]
test_metaGRS_perf <- test_metaGRS_perf[lambda == "lambda.min"]

# We can also see that in both the training and test datasets, prefiltering the input list of GRSs before running elasticnet
# doesn't improve performance, so we can conclude that including all the GRSs doesn't add noise and lead to overfitting
cvfit <- cvfit[prefilter == "none"]
train_metaGRS_perf <- train_metaGRS_perf[prefilter == "none"]
test_metaGRS_perf <- test_metaGRS_perf[prefilter == "none"]

# Lets take another look at control status
train_metaGRS_perf <- train_metaGRS_perf[order(-model_fit_estimate)][order(model_case_def)][order(model_type)]
train_metaGRS_perf[, rank := 1:.N, by=.(model_type, model_case_def)]
g <- ggplot(train_metaGRS_perf) +
  aes(y=model_fit_estimate, ymin=model_fit_L95, ymax=model_fit_U95,
      x=pgs_case_def, color=factor(rank)) +
  geom_errorbar(width=0, alpha=0.5, position=position_dodge(width=0.9)) +
  geom_point(shape=19, position=position_dodge(width=0.9)) +
  scale_colour_manual(name="Rank", values=structure(brewer.pal(name="Spectral", n=10), names=1:10)) + 
  guides(colour = guide_legend(nrow = 1)) +
  facet_wrap(~ model_type + model_case_def, scales="free_y", ncol=10) +
  xlab("") + 
  ylab("AUC or C-index (incident) (+/- 95% CI) for T2D") +
	theme_bw() +
	theme(
		legend.position="bottom",
		panel.grid.major.x=element_blank(),
		panel.grid.minor.x=element_blank(),
		axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=0.5),
		axis.text.y=element_text(size=8),
		axis.title=element_text(size=10),
		strip.text.x=element_text(size=8),
		strip.text.y=element_text(size=8),
    legend.title=element_text(size=8),
    legend.text=element_text(size=8)
	)
ggsave(g, width=30, height=10, units="in", file="output/metaGRS/hyperparam_selection/control_def_selection_filtered_trainingset.png")

test_metaGRS_perf <- test_metaGRS_perf[order(-model_fit_estimate)][order(model_case_def)][order(model_type)]
test_metaGRS_perf[, rank := 1:.N, by=.(model_type, model_case_def)]
g <- ggplot(test_metaGRS_perf) +
  aes(y=model_fit_estimate, ymin=model_fit_L95, ymax=model_fit_U95,
      x=pgs_case_def, color=factor(rank)) +
  geom_errorbar(width=0, alpha=0.5, position=position_dodge(width=0.9)) +
  geom_point(shape=19, position=position_dodge(width=0.9)) +
  scale_colour_manual(name="Rank", values=structure(brewer.pal(name="Spectral", n=10), names=1:10)) + 
  guides(colour = guide_legend(nrow = 1)) +
  facet_wrap(~ model_type + model_case_def, scales="free_y", ncol=10) +
  xlab("") + 
  ylab("AUC or C-index (incident) (+/- 95% CI) for T2D") +
	theme_bw() +
	theme(
		legend.position="bottom",
		panel.grid.major.x=element_blank(),
		panel.grid.minor.x=element_blank(),
		axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=0.5),
		axis.text.y=element_text(size=8),
		axis.title=element_text(size=10),
		strip.text.x=element_text(size=8),
		strip.text.y=element_text(size=8),
    legend.title=element_text(size=8),
    legend.text=element_text(size=8)
	)
ggsave(g, width=30, height=10, units="in", file="output/metaGRS/hyperparam_selection/control_def_selection_filtered_testset.png")

# Some observations:
#
# (1) Within prevalent or incident disease analyses, rank order is identical within training data no matter how you define control
#     status for the test analyses
# (2) Rank differs between prevalent and incident disease analysis
# (3) Differences in rank are not the same when looking in test dataset
# (4) Differences in performance are negligible in any case
#
# Therefore, best to default to simplest definition; i.e. no filtering of controls
cvfit <- cvfit[case_def == "basic_adjudicated"]
train_metaGRS_perf <- train_metaGRS_perf[pgs_case_def == "basic_adjudicated"]
test_metaGRS_perf <- test_metaGRS_perf[pgs_case_def == "basic_adjudicated"]

# Extract final metaGRS
varweights <- fread("output/metaGRS/train/combined_basic_adjudicated/prefilter_none/candidate_metaGRS.txt.gz")
varweights <- varweights[,.(rsid, chr, pos, effect_allele, other_allele, weight=alpha.1_lambda.min)]
varweights <- varweights[weight != 0] # drop 4,173 SNPs not present in contributing GWAS
fwrite(varweights, sep="\t", quote=FALSE, compress="gzip", file="output/metaGRS/hyperparam_selection/T2D_metaGRS.txt.gz")

metaGRS <- fread("output/metaGRS/all_candidate_metaGRS_lvls/collated_candidate_metaGRSs_prevalent_plus_incident_T2D/collated_scores.sscore.gz",
                 select=c("IID", "combined_basic_adjudicated_none_alpha_1_lambda_min"))
setnames(metaGRS, c("eid", "T2D_metaGRS"))
fwrite(metaGRS, sep="\t", quote=FALSE, file="data/UKB/T2D_metaGRS_levels.txt")

cvcoef <- fread("output/metaGRS/train/combined_basic_adjudicated/prefilter_none/enfit_best_coefs.txt")
cvcoef <- cvcoef[alpha == 1 & best.model == "lambda.min"]
gwass <- list.dirs("output/ldpred2/train/", recursive=FALSE, full.names=FALSE)
cvcoef <- cvcoef[coefficient %in% gwass]
cvcoef <- cvcoef[order(-abs(beta))]
cvcoef <- cvcoef[, .(GRS=coefficient, metaGRS_weight=beta)]
fwrite(cvcoef, sep="\t", quote=FALSE, file="output/metaGRS/hyperparam_selection/T2D_metaGRS_component_GRS_contributions.txt")
