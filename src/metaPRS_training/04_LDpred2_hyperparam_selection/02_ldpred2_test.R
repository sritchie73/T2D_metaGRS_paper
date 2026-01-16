library(data.table)
library(foreach)
library(ggplot2)
library(ggthemes)
library(cowplot)
source("src/functions/glm_test.R")

# Need to run script once per ancestry and GWAS
ancestry <- "EUR"
gwas <- "T2D_MVP_Multi"

# Create output directory
outdir <- sprintf("output/ldpred2/hyperparam_selection/%s/%s", ancestry, gwas)
system(sprintf("mkdir -p %s", outdir), wait=TRUE)

# Load phenotype data for this ancestry and filter to the metaPRS training samples 
pheno <- fread(sprintf("data/All_of_Us/%s/phenotypes.txt", ancestry))
pheno <- pheno[(metaPRS_train_samples)]

# Add in computed PRS levels for all possible LDpred2 hyperparameters
candidate_prs <- fread(sprintf("output/ldpred2/all_hyperparam_grs_lvls/%s/%s/collated_scores.sscore.gz", ancestry, gwas))

# Filter to metaPRS training cohort samples
candidate_prs <- candidate_prs[IID %in% pheno$person_id]

# Restrict phenotype data to people with genetics
pheno <- pheno[person_id %in% candidate_prs$IID]

# Test each candidate PRS for association with T2D
prs_assocs <- foreach(this_prs = unique(candidate_prs$score), .combine=rbind) %do% { 
	# Add specific PRS to phenotype data for testing
	pheno[candidate_prs[score == this_prs], on = .(person_id=IID), PRS := i.value]

	# Adjust for 20 PCs and standardise
	pheno[, PRS := scale(lm(scale(PRS) ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 +
		PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + 
		PC17 + PC18 + PC19 + PC20)$residuals)]

	# Test PRS for association with T2D in logistic regression
	suppressMessages(g1 <- glm.test(formula=T2D ~ PRS + age + sex, event_col="T2D", data=pheno, ci.method="wald"))
	cbind(LDpred2_model=this_prs, g1)
} 
fwrite(prs_assocs, sep="\t", quote=FALSE, file=sprintf("%s/all_model_performance.txt", outdir))

# Generate diagnostic plots
cat("Generating diagnositic plots...\n")

if (any(unique(candidate_prs$score) %like% "grid")) { # grid model missing where h2_est < 0
	# Load LDpred2 grid model parameters
	traindir <- sprintf("output/ldpred2/train/%s/%s", ancestry, gwas)
	params <- fread(sprintf("%s/grid_model_parameters.txt", traindir))
	params[, paramset := .I]

	# Plot performance against grid parameters
	grid_perf <- prs_assocs[LDpred2_model %like% "grid_" & coefficient == "PRS"]
	grid_perf[, paramset := as.integer(gsub("grid_", "", LDpred2_model))]
	grid_perf <- merge(prs_assocs, params, by="paramset", all.x=TRUE)

	g1 <- ggplot(grid_perf) + 
		aes(x=grid_param_p, y=AUC, ymin=AUC.L95, ymax=AUC.U95, color=as.factor(grid_param_h2)) +
		theme_bigstatsr() +
		geom_errorbar(width=0, alpha=0.5) +
		geom_point() +
		geom_line() + 
		scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
		facet_wrap(~ grid_param_sparse, labeller = label_both) +
		labs(y = "AUC (95% CI)", color = "h2") +
		theme_bw() +
		theme(legend.position = "top", panel.spacing = unit(1, "lines"))

	g2 <- ggplot(grid_perf) +
		aes(x=grid_param_p, y=OR, ymin=OR.L95, ymax=OR.U95, color=as.factor(grid_param_h2)) +
		theme_bigstatsr() +
		geom_hline(yintercept=1, linetype=2) +
		geom_errorbar(width=0, alpha=0.5) +
		geom_point() +
		geom_line() +
		scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
		facet_wrap(~ grid_param_sparse, labeller = label_both) +
		labs(y = "OR (95% CI)", color = "h2") +
		theme_bw() +
		theme(legend.position = "top", panel.spacing = unit(1, "lines"))

	g <- plot_grid(g1, g2, nrow = 2)
	ggsave(g, width=7.2, height=7.2, units="in", file=sprintf("%s/grid_parameter_performance.png", outdir))
}

# Get lasso parameters
params <- fread(sprintf("output/ldpred2/train/%s/lassosum2_model_parameters.txt", traindir))

# Plot performance against lasso parameters
lassosum2_perf <- prs_assocs[LDpred2_model %like% "lassosum_" & coefficient == "PRS"]
lassosum2_perf[, paramset := as.integer(gsub("lassosum_", "", LDpred2_model))]
lassosum2_perf <- merge(prs_assocs, params, by="paramset", all.x=TRUE)

g1 <- ggplot(lassosum2_perf) + 
	aes(x=lambda, y=AUC, ymin=AUC.L95, ymax=AUC.U95, color=as.factor(delta)) +
	theme_bigstatsr() +
	geom_errorbar(width=0, alpha=0.5) +
	geom_point() +
	geom_line() +
	scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
	labs(y = "AUC (95% CI)", color = "delta") +
	theme_bw() +
	theme(legend.position = "top", panel.spacing = unit(1, "lines"))

g2 <- ggplot(lassosum2_perf) +
	aes(x=lambda, y=OR, ymin=OR.L95, ymax=OR.U95, color=as.factor(delta)) +
	theme_bigstatsr() +
	geom_hline(yintercept=1, linetype=2) +
	geom_errorbar(width=0, alpha=0.5) +
	geom_point() +
	geom_line() +
	scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
	labs(y = "OR (95% CI)", color = "delta") +
	theme_bw() +
	theme(legend.position = "top", panel.spacing = unit(1, "lines"))

g <- plot_grid(g1, g2, nrow = 2)
ggsave(g, width=7.2, height=7.2, units="in", file=sprintf("%s/lassosum2_parameter_performance.png", outdir))

# Simplify the grid and lasso models down to the optimal parameters, then compare the four LDpred2 models (inf, grid, auto, lassosum)
model_comp <- prs_assocs[coefficient == "PRS"]
model_comp[, LDpred2_model := gsub("_.*", "", LDpred2_model)]
model_comp <- model_comp[, .SD[which.max(AUC)], by=LDpred2_model]

g1 <- ggplot(model_comp) +
  aes(x=AUC, xmin=AUC.L95, xmax=AUC.U95, y=LDpred2_model) +
  geom_errorbarh(height=0, alpha=0.5) +
  geom_point() +
  geom_line() +
  xlab("AUC (95% CI)") + 
  theme_bw()

g2 <- ggplot(model_comp) +
  aes(x=OR, xmin=OR.L95, xmax=OR.U95, y=LDpred2_model) +
  geom_vline(xintercept=1, linetype=2) +
  geom_errorbarh(height=0, alpha=0.5) +
  geom_point() +
  geom_line() +
  xlab("OR (95% CI)") +
  theme_bw()

g <- plot_grid(g1, g2, nrow = 1)
ggsave(g, width=7.2, height=3.5, units="in", file=sprintf("%s/ldpred2_model_comparison.png", outdir))
 
