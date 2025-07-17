library(data.table)
library(foreach)
library(caret)
library(glmnet)
library(ggplot2)

# Need to run script once per ancestry and GWAS
ancestry <- "EUR"

# Create output directory
outdir <- sprintf("output/metaPRS/train/%s/", ancestry)
system(sprintf("mkdir -p %s", outdir), wait=TRUE)

# Load in list of possible GWASs
gwas_list <- fread('data/gwas_summary_stats/curated_gwas_table.txt')

# Load phenotype data for this ancestry and filter to the metaPRS training samples 
pheno <- fread(sprintf("data/All_of_Us/%s/phenotypes.txt", ancestry))
pheno <- pheno[(metaPRS_train_samples)]

# Load in set of ancestry-specific LDpred2 trained PRS
input_prs <- fread(sprintf("output/ldpred2/PRS/%s/optimal_prs_levels.txt", ancestry))

# Filter to people in the metaPRS training samples 
input_prs <- input_prs[IID %in% pheno$person_id] 

# Restrict phenotype data to people with genetics
pheno <- pheno[person_id %in% input_prs$IID]

# Adjust each PRS for 20 PCs
pcs <- pheno[,.(person_id, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, 
                PC11, PC12, PC13, PC14, PC15, PC16, PC17, PC18, PC19, PC20)]

join <- function(l, r) {  l[r, on = .(IID)]  }
adj_prs <- foreach(this_gwas = gwas_list$PRS, .combine=join) %do% {
  this_prs <- input_prs[, c("IID", this_gwas), with=FALSE]
  setnames(this_prs, this_gwas, "PRS")
  this_prs <- this_prs[pcs, on = .(IID=person_id)]
  this_prs[, PRS := scale(lm(scale(PRS) ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 +
    PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 +
    PC17 + PC18 + PC19 + PC20)$residuals)]
  this_prs <- this_prs[, .(IID, PRS)]
  setnames(this_prs, "PRS", this_gwas)
  return(this_prs)
}

# Split into 10 folds for elasticnet cross-validation balancing T2D case status and sex
pheno[, foldid := createFolds(paste(T2D, sex), k=10, list=FALSE)]

# Write out fold assignments in case we need them later
fwrite(pheno[,.(person_id, foldid)], sep="\t", quote=FALSE, file=sprintf("%s/cross_validation_fold_assignments.txt", outdir))

# Code model covariates
pheno[, sex := factor(sex, levels=c("Female", "Male"))]
pheno[, age := scale(age)]

# Combine PRSs with covariates, T2D case status, and cross-validation fold assignment
pheno <- pheno[input_prs, on = .(person_id=IID)]

# Extract input matrix for model fitting
xmat <- model.matrix(~ 0 + ., pheno[, c(gwas_list$PRS, "age", "sex"), with=FALSE])

# Setup list of alpha mixing parameters to search across (controls balance of ridge vs. lasso)
alphas <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1) # 0 = ridge, 1 = lasso

# Run elasticnet
enfits <- foreach(alpha=alphas) %do% {
	cv.glmnet(y = pheno$T2D, x=xmat, foldid=pheno$foldid, type.measure="auc", family="binomial", alpha=alpha, trace.it=TRUE)
}
names(enfits) <- alphas
saveRDS(enfits, file=sprintf("%s/enfit.rds", out_dir))

# build table of AUCs 
fitdt <- foreach(idx = seq_along(alphas), .combine=rbind) %do% {
	data.table(alpha=alphas[idx], lambda=enfits[[idx]][["lambda"]],
						 mean_AUC=enfits[[idx]][["cvm"]], AUC_sd=enfits[[idx]][["cvsd"]],
						 mean_AUC_minus_sd=enfits[[idx]][["cvlo"]], mean_AUC_plus_sd=enfits[[idx]][["cvup"]],
						 nonzero=enfits[[idx]][["nzero"]],
						 lambda.min=enfits[[idx]][["lambda.min"]],
						 lambda.1se=enfits[[idx]][["lambda.1se"]])
}
fwrite(fitdt, sep="\t", quote=FALSE, file=sprintf("%s/enfit_model_fit.txt", out_dir))

# Build table containing information about the best fit
bestfit <- fitdt[lambda == lambda.min] # best fit for each alpha
bestfit <- bestfit[,.SD[which.max(mean_AUC)]] # overall best fit
fwrite(bestfit, sep="\t", quote=FALSE, file=sprintf("%s/enfit_best_model_fit.txt", out_dir))

# Plot model fits
g <- ggplot(fitdt) +
	aes(x=log(lambda), y=mean_AUC, ymin=mean_AUC_minus_sd, ymax=mean_AUC_plus_sd,
			fill=factor(alpha), colour=factor(alpha)) +
	geom_ribbon(colour="#00000000", alpha=0.3, show.legend=FALSE) +
	geom_line() +
	geom_vline(data=bestfit, aes(xintercept=log(lambda), colour=factor(alpha)), linetype="dashed", show.legend=FALSE) +
	scale_fill_manual(name="alpha", values=c("0"="#5e4fa2", "0.1"="#3288bd", "0.25"="#66c2a5",
										"0.5"="#ffff33", "0.75"="#f46d43", "0.9"="#d53e4f", "1"="#9e0142")) +
	scale_colour_manual(name="alpha", values=c("0"="#5e4fa2", "0.1"="#3288bd", "0.25"="#66c2a5",
										"0.5"="#ffff33", "0.75"="#f46d43", "0.9"="#d53e4f", "1"="#9e0142")) +
	xlab("Log(lambda)") +
  ylab("AUC (+/- SD)") +
	theme_bw()
ggsave(g, width=10, height=6, units="in", file=sprintf("%s/enfit.png", out_dir))

# Collate the coefficients for the best model (log odds for each PRS, age, and sex)
encoef <- as.data.table(as.matrix(coef(enfits[[which(alphas == bestfit$alpha)]], s="lambda.min")), keep.rownames="coefficient")
setnames(encoef, "1", "logOR")
fwrite(encoef, sep="\t", quote=FALSE, file=sprintf("%s/best_enfit_coefs.txt", out_dir))

# Load variant weights so we can derive the PRS weight file for the ancestry-specific metaPRS
var_weights <- fread(sprintf("output/ldpred2/PRS/%s/optimal_prs_var_weights.txt", ancestry))

# Remove PRS that were not selected by elasticnet
var_weights <- var_weights[PRS %in% encoef$coefficient]

# Multiply out LDpred2 weights by the log odds from elasticnet
var_weights[encoef, on = .(PRS=coefficient), weight := weight * logOR]

# Sum over all contributing PRS
metaPRS_weights <- var_weights[, .(weight=sum(weight)), by=.(AoU_varID, chr, pos, effect_allele, other_allele)]

# Write out ancestry-specific metaPRS
fwrite(metaPRS_weights, sep="\t", quote=FALSE, compress="gzip", file=sprintf("%s/T2D_%s_metaPRS.txt.gz", outdir, ancestry)

