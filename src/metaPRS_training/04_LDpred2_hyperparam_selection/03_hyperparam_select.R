library(data.table)
library(foreach)

# Need to run script once per ancestry and GWAS
ancestry <- "EUR"

# Load in list of possible GWASs
gwas_list <- fread('data/gwas_summary_stats/curated_gwas_table.txt')

# Create output directory
outdir <- sprintf("output/ldpred2/PRS/%s/", ancestry)
system(sprintf("mkdir -p %s", outdir), wait=TRUE)

# Loop through each GWAS and find the optimal LDpred2 hyperparameters
optimal_models <- foreach(this_gwas = gwas_list$PRS, .combine=rbind) %do% {
  # Load in association summary statistics
  ldpred_assocs <- fread(sprintf("output/ldpred2/hyperparam_selection/%s/%s/all_model_performance.txt", ancestry, this_gwas))

  # Determine the optimal model
	ldpred_assocs <- ldpred_assocs[coefficient == "PRS"]
	ldpred_assocs <- ldpred_assocs[, .SD[which.max(AUC)]]

  # Extract relevant information about the model
	ldpred_assocs[, LDpred2_model_group := gsub("_.*", "", LDpred2_model)]
  ldpred_assocs[(LDpred2_model_group %like% "_"), paramset := as.integer(gsub(".*_", "", LDpred2_model))]

  if (ldpred_assocs$LDpred2_model_group == "infinitesimal") {
    h2_est <- fread(sprintf("output/ldpred2/train/%s/%s/infinitesimal_model_parameters.txt", ancestry, this_gwas)) 
    ldpred_assocs[, parameters := sprintf("LD score regression heritability (h2): %.3f", h2_est$h2)]
  } else if (ldpred_assocs$LDpred2_model_group == "grid") {
    grid_params <- fread(sprintf("output/ldpred2/train/%s/%s/grid_model_parameters.txt", ancestry, this_gwas))
    grid_params <- grid_params[ldpred_assocs$paramset]
    ldpred_assocs[, parameters := sprintf("Proportion of causal variants (p): %.3f, heritability captured by variants (h2): %.3f, Sparse model: %s", grid_params$p, grid_params$h2, grid_params$sparse)]
    ldpred_assocs[, LDpred2_model_group := "grid-search"] 
  } else if (ldpred_assocs$LDpred2_model_group == "auto") {
    auto_params <- fread(sprintf("output/ldpred2/train/%s/%s/auto_model_parameters.txt", ancestry, this_gwas))
    auto_chain_qc <- fread(sprintf("output/ldpred2/hyperparam_selection/%s/%s/ldpred2_auto_chain_qc.txt", ancestry, this_gwas))
    auto_params <- auto_params[(auto_chain_qc$keep), .(p_est=mean(p_est), h2_est=mean(h2_est))] 
    ldpred_assocs[, parameters := sprintf("Proportion of causal variants (p): %.3f, heritability captured by variants (h2): %.3f", auto_params$p_est, auto_params$h2_est)]
    ldpred_assocs[, LDpred2_model_group := "automatic"]
  } else if (ldpred_assocs$LDpred2_model_group == "lassosum2") {
    lasso_params <- fread(sprintf("output/ldpred2/train/%s/%s/lassosum2_model_parameters.txt", ancestry, this_gwas))
    lasso_params <- lasso_params[ldpred_assocs$paramset]
    ldpred_assocs[, parameters := sprintf("Delta: %.3f, Lambda: %.3f", lasso_params$delta, lasso_params$lambda)]
    ldpred_assocs[, LDpred2_model_group := "lasso-sum"]
  }
  
  # Return
  return(cbind(PRS=this_gwas, ldpred_assocs))  
}

# Write out information about these models and their association with T2D
fwrite(optimal_models, sep="\t", quote=FALSE, file=sprintf("%s/optimal_prs_info.txt", outdir))

# Extract and collate the set of PRS levels for metaPRS training
join <- function(l, r) {  l[r, on = .(IID)]  }
prs_levels <- foreach(this_gwas = gwas_list$PRS, .combine=join) %do% {
  ldpred2_model <- optimal_models[PRS == this_gwas, LDpred2_model]
  if (ldpred2_model == "auto") {
    this_prs <- fread(sprintf("output/ldpred2/hyperparam_selection/%s/%s/ldpred2_auto_model_prs.txt.gz", ancestry, this_gwas))
    this_prs[, score := NULL]
    setnames(this_prs, "value", this_gwas)
  } else {
    this_prs <- fread(sprintf("output/ldpred2/all_hyperparam_grs_lvls/%s/%s/collated_scores.sscore.gz", ancestry, this_gwas))
    this_prs <- this_prs[,c("IID", ldpred2_model),with=FALSE]
    setnames(this_prs, ldpred2_model, this_gwas)
  }
  return(this_prs)
}
fwrite(prs_levels, sep="\t", quote=FALSE, file=sprintf("%s/optimal_prs_levels.txt", outdir))

# Extract and collate the set of variant weights for the optimal PRSs, so we can later combine them into the 
# ancestry-specific metaPRSs
var_weights <- foreach(this_gwas = gwas_list$PRS, .combine=rbind) %do% {
  ldpred2_model <- optimal_models[PRS == this_gwas, LDpred2_model]
  if (ldpred2_model == "auto") {
    this_weights <- fread(sprintf("output/ldpred2/hyperparam_selection/%s/%s/ldpred2_auto_model_varweights.txt.gz", ancestry, this_gwas))
    setnames(this_weights, "auto", "weight")
  } else {
    this_weights <- fread(sprintf("output/ldpred2/train/%s/%s/ldpred2_pgs_varweights.txt.gz", ancestry, this_gwas), 
      select=c("AoU_varID", "chr", "pos", "effect_allele", "other_allele", ldpred2_model))
    sentames(this_weights, ldpred2_model, "weight")
  }
  cbind(PRS=this_gwas, this_weights)
}
fwrite(var_weights, sep="\t", quote=FALSE, file=sprintf("%s/optimal_prs_var_weights.txt", outdir)) 

