library(data.table)
library(foreach)
library(ggplot2)

# Set up parallelisation (e.g. for fread/fwrite)
source("src/functions/par_setup.R")

# Determine which GWAS(s) we're working with
gwass <- list.dirs("data/filtered_sumstats", recursive=FALSE, full.names=FALSE)

# Create output directories
for (gwas in gwass) {
  system(sprintf("mkdir -p output/ldpred2/hyperparam_selection/%s", gwas), wait=TRUE)
}

# Load in all model performance
perf <- lapply(gwass, function(gwas) {
  fread(sprintf("output/ldpred2/test/%s/all_model_performance.txt", gwas), tmpdir="tmp")
})
names(perf) <- gwass
perf <- rbindlist(perf, idcol="gwas")

# Keep only one set of null models
null_models <- perf[pgs == "none"]
null_models[, gwas := NULL]
null_models <- unique(null_models)
perf <- perf[pgs != "none"]

# Drop auto chains, keep only final model
perf <- perf[pgs == "pgs_auto_final2" | !(pgs %like% "_auto_")]

# Flag LDpred2 model type
perf[, ldpred2_model := fcase(
  pgs == "pgs_inf", "infinitesimal",
  pgs == "pgs_auto_final2", "automatic",
  pgs %like% "_grid_", "grid-search",
  pgs %like% "_lassosum2", "lasso-sum"
)]

# Obtain best model by AUC (or C.index) for each LDpred2 model
best <- perf[metric %in% c("AUC", "C.index"), .SD[which.max(estimate)], by=.(gwas, ldpred2_model, case_def)]

# Flag instances where AUC (or C-index) overlaps with null:
best[null_models, on = .(case_def), pval := fcase(
  L95 > i.U95, 0.001, # 95% CI for model with PGS does not overlap with null model 95% CI
  L95 > i.estimate, 0.025, # 95% CI for model with PGS does not overal with null model point estimate (i.e. P < 0.05) but *does* overlap with null model 95% CI
  default = 0.5 # 95% CI for model with PGS overlaps null model point estimate
)]

# Get their corresponding OR / HR
toadd <- perf[metric %in% c("OR", "HR")][best[,.(gwas, case_def, pgs)], on = .(gwas, case_def, pgs)]
best <- rbind(best, toadd)

# Flag overall best model
best[, best_model := FALSE]
best_best <- best[metric %in% c("AUC", "C.index"), .SD[which.max(estimate)], by=.(gwas, case_def)]
best[best_best, on = .(gwas, case_def, pgs), best_model := TRUE]

# Write out best model(s)
fwrite(best, sep="\t", quote=FALSE, file="output/ldpred2/hyperparam_selection/best_hyperparams.txt")
fwrite(best[(best_model)], sep="\t", quote=FALSE, file="output/ldpred2/hyperparam_selection/best_models.txt")
fwrite(null_models, sep="\t", quote=FALSE, file="output/ldpred2/hyperparam_selection/null_models.txt")

# For each GWAS, plot performance of best models for each LDpred2 model
best[, case_group := gsub("^[a-z]*_", "", case_def)] # for plot facets
null_models[, case_group := gsub("^[a-z]*_", "", case_def)] # for plot facets  
best[, ldpred2_model := factor(ldpred2_model, levels=c("infinitesimal", "grid-search", "automatic", "lasso-sum"))]
best[, type := factor(type, levels=c("prevalent", "incident", "prevalent + incident"))]
for (this_gwas in gwass) {
  this_best <- best[gwas == this_gwas]

  # Plot AUC (or C-index)
  g <- ggplot(this_best[metric %in% c("AUC", "C.index")]) +
    aes(x=ldpred2_model, y=estimate, ymin=L95, ymax=U95, color=best_model) +
    geom_rect(inherit.aes=FALSE, data=null_models, aes(ymin=L95, ymax=U95), xmin=-Inf, xmax=Inf, fill="#bdbdbd", alpha=0.3) +
    geom_hline(data=null_models, aes(yintercept=estimate), linetype=2) +
    geom_errorbar(width=0, alpha=0.5) +
    geom_point(shape=19) + 
    scale_colour_manual(values=c("TRUE"="#810f7c", "FALSE"="#8c6bb1")) +
    facet_grid(type ~ case_group) +
		ylab("AUC or C-index (incident) (+/- 95% CI) for T2D") +
    xlab("LDpred2 method") +
		theme_bw() +
		theme(
      legend.position = "none", 
      panel.grid.major.x=element_blank(),
      panel.grid.minor.x=element_blank(),
      axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8),
      axis.text.y=element_text(size=8),
      axis.title=element_text(size=10), 
      strip.text.x=element_text(size=8)
    )
  ggsave(g, width=16, height=7.2, units="in", file=sprintf("output/ldpred2/hyperparam_selection/%s/AUC_Cindex_all_case_def.png", this_gwas))

  # Plot OR (or HR)
  g <- ggplot(this_best[metric %in% c("OR", "HR")]) +
    aes(x=ldpred2_model, y=estimate, ymin=L95, ymax=U95, color=best_model) +
    geom_hline(yintercept=1, linetype=2) +
    geom_errorbar(width=0, alpha=0.5) +
    geom_point(shape=19) + 
    scale_colour_manual(values=c("TRUE"="#810f7c", "FALSE"="#8c6bb1")) +
    facet_grid(type ~ case_group) +
		ylab("OR or HR (incident) (+/- 95% CI) for T2D") +
    xlab("LDpred2 method") +
		theme_bw() +
		theme(
      legend.position = "none", 
      panel.grid.major.x=element_blank(),
      panel.grid.minor.x=element_blank(),
      axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8),
      axis.text.y=element_text(size=8),
      axis.title=element_text(size=10), 
      strip.text.x=element_text(size=8)
    )
  ggsave(g, width=16, height=7.2, units="in", file=sprintf("output/ldpred2/hyperparam_selection/%s/OR_HR_all_case_def.png", this_gwas))
}

# Generate plots comparing all GWASs for each case-control definition set
for (cdef in unique(best$case_def)) {
  this_best <- best[case_def == cdef]
  this_null <- null_models[case_def == cdef]

  # Plot AUC (or C-index)
  g <- ggplot(this_best[metric %in% c("AUC", "C.index")]) +
    aes(x=ldpred2_model, y=estimate, ymin=L95, ymax=U95, color=best_model) +
    geom_rect(inherit.aes=FALSE, ymin=this_null$L95, ymax=this_null$U95, xmin=-Inf, xmax=Inf, fill="#bdbdbd", alpha=0.3) +
    geom_hline(yintercept=this_null$estimate, linetype=2) +
    geom_errorbar(width=0, alpha=0.5) +
    geom_point(shape=19) +
    scale_colour_manual(values=c("TRUE"="#810f7c", "FALSE"="#8c6bb1")) +
    facet_wrap(~ gwas, ncol=10, scales="free_y") +
    ylab(sprintf("%s (+/- 95%% CI) for %s T2D", this_null$metric, this_null$type)) +
    xlab("LDpred2 method") +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.major.x=element_blank(),
      panel.grid.minor.x=element_blank(),
      axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8),
      axis.text.y=element_text(size=8),
      axis.title=element_text(size=10),
      strip.text.x=element_text(size=8)
    )
  ggsave(g, width=16, height=7.2, units="in", file=sprintf("output/ldpred2/hyperparam_selection/%s_%s_compare.png", cdef, this_null$metric))

  # Plot OR (or HR)
  g <- ggplot(this_best[metric %in% c("OR", "HR")]) +
    aes(x=ldpred2_model, y=estimate, ymin=L95, ymax=U95, color=best_model) +
    geom_hline(yintercept=1, linetype=2) +
    geom_errorbar(width=0, alpha=0.5) +
    geom_point(shape=19) +
    scale_colour_manual(values=c("TRUE"="#810f7c", "FALSE"="#8c6bb1")) +
    facet_wrap(~ gwas, ncol=10, scales="free_y") +
    ylab(sprintf("%s (+/- 95%% CI) for %s T2D", this_best[metric %in% c("OR", "HR"), metric][1], this_null$type)) +
    xlab("LDpred2 method") +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.major.x=element_blank(),
      panel.grid.minor.x=element_blank(),
      axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8),
      axis.text.y=element_text(size=8),
      axis.title=element_text(size=10),
      strip.text.x=element_text(size=8)
    )
  ggsave(g, width=16, height=7.2, units="in", file=sprintf("output/ldpred2/hyperparam_selection/%s_%s_compare.png", cdef, this_best[metric %in% c("OR", "HR"), metric][1]))
}

