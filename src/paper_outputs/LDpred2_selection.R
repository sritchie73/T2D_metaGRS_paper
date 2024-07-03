library(data.table)
library(foreach)
library(ggplot2)
library(cowplot)

system("mkdir -p output/paper_outputs/{figures,tables}", wait=TRUE)

#####################################################################
# Show example of selection of LDpred2 models for the T2D Exome GWAS
#####################################################################

# Load LDpred2 per-model AUCs
ldpred_perf <- fread("output/ldpred2/hyperparam_selection/best_hyperparams.txt")
ldpred_perf <- ldpred_perf[case_def == "combined_basic_adjudicated" & gwas == "T2D_Exome" & metric == "AUC"]

# Load performance for all parameters tested
all_perf <- fread("output/ldpred2/test/T2D_Exome/all_model_performance.txt")
all_perf <- all_perf[type == "prevalent + incident" & metric == "AUC" & case_def == "combined_basic_adjudicated"]

# Curate grid model performance
ldsc_results <- readRDS("output/ldpred2/train/T2D_Exome/ldsc_results.rds")
h2_est <- ldsc_results["h2"]
h2_seq <- round(h2_est * c(0.3, 0.7, 1, 1.4), 4)
p_seq <- signif(bigsnpr::seq_log(1e-5, 1, length.out = 21), 2)
params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))
params <- as.data.table(params, keep.rownames="paramset")

grid_perf <- all_perf[pgs %like% "_grid_"]
grid_perf[, paramset := gsub(".*_", "", pgs)]
grid_perf <- grid_perf[params, on = .(paramset), nomatch=0]

# Curate lassosum model performance
params <- fread("output/ldpred2/train/T2D_Exome/ldpred2_lassosum2_parameters.txt")
lasso_perf <- all_perf[pgs %like% "_lassosum2_"]
lasso_perf[, paramset := as.integer(gsub(".*_", "", pgs))]
lasso_perf <- lasso_perf[params, on = .(paramset), nomatch=0]

# Create plot showing selection of LDpred2 model
g1 <- ggplot(ldpred_perf) + 
  aes(x=ldpred2_model, y=estimate, ymin=L95, ymax=U95) +
  geom_errorbar(width=0) +
  geom_point(shape=21, fill="white") +
  scale_y_continuous("AUC (95% CI)", limits=range(c(all_perf$L95, all_perf$U95))) +
  xlab("LDpred2 model") +
  ggtitle("UK Biobank", subtitle="T2D prevalence") + 
  theme_bw() +
  theme(
    axis.text.x=element_text(size=7, angle=90, hjust=1, vjust=0.5, color="black"),
    axis.text.y=element_text(size=6), axis.title=element_text(size=8),
    panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank()
  )

# Create plot showing grid-search parameter tuning
g2 <- ggplot(grid_perf) +
  aes(x=p, y=estimate, ymin=L95, ymax=U95, color=factor(h2), shape=sparse, linetype=sparse) +
  geom_errorbar(width=0, linetype="solid") +
  geom_line() + 
  geom_point(fill="white") +
  scale_shape_manual("Sparse model:", values=c("TRUE"=21, "FALSE"=19)) +
  scale_linetype_manual("Sparse model:", values=c("TRUE"=2, "FALSE"=1)) +
  scale_color_manual("Heritability (h2):", values=c("#c7e9b4", "#41b6c4", "#225ea8", "#081d58")) +
  scale_y_continuous("AUC (95% CI)", limits=range(c(all_perf$L95, all_perf$U95))) +
  xlab("Proportion of causal variations (p)") +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8),
    legend.position="bottom", legend.text=element_text(size=6), legend.title=element_text(size=7),
    legend.box="vertical"
  )

# Create plot showing lasso-sum paramaeter tuning
g3 <- ggplot(lasso_perf) +
  aes(x=lambda, y=estimate, ymin=L95, ymax=U95, color=factor(delta)) +
  geom_errorbar(width=0) +
  geom_line() +
  geom_point(fill="white", shape=21) +
  scale_color_manual("Delta", values=c("#9ebcda", "#8c96c6", "#8c6bb1", "#88419d", "#810f7c", "#4d004b")) +
  scale_y_continuous("AUC (95% CI)", limits=range(c(all_perf$L95, all_perf$U95))) +
  xlab("Lambda") +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8),
    legend.position="bottom", legend.text=element_text(size=6), legend.title=element_text(size=7)
  )

# Combine into single plot
g <- plot_grid(g1, g2, g3, nrow=1, align="h", axis="b")
ggsave(g, width=7.2, height=3, file="output/paper_outputs/figures/T2D_Exome_LDpred2_tuning.pdf")


##################################################################
# Compare AUC and OR for all trained PRSs
##################################################################

# Load curated LDpred2 performance results
ldpred_perf <- fread("output/ldpred2/hyperparam_selection/best_hyperparams.txt")
ldpred_perf <- ldpred_perf[case_def == "combined_basic_adjudicated"]

# Extract ORs and AUCs in wide format for  LDpred2 model and parameters that maximised AUC
ldpred_perf <- dcast(ldpred_perf[(best_model)], gwas + ldpred2_model + pgs ~ metric, value.var=c("estimate", "L95", "U95", "pval"))
ldpred_perf <- ldpred_perf[, .(gwas, ldpred2_model, parameters=pgs, OR=estimate_OR, OR.L95=L95_OR, OR.U95=U95_OR, OR.pval=pval_OR, AUC=estimate_AUC, AUC.L95=L95_AUC, AUC.U95=U95_AUC)]

# Determine ordering of GWASs by AUC:
ldpred_perf <- ldpred_perf[order(-AUC)]
ldpred_perf[, gwas := factor(gwas, levels=gwas)]

# Create plot showing both OR and AUC
ggdt <- rbind(idcol="metric",
  "AUC"=ldpred_perf[,.(gwas, estimate=AUC, L95=AUC.L95, U95=AUC.U95)],
  "Odds Ratio"=ldpred_perf[,.(gwas, estimate=OR, L95=OR.L95, U95=OR.U95)]
)

# Plot AUC
g <- ggplot(ggdt) +
  aes(x=gwas, y=estimate, ymin=L95, ymax=U95) +
  facet_grid(metric ~ ., scales="free_y") +
  geom_hline(data=data.table(metric="Odds Ratio", yintercept=1), aes(yintercept=yintercept), linetype=2) +
  geom_point(shape=23, fill="white") +
  geom_errorbar(width=0) +
  scale_y_continuous("AUC, OR (95% CI)") +
  xlab("") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    axis.text.x=element_text(size=7, angle=90, hjust=1, vjust=0.5, color="black"), 
    axis.text.y=element_text(size=6),
    axis.title=element_text(size=8),
    strip.background=element_blank(),
    strip.text=element_blank()
  )
ggsave(g, width=7.2, height=3, units="in", file="output/paper_outputs/figures/component_PRS_AUC_and_OR.pdf")

#####################################
# Get details for formatting table
#####################################

# Get number of SNPs in each PRS
ldpred_perf[, SNPs := NA]
n_snps <- foreach(this_gwas = ldpred_perf$gwas, .combine=c, .inorder=TRUE) %dopar% {
  this_pgs <- ldpred_perf[gwas == this_gwas, parameters]
  if (this_pgs == "pgs_auto_final2") {
    this_snps <- fread(sprintf("output/ldpred2/test/%s/ldpred2_auto_final2_varweights_combined_basic_adjudicated.txt.gz", this_gwas))
  } else {
    this_snps <- fread(sprintf("output/ldpred2/train/%s/ldpred2_pgs_varweights.txt.gz", this_gwas),
      select=c("rsid", "chr", "pos", "effect_allele", "other_allele", gsub("^pgs_", "beta_", this_pgs)), tmpdir="tmp")
  }
  setnames(this_snps, c("rsid", "chr", "pos", "effect_allele", "other_allele", "weight"))
  this_snps[weight != 0, .N]
}
ldpred_perf[, SNPs := n_snps]

# Recover information about parameters for optimal model
for (this_gwas in ldpred_perf$gwas) {
  this_model <- ldpred_perf[gwas == this_gwas, ldpred2_model]
  
  ldsc_results <- readRDS(sprintf("output/ldpred2/train/%s/ldsc_results.rds", this_gwas))
  h2_est <- ldsc_results["h2"]
  if (this_model == "infinitesimal") {
    ldpred_perf[gwas == this_gwas, parameters := sprintf("LD score regression heritability (h2): %s", format(h2_est, digits=3))]
  } else if (this_model == "grid-search") {
    h2_seq <- round(h2_est * c(0.3, 0.7, 1, 1.4), 4)
    p_seq <- signif(bigsnpr::seq_log(1e-5, 1, length.out = 21), 2)
    params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))

    param_row <- ldpred_perf[gwas == this_gwas, as.integer(gsub(".*_", "", parameters))]
    params <- params[param_row,]

    ldpred_perf[gwas == this_gwas, parameters := sprintf("Proportion of causal variants (p): %.3f, heritability captured by variants (h2): %.3f, Sparse model: %s",
      round(params$p, digits=3), round(params$h2, digits=3), params$sparse)]
  } else if (this_model == "automatic") {
    params <- fread(sprintf("output/ldpred2/train/%s/ldpred2_auto_parameters.txt", this_gwas))    
    ldpred_perf[gwas == this_gwas, parameters := sprintf("Proportion of causal variants (p): %.3f, heritability captured by variants (h2): %.3f", 
      params[(pass_chain_filter), round(mean(p_est), digits=3)], params[(pass_chain_filter), round(mean(h2_est), digits=3)])]
  } else if (this_model == "lasso-sum") {
    params <- fread(sprintf("output/ldpred2/train/%s/ldpred2_lassosum2_parameters.txt", this_gwas))
    param_row <- ldpred_perf[gwas == this_gwas, as.integer(gsub(".*_", "", parameters))]
    params <- params[paramset == param_row]

    ldpred_perf[gwas == this_gwas, parameters := sprintf("Delta: %s, Lambda: %.3f", prettyNum(params$delta), params$lambda)]
  }
}

# Order table for excel
ldpred_perf <- ldpred_perf[order(-AUC), .(gwas, ldpred2_model, parameters, SNPs, OR, OR.L95, OR.U95, OR.pval, AUC, AUC.L95, AUC.U95)]
fwrite(ldpred_perf, sep="\t", quote=FALSE, file="output/paper_outputs/tables/component_PRS_details.txt")


