# Code adapted from LDpred2 tutorial https://privefl.github.io/bigsnpr/articles/LDpred2.html
# Several key differences/modifications:
#
#  - genotype data are stored per-chromosome, so matching SNPs and computing PGS from
#    fitted betas requires looping
#
#  - SNP-SNP correlation matrix has been precomputed in previous script
#
#  - The 'auto', 'grid', and 'lassosum2' models all return betas/PGS for a variety of
#    of hyperparameters. In the tutorial they split their data into training and test
#    sets to determine the best hyperparameter. Here, we don't split the training set
#    and instead return betas for all hyperparameters so that we can then compute them
#    in the metaPRS training samples to  determine the best hyperparameter and model
#    the will then be fed into elasticnet to derive the ancestry-specific metaPRS




library(data.table)
library(foreach)
library(bigsnpr)
library(ggplot2)
library(cowplot)
library(stringr)
library(R.utils)

bigparallelr::set_blas_ncores(1)#Needed to allow multiple cores without error
options(bigstatsr.check.parallel.blas = FALSE) #Needed to allow multiple cores without error
#my_bucket <- Sys.getenv("WORKSPACE_BUCKET")


#Example inputs
# gwas_info_fname=/home/jupyter/workspaces/t2dmetaprs/data/gwas/curated_gwas_table.txt
# variants_to_exclude_fname=/home/jupyter/workspaces/t2dmetaprs/data/variant_lists/NA_variants_snp_cor.csv
# gwas_file_name=/home/jupyter/workspaces/t2dmetaprs/data/gwas/filtered_gwas/PCOS_2018.txt.gz
# genotype_pattern_rds=/home/jupyter/workspaces/t2dmetaprs/data/genotype/ldpred2/sas/*.rds
# root_output_dir=test_output
# nCores=8
#
# export gwas_info_fname
# export variants_to_exclude_fname
# export gwas_file_name
# export genotype_pattern_rds
# export root_output_dir
# export nCores





#Inputs:
nCores <- as.integer(Sys.getenv("nCores"))
gwas_info_fname=Sys.getenv("gwas_info_fname")
variants_to_exclude_fname=Sys.getenv("variants_to_exclude_fname")
gwas_file=Sys.getenv("gwas_file_name")
genotype_dir=dirname(Sys.getenv("genotype_pattern_rds"))
#Output:
root_output_dir=Sys.getenv("root_output_dir")

print(nCores)
print(gwas_info_fname)
print(variants_to_exclude_fname)
print(gwas_file)
print(genotype_dir)
print(root_output_dir)

#Set tmpdir
tmpdir="tmp_dir"
dir.create(tmpdir,recursive = T)

#Get ancestry and gwas labels
ancestry=basename(genotype_dir)
gwas=str_remove(basename(gwas_file),".txt.gz")
print(ancestry)
print(gwas)
#####Set output directories

out_dir_results_root=paste0(root_output_dir,"/results")
out_dir_qc_root=paste0(root_output_dir,"/qc")
out_dir_plot_root=paste0(root_output_dir,"/plots")

#Create directores
# Create output directory
outdir =sprintf("%s/%s/%s",out_dir_results_root, ancestry, gwas)
dir.create(outdir,recursive = T,showWarnings = T)

#Create qc directories
out_dir_qc_stats=sprintf("%s/stats/%s",out_dir_qc_root,ancestry)
out_dir_qc_summary=sprintf("%s/qc_summary/%s",out_dir_qc_root,ancestry)
out_dir_sd_y_est=sprintf("%s/sd_y_est/%s",out_dir_qc_root,ancestry)
dir.create(out_dir_sd_y_est,recursive = T,showWarnings = T)
dir.create(out_dir_qc_stats,recursive = T,showWarnings =T)
dir.create(out_dir_qc_summary,recursive = T,showWarnings =T)

#Create plot directories
out_dir_auto_model_chain_convergence_plot=sprintf("%s/auto_model_chain_conv/%s",out_dir_plot_root,gwas)
out_dir_qc_plots=sprintf("%s/qc_plots/%s",out_dir_plot_root,gwas)
dir.create(out_dir_qc_plots,recursive = T,showWarnings = T)
dir.create(out_dir_auto_model_chain_convergence_plot,recursive = T,showWarnings = T)

#Print output dirs for debugging
print(outdir)
print(out_dir_qc_stats)
print(out_dir_sd_y_est)
print(out_dir_auto_model_chain_convergence_plot)
print(out_dir_qc_plots)

# Load in the curated GWAS information so we can determine whether this GWAS is
# for a continuous or binary trait
gwas_info <- fread(gwas_info_fname, na.strings=c("-"))
gwas_type <- ifelse(gwas_info[PRS == gwas, is.na(Cases)], "continuous", "binary")

#########
# Setup
#########

# Create output directory

cat("Starting",gwas,ancestry,format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"\n")

# Set up parallism
#doMC::registerDoMC(1) #Multitasking can fail with just a warning which is dangerous

###############################################
# Do per-SNP QC of GWAS summary statistics
#
# See:
# https://privefl.github.io/bigsnpr/articles/LDpred2.html#quality-control-of-gwas-summary-statistics
# https://github.com/privefl/paper-misspec/blob/main/code/prepare-sumstats/CAD.R
# https://github.com/privefl/paper-misspec/blob/main/code/prepare-sumstats/vitaminD.R
# https://pubmed.ncbi.nlm.nih.gov/36105883/
# https://privefl.github.io/bigsnpr-extdoc/polygenic-scores-pgs.html
###############################################

# Load the gwas summary statistics, match to the genotype data, and obtain per-SNP standard deviations and
# allele frequencies for downstream SNP QC.
gwas_ss <- fread(gwas_file)
setnames(gwas_ss, "AoU_varID", "rsid")
setnames(gwas_ss, "a1freq", "a1freq_gwas")

# Add flag column noting QC exclusion reason
gwas_ss[, snp_qc := "Passes LDpred2 SNP QC"] # placeholder, overwritten on QC fail

#Load the ancestry-specific variants to omit
variants_to_exclude=fread(variants_to_exclude_fname)
setnames(variants_to_exclude,"ancestry","variants_to_exclude_ancestry")
variants_to_exclude=variants_to_exclude[ancestry==variants_to_exclude_ancestry,]
gwas_ss[rsid %in% variants_to_exclude$AoU_varID, snp_qc := "Did not pass ancestry-specific MAF threshold"]

cat("Matching SNPS",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"\n")

##This fast enough that dosen't need to be paralelized probably
gwas_ss <- foreach(this_chr = unique(gwas_ss$chr), .combine=rbind) %do% { #%dopar% with %do% if you don't want to parallelise it
  # Attach file-backed genotype data
  geno <- snp_attach(sprintf("%s/filtered_chr%s.rds",genotype_dir, this_chr))

  # Match summary stats to genotype data. A few notes:
  # - Summary stats for all GWAS have already been filtered to a common variant
  #   set that intersects with the variants with MAF > 1% in any AoU ancestry
  # - Strand orientation of alleles has also already been harmonized to AoU for
  #   all GWAS
  # - Effect alleles have also all been harmonized so that the effect allele is
  #   the second allele in the AoU variant ID
  # - The extracted dosages in snp_readBed (and subsequent snp_attach here)
  #   correspond to counts of 'a1'.
  # - In 'matched_snps' the 'beta' estimates are harmonized to the 'a1' column
  # - The '_NUM_ID_' column corresponds to the column index of the SNP in the
  #   genotype data, while '_NUM_ID_.ss' corresponds to the row index in the
  #   summary stats file
  map <- geno$map[-2]
  names(map) <- c("chr", "rsid", "pos", "a1", "a0")
  matched_snps <- snp_match(as.data.frame(gwas_ss[chr == this_chr,drop=F]), map, strand_flip=FALSE)
  setDT(matched_snps)

  # Obtain allele frequencies in the training data - dosages count 'a0'
  matched_snps[, a1freq_AoU := big_colstats(geno$genotypes, ind.col=`_NUM_ID_`)$sum / (2 * nrow(geno$genotypes))]

  # Chromosomes with only 1 SNP cause problems in downstream code, so are excluded.
  if(nrow(matched_snps) == 1) {
    matched_snps[, snp_qc := "Only SNP on chromosome"]
  }

  return(matched_snps)
}

#For debugging
print(head(gwas_ss))

# First, remove low-powered SNPs in the GWAS.
# Filter at 70% of max effective sample size where per-SNP n_eff reported
# If n_eff is for the overall study, no SNPs will be excluded
n_eff_threshold <- gwas_ss[, 0.7 * max(n_eff, na.rm=TRUE)]
gwas_ss[snp_qc == "Passes LDpred2 SNP QC" & # Only apply each QC filter to SNPs not already excluded by previous steps
          n_eff < n_eff_threshold, snp_qc := "Low power SNP in GWAS (N_eff < 70% max N_eff)"]

# Check whether standard deviations of genotypes in the GWAS summary statistics
# are consistent with the ones in All of Us, filtering out SNPs which are too
# different
gwas_ss[, sd_AoU := sqrt(2 * a1freq_AoU * (1 - a1freq_AoU))] # Standard deviation of Allele frequency in All of Us

# Infer standard deviation of the genotypes in the gwas
if (gwas_type == "binary") {
  sd_y_est <- 2
} else {
  sd_y_est <- gwas_ss[, quantile(sqrt(0.5 * (n_eff * beta_se^2 + beta^2)), 0.01)]
  cat(sprintf("SD(y) estimated as: %0.3f\n", sd_y_est), file=sprintf("%s/%s__%s__ldpred2_gwasqc_sd_y_est.txt", out_dir_sd_y_est,gwas,ancestry))
}
gwas_ss[, sd_gwas := sd_y_est / sqrt(n_eff * beta_se^2 + beta^2)]

# Remove SNPs where the estimated genotype dosage SD in the GWAS diverges too far from the SD of the AoU allele frequencies
gwas_ss[snp_qc == "Passes LDpred2 SNP QC" & sd_gwas < 0.05, snp_qc := "Low MAF in GWAS (Estimated genotype dosage SD < 0.05)"]
gwas_ss[snp_qc == "Passes LDpred2 SNP QC" & sd_AoU < 0.05, snp_qc := "Low MAF in All of Us (Genotype dosage SD < 0.05)"]
gwas_ss[snp_qc == "Passes LDpred2 SNP QC" & (sd_gwas < (0.5 * sd_AoU) | sd_gwas > (sd_AoU + 0.1)),
        snp_qc := "Estimated genotype dosage SD in GWAS too different from genotype dosage in All of Us"]

# Where the GWAS summary statistics include allele frequencies, also remove SNPs that are too divergent based
# directly on the observed allele frequencies
gwas_ss[, sd_gwas_af := sqrt(2 * a1freq_gwas * (1 - a1freq_gwas))]
gwas_ss[, a1freq_diff := a1freq_AoU - a1freq_gwas]
gwas_ss[snp_qc == "Passes LDpred2 SNP QC" & sd_gwas_af < 0.05, snp_qc := "Low MAF in GWAS (Reported genotype dosage SD < 0.05)"]
gwas_ss[snp_qc == "Passes LDpred2 SNP QC" & (sd_gwas < (0.7 * sd_gwas_af) | sd_gwas > (sd_gwas_af + 0.1)),
        snp_qc := "Estimated genotype dosage SD in GWAS diverged from SD calculated from reported allele frequencies"]
gwas_ss[snp_qc == "Passes LDpred2 SNP QC" & abs(a1freq_diff) > 0.07,
        snp_qc := "GWAS allele frequency too divergent from All of Us allele frequency (abs diff > 0.07)"]

# Order the SNP QC messages by their sequence for plotting/tabulation
gwas_ss[, snp_qc := factor(snp_qc, levels=c(
  "Passes LDpred2 SNP QC",
  "Did not pass ancestry-specific MAF threshold",
  "Only SNP on chromosome",
  "Low power SNP in GWAS (N_eff < 70% max N_eff)",
  "Low MAF in GWAS (Estimated genotype dosage SD < 0.05)",
  "Low MAF in All of Us (Genotype dosage SD < 0.05)",
  "Estimated genotype dosage SD in GWAS too different from genotype dosage in All of Us",
  "Low MAF in GWAS (Reported genotype dosage SD < 0.05)",
  "Estimated genotype dosage SD in GWAS diverged from SD calculated from reported allele frequencies",
  "GWAS allele frequency too divergent from All of Us allele frequency (abs diff > 0.07)"
))]

# Tabulate QC exclusion reasons:
gwas_qc_stats <- gwas_ss[, .(SNPs=.N, pct=.N/gwas_ss[,.N]*100), by=snp_qc]
gwas_qc_stats <- gwas_qc_stats[order(snp_qc)]
fwrite(gwas_qc_stats, sep="\t", quote=FALSE, file=sprintf("%s/%s__%s__ldpred2_gwasqc_summary.txt", out_dir_qc_summary,gwas,ancestry))

# Write out GWAS QC details
#gwas_qc <- gwas_ss[, .(AoU_varID=rsid.ss, chr, pos, effect_allele=a1, other_allele=a0, gwas_beta=beta, gwas_se=beta_se,
#                       n_eff, gwas_eaf=a1freq_gwas, AoU_eaf=a1freq, eaf_diff=a1freq_diff, sd_gwas, sd_AoU, sd_gwas_af, snp_qc)]

gwas_qc <- gwas_ss[, .(AoU_varID=rsid.ss, chr, pos, effect_allele=a1, other_allele=a0, gwas_beta=beta, gwas_se=beta_se,
                       n_eff, gwas_eaf=a1freq_gwas, AoU_eaf=a1freq_AoU, eaf_diff=a1freq_diff, sd_gwas, sd_AoU, sd_gwas_af, snp_qc)]

fwrite(gwas_qc, sep="\t", quote=FALSE, compress="gzip", file=sprintf("%s/%s__%s__ldpred2_gwasqc.txt.gz", out_dir_qc_stats,gwas,ancestry))

# Plot histogram of per SNP effective sample size used in first QC step
g <- ggplot(gwas_ss, aes(x=n_eff)) +
  geom_histogram() + #geom_hist() +
  geom_vline(xintercept=n_eff_threshold, linetype=2, color="red") +
  xlab("Per-SNP effective sample size") +
  ylab("Count") +
  theme_bw() +
  theme(
    axis.title=element_text(size=8), axis.text=element_text(size=6)
  )
ggsave(g, width=7.2, height=3.5, file=sprintf("%s/%s__%s__ldpred2_gwasqc_n_eff_hist.png", out_dir_qc_plots, gwas, ancestry))

# Scatterplot comparing the SD of genotype dosages between the GWAS and All of Us
g <- ggplot(gwas_qc, aes(x=sd_AoU, y=sd_gwas, color=snp_qc)) +
  geom_point(shape=19, size=0.5, alpha=0.5) +
  geom_abline(intercept=0, slope=1, linetype=2, colour="red") +
  scale_colour_manual(name="LDpred2 SNP QC", values=c(
    "Passes LDpred2 SNP QC"="black",
    "Did not pass ancestry-specific MAF threshold"="#fdbf6f",
    "Only SNP on chromosome"="#e31a1c",
    "Low power SNP in GWAS (N_eff < 70% max N_eff)"="#cab2d6",
    "Low MAF in GWAS (Estimated genotype dosage SD < 0.05)"="#b2df8a",
    "Low MAF in All of Us (Genotype dosage SD < 0.05)"="#ff7f00",
    "Estimated genotype dosage SD in GWAS too different from genotype dosage in All of Us"="#a6cee3",
    "Low MAF in GWAS (Reported genotype dosage SD < 0.05)"="#33a02c",
    "Estimated genotype dosage SD in GWAS diverged from SD calculated from reported allele frequencies"="#1f78b4",
    "GWAS allele frequency too divergent from All of Us allele frequency (abs diff > 0.07)"="#253494"
  )) +
  guides(colour = guide_legend(ncol = 1)) +
  xlab(sprintf("SNP dosage SD in AoU %s ancestry", ancestry)) +
  ylab(sprintf("SNP dosage SD in %s GWAS", gwas)) +
  theme_bw() +
  theme(
    axis.title=element_text(size=8), axis.text=element_text(size=6),
    legend.title=element_text(size=8), legend.text=element_text(size=6)
  )
ggsave(g, width=7.2, height=4, file=sprintf("%s/%s__%s__ldpred2_gwasqc_sd_dosage_compare.png", out_dir_qc_plots,gwas,ancestry))

# Scatterplot comparing the SD of genotype dosages derived from summary stats vs. allele frequencies
g <- ggplot(gwas_qc, aes(x=sd_gwas_af, y=sd_gwas, color=snp_qc)) +
  geom_point(shape=19, size=0.5, alpha=0.5) +
  geom_abline(intercept=0, slope=1, linetype=2, colour="red") +
  scale_colour_manual(name="LDpred2 SNP QC", values=c(
    "Passes LDpred2 SNP QC"="black",
    "Did not pass ancestry-specific MAF threshold"="#fdbf6f",
    "Only SNP on chromosome"="#e31a1c",
    "Low power SNP in GWAS (N_eff < 70% max N_eff)"="#cab2d6",
    "Low MAF in GWAS (Estimated genotype dosage SD < 0.05)"="#b2df8a",
    "Low MAF in All of Us (Genotype dosage SD < 0.05)"="#ff7f00",
    "Estimated genotype dosage SD in GWAS too different from genotype dosage in All of Us"="#a6cee3",
    "Low MAF in GWAS (Reported genotype dosage SD < 0.05)"="#33a02c",
    "Estimated genotype dosage SD in GWAS diverged from SD calculated from reported allele frequencies"="#1f78b4",
    "GWAS allele frequency too divergent from All of Us allele frequency (abs diff > 0.07)"="#253494"
  )) +
  guides(colour = guide_legend(ncol = 1)) +
  xlab(sprintf("SNP dosage SD derived from frequencies in %s GWAS", gwas)) +
  ylab(sprintf("SNP dosage SD derived from effect sizes in %s GWAS", gwas)) +
  theme_bw() +
  theme(
    axis.title=element_text(size=8), axis.text=element_text(size=6),
    legend.title=element_text(size=8), legend.text=element_text(size=6)
  )
ggsave(g, width=7.2, height=4, file=sprintf("%s/%s__%s__ldpred2_gwasqc_sd_from_dosage_vs_frequencies.png", out_dir_qc_plots,gwas,ancestry))

# Plot histogram of allele frequency differences
g <- ggplot(gwas_ss, aes(x=a1freq_diff)) +
  geom_histogram() +#geom_hist() +
  geom_vline(xintercept=0.07, linetype=2, color="red") +
  geom_vline(xintercept=-0.07, linetype=2, color="red") +
  xlab(sprintf("EAF in AoU %s ancestry - EAF in %s GWAS", ancestry, gwas)) +
  ylab("Count") +
  theme_bw() +
  theme(
    axis.title=element_text(size=8), axis.text=element_text(size=6)
  )
ggsave(g, width=7.2, height=3.5, file=sprintf("%s/%s__%s__ldpred2_gwasqc_n_eff_hist.png", out_dir_qc_plots, gwas, ancestry))

# Filter to variants passing qc
gwas_ss <- gwas_ss[snp_qc == "Passes LDpred2 SNP QC"]

# Error if no variants pass QC
if (gwas_ss[,.N] == 0) {
  system(sprintf("rm -rf %s", tmpdir), wait=TRUE)
  stop(sprintf("No variants passing LDpred2 QC for GWAS %s in AoU %s training data"), gwas, ancestry)
}

#######################################################
# Load and aggregate pre-computed SNP-SNP correlations
#######################################################

# Load correlation matrix for each chromosome, filter to variants passing QC,
# compute LD score, and aggregate into single large sparse file-backed big
# matrix
cat("Loading Correlation Data",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"\n")

for (this_chr in 1:22) {
  cat("Loading correlation matrix for chromosome ", this_chr, "\n")
  # Load correlation matrix precomputed for all candidate variants
  corr0 <- readRDS(sprintf("%s/filtered_chr%s_ldcorr.rds",genotype_dir, this_chr))

  # Filter to those passing QC for this GWAS
  corr0 <- gwas_ss[chr == this_chr, corr0[`_NUM_ID_`, `_NUM_ID_`]]
  if(is.null(dim(corr0))){next}
  if(any(is.na(corr0@x))){stop("Correlation Matrix has NA!")}
  # Compute LD score
  gwas_ss[chr == this_chr, LDsum := Matrix::colSums(corr0^2)]

  # Aggregate into a single sparse big matrix
  if (this_chr == 1) {
    cat("Initialized SFBM\n")
    #unlink("/home/jupyter/workspaces/t2dmetaprs/tmp_wdir/ldpred2//eur/T2D_MVP_Multi/ldcorr_passqc.sbk")
    genocorr <- as_SFBM(corr0, backingfile=sprintf("%s/ldcorr_passqc", tmpdir), compact = TRUE)
  } else {
    cat("Adding matrix to SFBM\n")
    genocorr$add_columns(corr0, nrow(genocorr))
  }
}

##############################
# Run infinitesimal model
##############################
# table(is.na(gwas_ss$LDsum))
#	#gwas_ss=gwas_ss[!is.na(LDsum),]

cat("Running infinitesimal model",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"\n")

# Compute LD score regression to estimate SNP-based heritability
ldsc <- gwas_ss[, snp_ldsc(
  ld_score = LDsum, ld_size = .N,
  chi2 = (beta / beta_se)^2,
  sample_size = n_eff,
  ncores = nCores
)]

# Write out LDSC results
saveRDS(ldsc, sprintf("%s/%s__%s__ldsc_results.rds", outdir,gwas,ancestry))

# Extract estimated heritability (used by other models downstream)
h2_est <- ldsc[["h2"]]

# Compute PGS per SNP weights assuming infinitesimal model
if (h2_est > 0) {
  beta_inf <- snp_ldpred2_inf(corr=genocorr, df_beta=gwas_ss, h2=h2_est)

  # Write out model parameters (in this case, just estimated heritability)
  inf_params <- data.table(h2=h2_est)
  fwrite(inf_params, sep="\t", quote=FALSE, file=sprintf("%s/%s__%s__infinitesimal_model_parameters.txt", outdir,gwas,ancestry))
}

###########################
# Run Grid model
###########################
if (h2_est > 0) {
  cat("Running grid model",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"\n")

  h2_seq <- round(h2_est * c(0.3, 0.7, 1, 1.4), 4)
  p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2)
  params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))

  beta_grid <- snp_ldpred2_grid(genocorr, gwas_ss, params, ncores = nCores)

  # Write out model parameters
  setDT(params)
  fwrite(params, sep="\t", quote=FALSE, file=sprintf("%s/%s__%s__grid_model_parameters.txt", outdir,gwas,ancestry))
}

##########################
# Run lassosum2
##########################
cat("Running lassosum2 model",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"\n")

beta_lassosum2 <- snp_lassosum2(genocorr, gwas_ss, ncores = nCores)

# Extract and write out model parameters
params <- attr(beta_lassosum2, "grid_param")
setDT(params)
fwrite(params, sep="\t", quote=FALSE, file=sprintf("%s/%s__%s__lassosum2_model_parameters.txt", outdir,gwas,ancestry))

###########################################################################
# Write out file with all LDpred2 per-SNP betas so far required to compute
# PRSs for hyperparameter selection in the metaPRS training samples
###########################################################################
cat("Writting Final Output",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"\n")

pgs_betas <- gwas_ss[, .(AoU_varID=rsid.ss, chr, pos, effect_allele=a1, other_allele=a0)]

if (h2_est > 0) {
  pgs_betas[, infinitesmial := beta_inf]

  beta_grid <- as.data.table(beta_grid)
  setnames(beta_grid, gsub("V", "grid_", names(beta_grid)))
  pgs_betas <- cbind(pgs_betas, beta_grid)
}

beta_lassosum2 <- as.data.table(beta_lassosum2)
setnames(beta_lassosum2, gsub("V", "lassosum2_", names(beta_lassosum2)))
pgs_betas <- cbind(pgs_betas, beta_lassosum2)

fwrite(pgs_betas, sep="\t", quote=FALSE, compress="gzip", file=sprintf("%s/%s__%s__ldpred2_pgs_varweights.txt.gz", outdir,gwas,ancestry))

cat("Done!",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"\n")

#######################
# Run auto model - done last as sometimes it takes too long and needs to be killed
#######################
if (h2_est > 0) {
  cat("Running auto model",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"\n")

  # Run through possible combinations of shrink_cor and use_MLE until >50% of
  # chains converge; see: https://github.com/privefl/bigsnpr/issues/554
  hyper_params <- data.table(expand.grid(use_MLE=c(TRUE, FALSE), shrink_corr=seq(0.95, 0.4, by=-0.05))) # Worst case, runs 24 times
  for (ii in hyper_params[,.I]) {
    cat(sprintf("Running with hyper-parameters use_MLE = %s and shrink_corr = %0.2f\n", hyper_params[ii, use_MLE], hyper_params[ii, shrink_corr]))

    finished <- withTimeout(timeout=9000, onTimeout = "silent", { # skip if takes more than 2.5 hours
      multi_auto <- snp_ldpred2_auto(
        genocorr, gwas_ss, h2_init = h2_est,
        vec_p_init = seq_log(1e-4, 0.5, length.out = 30),
        allow_jump_sign = FALSE,
        use_MLE = hyper_params[ii, use_MLE],
        shrink_corr = hyper_params[ii, shrink_corr],
        ncores = nCores
      )
    })
    if (is.null(finished)) {
      next
    }

    # Check how many chains converged
    range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
    keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)))
    converged <- length(keep) / length(range)

    # Plot convergence of all chains for debugging purposes
    auto_params <- foreach(pIdx = seq_along(multi_auto), .combine=rbind) %do% {
      x <- multi_auto[[pIdx]]
      data.table(chain = pIdx, p_init = x$p_init, h2_init = x$h2_init, p_est = x$p_est, h2_est = x$h2_est)
    }
    auto_params[, converged := chain %in% keep]

    auto_path <- foreach(pIdx = seq_along(multi_auto), .combine=rbind) %do% {
      auto = multi_auto[[pIdx]]
      data.table(chain = pIdx, path_iter = seq_along(auto$path_p_est),
                 p_est = auto$path_p_est, h2_est = auto$path_h2_est)
    }
    auto_path[, converged := chain %in% keep]

    g1 <- ggplot(auto_path) + aes(x = path_iter, y=p_est) +
      geom_hline(data = auto_params, aes(yintercept=p_est), col="blue") +
      geom_point(shape=19, size=0.5, aes(color=converged)) +
      scale_color_manual("Chain converged", values=c("TRUE"="black", "FALSE"="#a50f15")) +
      scale_y_log10(name="p") + xlab("") +
      facet_wrap(~ chain, ncol=10, labeller = label_both) +
      theme_bw() +
      theme(legend.position="bottom",
            legend.title=element_text(size=10), legend.text=element_text(size=8),
            strip.background=element_blank(), strip.text=element_text(size=6),
            axis.text=element_text(size=6), axis.title=element_text(size=10))

    g2 <- ggplot(auto_path) + aes(x = path_iter, y=h2_est) +
      geom_hline(data = auto_params, aes(yintercept=h2_est), col="blue") +
      geom_point(shape=19, size=0.5, aes(color=converged)) +
      scale_color_manual("Chain converged", values=c("TRUE"="black", "FALSE"="#a50f15")) +
      ylab("h2") + xlab("") +
      facet_wrap(~ chain, ncol=10, labeller = label_both) +
      theme_bw() +
      theme(legend.position="bottom",
            legend.title=element_text(size=10), legend.text=element_text(size=8),
            strip.background=element_blank(), strip.text=element_text(size=6),
            axis.text=element_text(size=6), axis.title=element_text(size=10))

    g <- plot_grid(g1, g2, nrow=2)
    ggsave(g, width=20, height=10, units="in",
           file=sprintf("%s/%s__%s__auto_model_chain_convergence_use_MLE__%s__shrink_corr__%.2f.png",
                        out_dir_auto_model_chain_convergence_plot, gwas, ancestry, hyper_params[ii, use_MLE], hyper_params[ii, shrink_corr])) # Carles: change outdir?

    # Write out extended parameters in case we need to debug
    fwrite(auto_params, sep="\t", quote=FALSE,
           file=sprintf("%s/%s__%s__auto_model_parameters__%s__shrink_corr__%.2f.txt",
                        outdir, gwas, ancestry, hyper_params[ii, use_MLE], hyper_params[ii, shrink_corr])) # Carles: change outdir?

    # If more than 50% chains converged, no need to decrease the shrink_corr or
    # set use_MLE to FALSE
    cat (sprintf("%d percent chains converged\n", round(converged*100)))
    if (converged > 0.5) {
      break
    }
  }

  if (converged > 0.5) {
    # Take mean of the converged chains to get the final auto model
    beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
    auto_params_final <- auto_params[(keep),.("converged"=converged, use_MLE=hyper_params[ii, use_MLE],
                                              shrink_corr=hyper_params[ii, shrink_corr],  p_est = mean(p_est), h2_est = mean(h2_est))]

    # Write out final parameter information
    fwrite(auto_params_final, sep="\t", quote=FALSE, file=sprintf("%s/%s__%s__auto_model_parameters.txt",outdir,gwas,ancestry))

    # Add the PGS betas to the output file
    pgs_betas <- fread(file=sprintf("%s/%s__%s__ldpred2_pgs_varweights.txt.gz", outdir,gwas,ancestry))
    pgs_betas[, auto := beta_auto]
    fwrite(pgs_betas, sep="\t", quote=FALSE, compress="gzip", file=sprintf("%s/%s__%s__ldpred2_pgs_varweights.txt.gz", outdir,gwas,ancestry))
  }
}

#######################
# Cleanup
#######################
#Not needed as the VM is deleted after
#system(sprintf("rm -rf %s", tmpdir), wait=TRUE)
