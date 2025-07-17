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

# What ancestry are we currently extracting the genotype data for? 
# What GWAS are we running LDpred2 on?
ancestry <- "EUR"
gwas <- "T2D_MVP_Multi"

# Load in the curated GWAS information so we can determine whether this GWAS is
# for a continuous or binary trait
gwas_info <- fread("data/gwas_summary_stats/curated_gwas_table.txt", na.strings=c("-"))
gwas_type <- ifelse(gwas_info[PRS == gwas, is.na(Cases)], "continuous", "binary")

#########
# Setup
#########

# Create output directory
outdir <- sprintf("output/ldpred2/train/%s/%s", ancestry, gwas)
system(sprintf("mkdir -p %s", outdir), wait=TRUE)

# Set path to harmonized and filtered GWAS summary statistics 
gwas_file <- sprintf("data/filtered_sumstats/filtered_gwas/%s.txt.gz", gwas)

# Set up temporary directory - clean up if already exists
tmpdir <- sprintf("tmp/ldpred2/%s/%s", ancestry, gwas)
if (dir.exists(tmpdir)) system(sprintf("rm -rf %s", tmpdir), wait=TRUE)
system(sprintf("mkdir -p %s", tmpdir), wait=TRUE)

# Set up parallism 
nCores <- 22 # N.b. 'nCores' variable is expected elsewhere in script
doMC::registerDoMC(nCores)

###############################################
# Do per-SNP QC of GWAS summary statistics
############################################

# Load the gwas summary statistics, match to the genotype data, and obtain per-SNP standard deviations and 
# allele frequencies for downstream SNP QC.
gwas_ss <- fread(gwas_file)
setnames(gwas_ss, "AoU_varID", "rsid")
gwas_ss <- foreach(this_chr = 1:22, .combine=rbind) %dopar% {
	# Attach file-backed genotype data
  geno <- snp_attach(sprintf("data/ldpred2/%s/filtered_chr%s.bed", ancestry, this_chr))

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
	matched_snps <- snp_match(as.data.frame(gwas_ss[chr == this_chr]), map, strand_flip=FALSE)
	setDT(matched_snps)

	# Obtain allele frequencies in the training data - dosages count 'a0'
	matched_snps[, a1freq := Matrix::colSums(geno$genotype[, `_NUM_ID_`], na.rm=TRUE) / (Matrix::colSums(!is.na(geno$genotype[, `_NUM_ID_`]))*2)]

	# Compute the standard deviation of the allele frequency
	# See https://privefl.github.io/bigsnpr-extdoc/polygenic-scores-pgs.html
	matched_snps[, sd_val := sqrt(2 * a1freq * (1 - a1freq))]

	# Return
	return(matched_snps)
}

# Do Per SNP QC of the summary stats.
# See https://privefl.github.io/bigsnpr-extdoc/polygenic-scores-pgs.html and
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8016455/
if (gwas_type == "binary") {
	gwas_ss[sd_val > 0, sd_ss := 2 / (beta_se * sqrt(n_eff))]
} else {
	gwas_ss[sd_val > 0, sd_y_est := median(sd_val * beta_se * sqrt(n_eff))]
	gwas_ss[sd_val > 0, sd_ss := sd_y_est / (beta_se * sqrt(n_eff))]
} 
gwas_ss[, fail_qc := sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05]

# Write out GWAS QC details
if (gwas_type == "continuous") {
	gwas_qc <- gwas_ss[, .(AoU_varID=rsid.ss, chr, pos, effect_allele=a1, other_allele=a0, gwas_beta=beta, gwas_se=beta_se, 
												 n_eff, trainingset_EAF=a1freq, sd_val, sd_y_est, sd_ss, fail_qc)]
} else {
	gwas_qc <- gwas_ss[, .(AoU_varID=rsid.ss, chr, pos, effect_allele=a1, other_allele=a0, gwas_beta=beta, gwas_se=beta_se, 
												 n_eff, trainingset_EAF=a1freq, sd_val, sd_ss, fail_qc)]
}
fwrite(gwas_qc, sep="\t", quote=FALSE, compress="gzip", file=sprintf("%s/ldpred2_gwasqc.txt.gz", outdir))

# Diagnostic plot
g <- ggplot(gwas_qc, aes(x=sd_val, y=sd_ss, color=fail_qc)) +
	theme_bigstatsr() +
	geom_point(shape=19, size=0.5, alpha=0.5) +
	geom_abline(intercept=0, slope=1, linetype=2, colour="red") +
	scale_colour_manual(name="SNP failed ldpred2 QC", values=c("TRUE"="purple", "FALSE"="yellow")) +
	xlab(sprintf("SNP dosage SD in AoU %s ancestry", ancestry)) + 
	ylab(sprintf("SNP dosage SD in %s GWAS", gwas)) +
	theme(legend.position="bottom")
ggsave(g, width=7.2, height=6, file=sprintf("%s/ldpred2_gwasqc.png", outdir))

cat(sprintf("%s of %s (%s%%) SNPs failed LDpred2 QC.\n", format(gwas_ss[(fail_qc), .N], big.mark=","), 
						format(gwas_ss[, .N], big.mark=","), round(gwas_ss[(fail_qc), .N]/gwas_ss[, .N]*100*100)/100),
						file=sprintf("%s/ldpred2_gwasqc_fail_rate.txt", outdir))

# Error if no variants pass QC 
if (all(gwas_ss$fail_qc)) {
	system(sprintf("rm -rf %s", tmpdir), wait=TRUE)
	stop(sprintf("No variants passing LDpred2 QC for GWAS %s in AoU %s training data"), gwas, ancestry)
}

#######################################################
# Load and aggregate pre-computed SNP-SNP correlations
#######################################################

# Load correlation matrix for each chromosome, filter to variants passing QC,
# compute LD score, and aggregate into single large sparse file-backed big
# matrix
for (this_chr in 1:22) {
	cat("Loading correlation matrix for chromosome ", this_chr, "\n")
	# Load correlation matrix precomputed for all candidate variants
  corr0 <- readRDS(sprintf("data/ldpred2/%s/filtered_chr%s_ldcorr.rds", ancestry, this_chr))

	# Filter to those passing QC for this GWAS
	corr0 <- gwas_ss[chr == this_chr & !(fail_qc), corr0[`_NUM_ID_`, `_NUM_ID_`]]

	# Compute LD score
	gwas_ss[chr == this_chr & !(fail_qc), LDsum := Matrix::colSums(corr0^2)]

	# Aggregate into a single sparse big matrix
	if (this_chr == 1) {
		cat("Initialized SFBM\n")
		genocorr <- as_SFBM(corr0, backingfile=sprintf("%s/ldcorr_passqc", tmpdir), compact = TRUE)
	} else {
		cat("Adding matrix to SFBM\n")
		genocorr$add_columns(corr0, nrow(genocorr))
	}
}

##############################
# Run infinitesimal model
##############################
cat("Running infinitesimal model\n")
ldsc <- gwas_ss[!(fail_qc), snp_ldsc(
	ld_score = LDsum, ld_size = .N, 
	chi2 = (beta / beta_se)^2,
	sample_size = n_eff,
	ncores = nCores 
)]

# Write out LDSC results
saveRDS(ldsc, sprintf("%s/ldsc_results.rds", outdir))

# Write out model parameters (in this case, just estimated heritability)
inf_params <- data.table(h2=ldsc[["h2"]])
fwrite(inf_params, sep="\t", quote=FALSE, file=sprintf("%s/infinitesimal_model_parameters.txt", outdir))

###########################
# Run Grid model
###########################
cat("Running grid model\n")
h2_est <- ldsc[["h2"]]
h2_seq <- round(h2_est * c(0.3, 0.7, 1, 1.4), 4)
p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2)
params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))

beta_grid <- snp_ldpred2_grid(genocorr, gwas_ss[!(fail_qc)], params, ncores = nCores)

# Write out model parameters
setDT(params)
fwrite(params, sep="\t", quote=FALSE, file=sprintf("%s/grid_model_parameters.txt", outdir))

#######################
# Run auto model
#######################
cat("Running auto model\n")

multi_auto <- snp_ldpred2_auto(
	genocorr, gwas_ss[!(fail_qc)], h2_init = h2_est,
	vec_p_init = seq_log(1e-4, 0.5, length.out = 30),
	ncores = nCores
)

# check for "chain" convergence
auto_params <- rbindlist(lapply(multi_auto, function(x) {
	data.table(p_init = x$p_init, h2_init = x$h2_init, p_est = x$p_est, h2_est = x$h2_est)
}))
auto_params[, paramset := .I]

auto_path <- foreach(pIdx = seq_along(multi_auto), .combine=rbind) %do% {
	auto = multi_auto[[pIdx]]
	data.table(paramset = pIdx, path_iter = seq_along(auto$path_p_est), 
						 p_est = auto$path_p_est, h2_est = auto$path_h2_est)
}

g1 <- ggplot(auto_path) + aes(x = path_iter, y=p_est) +
	theme_bigstatsr() + 
	geom_hline(data = auto_params, aes(yintercept=p_est), col="blue") +
	geom_point(shape=19, size=0.5) +
	scale_y_log10(name="p") + xlab("") +
	facet_wrap(~ paramset, ncol=10, labeller = label_both) + 
	theme(strip.background=element_blank(), strip.text=element_text(size=6), 
				axis.text=element_text(size=6), axis.title=element_text(size=10))

g2 <- ggplot(auto_path) + aes(x = path_iter, y=h2_est) +
	theme_bigstatsr() + 
	geom_hline(data = auto_params, aes(yintercept=h2_est), col="blue") +
	geom_point(shape=19, size=0.5) +
	ylab("h2") + xlab("") +
	facet_wrap(~ paramset, ncol=10, labeller = label_both) +
	theme(strip.background=element_blank(), strip.text=element_text(size=6), 
				axis.text=element_text(size=6), axis.title=element_text(size=10))

g <- plot_grid(g1, g2, nrow=2) 
ggsave(g, width=20, height=10, units="in", file=sprintf("%s/auto_model_chain_convergence.png", outdir))

# Extract and write out parameters
auto_params[, paramset := NULL]
fwrite(auto_params, sep="\t", quote=FALSE, file=sprintf("%s/auto_model_parameters.txt", outdir))

##########################
# Run lassosum2
##########################
cat("Running lassosum2 model\n")
beta_lassosum2 <- snp_lassosum2(genocorr, gwas_ss[!(fail_qc)], ncores = nCores)

# Extract and write out model parameters
params <- attr(beta_lassosum2, "grid_param")
setDT(params)
fwrite(params, sep="\t", quote=FALSE, file=sprintf("%s/lassosum2_model_parameters.txt", outdir))

####################################################################
# Write out file with all LDpred2 per-SNP betas required to compute
# PRSs for hyperparameter selection in the metaPRS training samples
####################################################################

pgs_betas <- gwas_ss[!(fail_qc), .(AoU_varID=rsid.ss, chr, pos, effect_allele=a1, other_allele=a0)]

if ("beta_inf" %in% names(gwas_ss)) { 
	pgs_betas <- cbind(pgs_betas, gwas_ss[!(fail_qc), .(infinitesimal=beta_inf)])
}

beta_grid <- as.data.table(beta_grid)
setnames(beta_grid, gsub("V", "grid_", names(beta_grid)))
pgs_betas <- cbind(pgs_betas, beta_grid)

beta_auto <- as.data.table(beta_auto)
setnames(beta_auto, gsub("V", "auto_", names(beta_auto)))
pgs_betas <- cbind(pgs_betas, beta_auto)

beta_lassosum2 <- as.data.table(beta_lassosum2)
setnames(beta_lassosum2, gsub("V", "lassosum2_", names(beta_lassosum2)))
pgs_betas <- cbind(pgs_betas, beta_lassosum2)

fwrite(pgs_betas, sep="\t", quote=FALSE, compress="gzip", file=sprintf("%s/ldpred2_pgs_varweights.txt.gz", outdir))

#######################
# Cleanup
#######################
system(sprintf("rm -rf %s", tmpdir), wait=TRUE)

