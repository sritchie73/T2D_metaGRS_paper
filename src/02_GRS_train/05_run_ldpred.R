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
#    (8,616 participants), and instead return betas for all hyperparameters. These will
#   Then be computed in the separate training set of 103,832 participants to determine
#   the best hyperparameter per model and best model, to be fed into elasticnet for
#   metaGRS derivation.

library(data.table)
library(foreach)
library(bigsnpr)
library(ggplot2)
library(cowplot)
source("src/functions/glm_test.R")

# Need to know which GWAS are continuous traits vs. binary traits:
cont.GWAS <- c("BMI", "DBP", "EduYears", "FastingGlucose", "FastingInsulin", "HBA1C",
               "HDL", "HipCircumference", "HipCircumferenceAdjBMI", "HOMAB", "HOMAIR",
               "LDL", "Leptin", "LeptinAdjBMI", "SBP", "SmokingCigsPerDay", "TChol",
               "Triglycerides", "WaistCircumference", "WaistCircumferenceAdjBMI",
               "WHR", "WHRadjBMI")
bin.GWAS <- c("Bipolar", "CAD", "CAD_Japan", "Schizophrenia", "SmokingStatus", "StrokeAS",
              "StrokeCES", "StrokeIS", "StrokeLAS", "StrokeSVS", "T2D_2017", "T2D_2018_no_UKB",
              "T2D_Africa", "T2D_EastAsia", "T2D_Exome", "T2D_Exome_Adj_BMI", "T2D_FinnGen_v6",
              "T2D_GERA", "T2D_Japan", "T2D_PAGE", "SmokingStatusMale", "SmokingStatusFemale")


# Make sure we're running on the appropriate compute node to make use of ramdisks partition
#
# File back shared memory obejcts can only effectively be used when using the /ramdisks/ partition
# as a temporary working directory, which is only available here on the skylake partitions.
# Trying to use the lustre filesystem means LDpred2 grinds to a halt unless run on a single core
# (and only one instance across the whole cluster) due to the consistency checks made by the lustre
# filesystem on shared memory objects.
if (!(Sys.getenv("SLURM_JOB_PARTITION") %like% "skylake")) {
	stop("Script must be run on compute node on skylake or skylake-himem partitions")
}

# Set up parallelisation
source("src/functions/par_setup.R")

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

bad <- setdiff(gwass, c(cont.GWAS, bin.GWAS))
if (length(bad) > 0) {
  stop("GWAS(s) ",  paste(paste0("'", bad, "'"), collapse=", "), 
       " not found in 'cont.GWAS' or 'bin.GWAS' - you need to add", 
       " it to one of these in the script so the script knows whether",
       " its a binary or continuous trait")
}

###############################################
# Loop through GWASs
###############################################
for (gwas in gwass) {
  #########
  # Setup
  #########
  # Clean up memory from previous GWAS 
  gc()

  # Extract gwas name and path to summary stats
	gwas <- basename(gwas) # in case specified via path
	gwas_file <- sprintf("data/filtered_sumstats/%s/gwas_filtered_oriented_SNPs.txt.gz", gwas)
	stopifnot(file.exists(gwas_file))

	# Setup output directory
	outdir <- sprintf("output/ldpred2/train/%s", gwas)
	system(sprintf("mkdir -p %s", outdir), wait=TRUE)
	
	# Set up temporary directories - clean up if already exists
  tmpdir <- sprintf("tmp/ldpred2/%s", gwas)
  if (dir.exists(tmpdir)) system(sprintf("rm -rf %s", tmpdir), wait=TRUE)
  system(sprintf("mkdir -p %s", tmpdir), wait=TRUE)

  ramdir <- sprintf("/ramdisks/ldpred2/%s", gwas)
  if (dir.exists(ramdir)) system(sprintf("rm -rf %s", ramdir), wait=TRUE)
  system(sprintf("mkdir -p %s", ramdir), wait=TRUE)

  # Copy genotype data to ramdisks
  system(sprintf("cp data/ldpred2/filtered_ukb_chr*.{rds,bk} %s", ramdir), wait=TRUE)

	###############################################
	# Do per-SNP QC of GWAS summary statistics
	############################################

	# Load the gwas_ss, match to the genotype data, and obtain per-SNP standard deviations and 
	# allele frequencies for downstream SNP QC.
  registerDoMC(nCores/2) # re-register parallel backend in case its been clobbered by LDpred2
	gwas_ss <- fread(gwas_file)
	gwas_ss <- foreach(this_chr = 1:22, .combine=rbind) %dopar% {
		# Attach file-backed genotype data
		geno <- snp_attach(sprintf("%s/filtered_ukb_chr%s.rds", ramdir, this_chr))

		# Match summary stats to genotype data. A few notes:
		# - Summary stats for all GWAS have already been filtered to a common variant
		#   set that intersects with the variants in the UKB training data 
		# - Strand orientation of alleles has also already been harmonized to UKB for
		#   all GWAS
		# - Effect alleles have also all been harmonized to the minor allele in UKB
		# - Below, "a1" and "a0" correspond to the first and second occurring alleles
		#   in the .bim files for UKB.
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
	if (gwas %in% bin.GWAS) {
		gwas_ss[sd_val > 0, sd_ss := 2 / (beta_se * sqrt(n_eff))]
	} else if (gwas %in% cont.GWAS) {
		gwas_ss[sd_val > 0, sd_y_est := median(sd_val * beta_se * sqrt(n_eff))]
		gwas_ss[sd_val > 0, sd_ss := sd_y_est / (beta_se * sqrt(n_eff))]
	} 
	gwas_ss[, fail_qc := sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05]

	# Write out GWAS QC details
	if(gwas %in% cont.GWAS) {
		gwas_qc <- gwas_ss[, .(rsid=rsid.ss, chr, pos, effect_allele=a1, other_allele=a0, gwas_beta=beta, gwas_se=beta_se, 
													 n_eff, trainingset_EAF=a1freq, sd_val, sd_y_est, sd_ss, fail_qc)]
	} else if (gwas %in% bin.GWAS) {
		gwas_qc <- gwas_ss[, .(rsid=rsid.ss, chr, pos, effect_allele=a1, other_allele=a0, gwas_beta=beta, gwas_se=beta_se, 
													 n_eff, trainingset_EAF=a1freq, sd_val, sd_ss, fail_qc)]
	}
	fwrite(gwas_qc, sep="\t", quote=FALSE, compress="gzip", file=sprintf("%s/ldpred2_gwasqc.txt.gz", outdir))

  # Diagnostic plot
  g <- ggplot(gwas_qc, aes(x=sd_val, y=sd_ss, color=fail_qc)) +
    theme_bigstatsr() +
    geom_point(shape=19, size=0.5, alpha=0.5) +
    geom_abline(intercept=0, slope=1, linetype=2, colour="red") +
    scale_colour_manual(name="SNP failed ldpred2 QC", values=c("TRUE"="purple", "FALSE"="yellow")) +
    xlab("SNP dosage SD in training dataset") + 
    ylab("SNP dosage SD in GWAS") +
    theme(legend.position="bottom")
  ggsave(g, width=7.2, height=6, file=sprintf("%s/ldpred2_gwasqc.png", outdir))

  cat(sprintf("%s of %s (%s%%) SNPs failed LDpred2 QC.\n", format(gwas_ss[(fail_qc), .N], big.mark=","), 
              format(gwas_ss[, .N], big.mark=","), round(gwas_ss[(fail_qc), .N]/gwas_ss[, .N]*100*100)/100),
              file=sprintf("%s/ldpred2_gwasqc_fail_rate.txt", outdir))

  # If all variants fail QC, can move to next iteration
	if (all(gwas_ss$fail_qc)) {
		system(sprintf("rm -rf %s", tmpdir), wait=TRUE)
		system(sprintf("rm -rf %s", ramdir), wait=TRUE)
    next
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
		corr0 <- readRDS(sprintf("data/ldpred2/filtered_ukb_ldcorr_chr%s.rds", this_chr))

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

  # Relocate backing file to /ramdisks/
  cat("Moving SFBM backing file to /ramdisks\n")
  system(sprintf("cp %s/ldcorr_passqc.sbk %s/", tmpdir, ramdir), wait=TRUE)
  system(sprintf("rm %s/ldcorr_passqc.sbk", tmpdir), wait=TRUE)
  system(sprintf("ln -s %s/ldcorr_passqc.sbk %s/", ramdir, tmpdir), wait=TRUE)

	###################################################################
	# Load phenotype data and setup T2D prediction model comparisons
	###################################################################
	pheno <- fread("data/ukb/collated_curated_data.txt")
	pheno <- pheno[(ldpred2_samples) & visit_index == 0]

	# Ensure samples are ordered the same as genotype data
	fam <- snp_attach(sprintf("%s/filtered_ukb_chr22.rds", ramdir))$fam
	setDT(fam)
	pheno <- pheno[fam[,.(sample.ID)], on = .(eid=sample.ID)]

	# Fit null model (no PRS) in the training data
	null_model <- glm.test(
		type_2_diabetes ~ genetic_sex + age + PC1 + PC2 + PC3 + PC4 + PC5 + 
		PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + 
		PC17 + PC18 + PC19 + PC20,
		event_col = "type_2_diabetes", data=pheno, skip.ci=TRUE
	)

	# Set up model for testing individual pgs
	mf <- "type_2_diabetes ~ scale(%s) + genetic_sex + age"
	mf <- paste(mf, "+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
	mf <- paste(mf, "+ PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20")

	# Start table of model performance - we only care about whole model performance (AUC)
	# and odds ratios for the individual PGS.
	model_perf <- data.table(model="No PGS", AUC=null_model$AUC[1])

	##############################
	# Run infinitesimal model
	##############################
  cat("Running infinitesimal model\n")
	ldsc <- gwas_ss[!(fail_qc), snp_ldsc(
		ld_score = LDsum, ld_size = .N, 
		chi2 = (beta / beta_se)^2,
		sample_size = n_eff,
		ncores = nCores # number of cores detected in src/functions/par_setup.R
	)]

  # Write out LDSC results
  saveRDS(ldsc, sprintf("%s/ldsc_results.rds", outdir))

  # Extract estimated heritability
	h2_est = ldsc[["h2"]]
  
  # Infinitesimal model only valid where h2 > 0
  if (h2_est >= 0) {
		# Extract weights for corresponding PGS:
		gwas_ss[!(fail_qc), beta_inf := snp_ldpred2_inf(
			corr = genocorr, 
			df_beta = gwas_ss[!(fail_qc)], 
			h2 = h2_est
		)]

		# Get PGS levels in phenotype data:
		pheno[, pgs_inf := 0]
		for (this_chr in 1:22) {
			geno <- snp_attach(sprintf("%s/filtered_ukb_chr%s.rds", ramdir, this_chr))
			geno <- snp_fastImputeSimple(geno$genotypes, ncores=nCores) # doesn't work with missing genotypes, so need to impute as median
			pheno[, pgs_inf := pgs_inf + big_prodVec(
				X = geno,
				y.col = gwas_ss[!(fail_qc) & chr == this_chr, beta_inf],
				ind.col = gwas_ss[!(fail_qc) & chr == this_chr, `_NUM_ID_`],
				ncores = nCores
			)]
		}

		# Test performance
		inf_perf <- glm.test(sprintf(mf, "pgs_inf"), "type_2_diabetes", pheno, skip.ci=TRUE)

		# Add whole model performance (AUC) and odds ratio for the PGS (discarding OR estimates for covariates)
		model_perf <- rbind(model_perf, inf_perf[coefficient %like% "scale",.(model="inf", OR, SE=logOR.SE, P.value, AUC)], fill=TRUE)
  } else {
    model_perf <- rbind(model_perf, data.table(model="inf"), fill=TRUE)
  }

	###########################
	# Run Grid model
	###########################
  cat("Running grid model\n")
	h2_seq <- round(h2_est * c(0.3, 0.7, 1, 1.4), 4)
	p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2)
	params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))

	beta_grid <- snp_ldpred2_grid(genocorr, gwas_ss[!(fail_qc)], params, ncores = nCores)

	# Get PGS levels for all models
	pgs_grid <- foreach(this_chr = 1:22, .combine=`+`) %do% {
		geno <- snp_attach(sprintf("%s/filtered_ukb_chr%s.rds", ramdir, this_chr))
		geno <- snp_fastImputeSimple(geno$genotypes, ncores=nCores) # doesn't work with missing genotypes, so need to impute as median
		big_prodMat(
			X = geno, 
			A.col = beta_grid[gwas_ss[!(fail_qc), which(chr == this_chr)],],
			ind.col = gwas_ss[!(fail_qc) & chr == this_chr, `_NUM_ID_`],
			ncores = nCores
		)
	}

	# Assess model performance
  registerDoMC(nCores/2) # big_prodMat clobbers existing backend
	pgs_grid_perf <- foreach(gi = seq_len(nrow(params)), .combine=rbind) %dopar% {
		this_pgs <- pgs_grid[, gi]
		if (all(is.na(this_pgs))) {
			return(data.table(model="grid", paramset=gi, OR=NA_real_, SE=NA_real_, P.value=NA_real_, AUC=NA_real_))
		}
    this_dat <- copy(pheno)
	  this_dat[, pgs_grid := this_pgs]
		this_perf <- glm.test(sprintf(mf, "pgs_grid"), "type_2_diabetes", this_dat, skip.ci=TRUE)
		this_perf[coefficient %like% "scale", .(model="grid", paramset=gi, OR, SE=logOR.SE, P.value, AUC)]
	}
	pgs_grid_perf[, grid_param_p := params$p]
	pgs_grid_perf[, grid_param_h2 := params$h2]
	pgs_grid_perf[, grid_param_sparse := params$sparse]

	# Diagnostic plot
	g1 <- ggplot(pgs_grid_perf, aes(x=grid_param_p, y=AUC, color=as.factor(grid_param_h2))) +
		theme_bigstatsr() + 
		geom_point() + 
		geom_line() + 
		scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
		facet_wrap(~ grid_param_sparse, labeller = label_both) +
		labs(y = "AUC in GLM for prevalent T2D", color = "h2") +
		theme(legend.position = "top", panel.spacing = unit(1, "lines"))

	g2 <- ggplot(pgs_grid_perf) + 
		aes(x=grid_param_p, color=as.factor(grid_param_h2), 
				y=OR, ymin = exp(log(OR) - SE), ymax = exp(log(OR) + SE)) +
		theme_bigstatsr() + 
		geom_hline(yintercept=1, linetype=2) +
		geom_errorbar(width=0, alpha=0.5) +
		geom_point() + 
		scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
		facet_wrap(~ grid_param_sparse, labeller = label_both) +
		labs(y = "Odds Ratio (+/- SE) for PGS", color = "h2") +
		theme(legend.position = "top", panel.spacing = unit(1, "lines"))

	g <- plot_grid(g1, g2, nrow = 2)
	ggsave(g, width=7.2, height=7.2, units="in", file=sprintf("%s/ldpred2_grid_performance.png", outdir))

  # Add to model performance
  model_perf <- rbind(model_perf, pgs_grid_perf[, .(model="grid", paramset, AUC, OR, SE, P.value)], fill=TRUE)

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
	ggsave(g, width=20, height=10, units="in", file=sprintf("%s/ldpred2_auto_chain_convergence.png", outdir))

	# Extract PGS weights
	beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)

	# Get PGS levels for all models
	pgs_auto <- foreach(this_chr = 1:22, .combine=`+`) %do% {
		geno <- snp_attach(sprintf("%s/filtered_ukb_chr%s.rds", ramdir, this_chr))
		geno <- snp_fastImputeSimple(geno$genotypes, ncores=nCores) # doesn't work with missing genotypes, so need to impute as median
		big_prodMat(
			X = geno, 
			A.col = beta_auto[gwas_ss[!(fail_qc), which(chr == this_chr)],], 
			ind.col = gwas_ss[!(fail_qc) & chr == this_chr, `_NUM_ID_`],
			ncores = nCores
		)
	}

	# Assess model performance
  registerDoMC(nCores/2) # big_prodMat clobbers existing backend
	pgs_auto_perf <- foreach(gi = seq_along(multi_auto), .combine=rbind) %dopar% {
		this_pgs <- pgs_auto[, gi]
		if (all(is.na(this_pgs))) {
			return(data.table(model="auto", paramset=gi, OR=NA_real_, SE=NA_real_, P.value=NA_real_, AUC=NA_real_))
		}
    this_dat <- copy(pheno)
		this_dat[, pgs_auto := this_pgs]
		this_perf <- glm.test(sprintf(mf, "pgs_auto"), "type_2_diabetes", this_dat, skip.ci=TRUE)
		this_perf[coefficient %like% "scale", .(model="auto", paramset=gi, OR, SE=logOR.SE, P.value, AUC)]
	}

	# Add whole model performance (AUC) and odds ratio for the PGS (discarding OR estimates for covariates)
	model_perf <- rbind(model_perf, pgs_auto_perf, fill=TRUE)

	# Flag bad chains
	sc <- apply(pgs_auto, 2, sd)
	keep <- abs(sc - median(sc)) < 3 * mad(sc)

	# Average to get best auto pgs
	if (sum(keep) > 0) {
		final_auto_pgs <- rowMeans(pgs_auto[, keep])
		final_auto_betas <- rowMeans(beta_auto[, keep])
	}

	# Test performance
  pheno[, pgs_auto := final_auto_pgs]
	final_auto_perf <- glm.test(sprintf(mf, "pgs_auto"), "type_2_diabetes", pheno, skip.ci=TRUE)

	# Add whole model performance (AUC) and odds ratio for the PGS (discarding OR estimates for covariates)
	model_perf <- rbind(model_perf, final_auto_perf[coefficient %like% "scale",.(model="auto", paramset="mean", OR, SE=logOR.SE, P.value, AUC)], fill=TRUE)

	##########################
	# Run lassosum2
	##########################
  cat("Running lassosum2 model\n")
	beta_lassosum2 <- snp_lassosum2(genocorr, gwas_ss[!(fail_qc)], ncores = nCores)
	params2 <- attr(beta_lassosum2, "grid_param")

	# Get PGS levels for all parameters
	pgs_lassosum2 <- foreach(this_chr = 1:22, .combine=`+`) %do% {
		geno <- snp_attach(sprintf("%s/filtered_ukb_chr%s.rds", ramdir, this_chr))
		geno <- snp_fastImputeSimple(geno$genotypes, ncores=nCores) # doesn't work with missing genotypes, so need to impute as median
		big_prodMat(
			X = geno, 
			A.col = beta_lassosum2[gwas_ss[!(fail_qc), which(chr == this_chr)],],
			ind.col = gwas_ss[!(fail_qc) & chr == this_chr, `_NUM_ID_`],
			ncores = nCores
		)
	}

	# Assess model performance for all lassosum parameters
  registerDoMC(nCores/2) # big_prodMat clobbers existing backend
	pgs_lassosum2_perf <- foreach(gi = seq_len(nrow(params2)), .combine=rbind) %dopar% {
		this_pgs <- pgs_lassosum2[, gi]
		if (all(is.na(this_pgs))) {
			return(data.table(model="lassosum2", paramset=gi, OR=NA_real_, SE=NA_real_, P.value=NA_real_, AUC=NA_real_))
		}
    this_dat <- copy(pheno)
		this_dat[, pgs_lassosum2 := this_pgs]
		this_perf <- glm.test(sprintf(mf, "pgs_lassosum2"), "type_2_diabetes", this_dat, skip.ci=TRUE)
		this_perf[coefficient %like% "scale", .(model="lassosum2", paramset=gi, OR, SE=logOR.SE, P.value, AUC)]
	}
	pgs_lassosum2_perf <- cbind(params2, pgs_lassosum2_perf)
	setDT(pgs_lassosum2_perf)
	model_perf <- rbind(model_perf, pgs_lassosum2_perf[, .(model="lassosum2", paramset=.I, AUC, OR, SE, P.value)])

	# Diagnostic plot
	g1 <- ggplot(pgs_lassosum2_perf, aes(x=lambda, y=AUC, color=as.factor(delta))) +
		theme_bigstatsr() + 
		geom_point() + 
		geom_line() + 
		scale_x_log10(breaks = 10^(-5:0)) +
		labs(y = "AUC in GLM for prevalent T2D", color = "delta") +
		theme(legend.position = "top", panel.spacing = unit(1, "lines")) +
		guides(colour = guide_legend(nrow = 1))

	g2 <- ggplot(pgs_lassosum2_perf) + 
		aes(x=lambda, color=as.factor(delta), 
				y=OR, ymin = exp(log(OR) - SE), ymax = exp(log(OR) + SE)) +
		theme_bigstatsr() + 
		geom_hline(yintercept=1, linetype=2) +
		geom_errorbar(width=0, alpha=0.5) +
		geom_point() + 
		scale_x_log10(breaks = 10^(-5:0)) +
		labs(y = "Odds Ratio (+/- SE) for PGS", color = "delta") +
		theme(legend.position = "top", panel.spacing = unit(1, "lines")) +
		guides(colour = guide_legend(nrow = 1))

	g <- plot_grid(g1, g2, nrow = 2)
	ggsave(g, width=7.2, height=7.2, units="in", file=sprintf("%s/ldpred2_lassosum_performance.png", outdir))

  # Write out lasso parameters
  fwrite(pgs_lassosum2_perf[,.(paramset=.I, lambda, delta)], sep="\t", quote=FALSE, file=sprintf("%s/ldpred2_lassosum2_parameters.txt", outdir))

	###############################
	# Write out results
	###############################
	# Model performance
	fwrite(model_perf, sep="\t", quote=FALSE, file=sprintf("%s/model_performance_comparison.txt", outdir))

  # PGS variant weights so we can compute the PGS in new samples and perform model selection
  # Models that failed to be fit (or where betas were all 0) are removed
  pgs_betas <- gwas_ss[!(fail_qc), .(rsid=rsid.ss, chr, pos, effect_allele=a1, other_allele=a0)]

  if ("beta_inf" %in% names(gwas_ss)) { 
    pgs_betas <- cbind(pgs_betas, gwas_ss[!(fail_qc), .(beta_inf)])
  }

  beta_grid <- as.data.table(beta_grid)
  setnames(beta_grid, gsub("V", "beta_grid_", names(beta_grid)))
  beta_grid <- beta_grid[, as.integer(model_perf[model == "grid" & !is.na(AUC), paramset]), with=FALSE]
  pgs_betas <- cbind(pgs_betas, beta_grid)

  beta_auto <- as.data.table(beta_auto)
  setnames(beta_auto, gsub("V", "beta_auto_", names(beta_auto)))
  beta_auto <- beta_auto[, as.integer(model_perf[model == "auto" & paramset != "mean" & !is.na(AUC), paramset]), with=FALSE]
  pgs_betas <- cbind(pgs_betas, beta_auto)

  if (!is.na(model_perf[model == "auto" & paramset == "mean", AUC])) {
    pgs_betas[, beta_auto_final := final_auto_betas]
  }

  beta_lassosum2 <- as.data.table(beta_lassosum2)
  setnames(beta_lassosum2, gsub("V", "beta_lassosum2_", names(beta_lassosum2)))
  beta_lassosum2 <- beta_lassosum2[, as.integer(model_perf[model == "lassosum2" & !is.na(AUC), paramset]), with=FALSE]
  pgs_betas <- cbind(pgs_betas, beta_lassosum2)

	fwrite(pgs_betas, sep="\t", quote=FALSE, compress="gzip", file=sprintf("%s/ldpred2_pgs_varweights.txt.gz", outdir))

	# Predicted PGS levels (for sanity checking)
  # Models that failed to be fit (or where betas were all 0) are removed
	pgs_levels <- pheno[,.(eid)]

  if ("pgs_inf" %in% names(pheno)) {
    pgs_levels <- cbind(pgs_levels, pheno[,.(pgs_inf)])
  }

  pgs_grid <- as.data.table(pgs_grid)
  setnames(pgs_grid, gsub("V", "pgs_grid_", names(pgs_grid)))
  pgs_grid <- pgs_grid[, as.integer(model_perf[model == "grid" & !is.na(AUC), paramset]), with=FALSE]
  pgs_levels <- cbind(pgs_levels, pgs_grid)

  pgs_auto <- as.data.table(pgs_auto)
  setnames(pgs_auto, gsub("V", "pgs_auto_", names(pgs_auto)))
  pgs_auto <- pgs_auto[, as.integer(model_perf[model == "auto" & paramset != "mean" & !is.na(AUC), paramset]), with=FALSE]
  pgs_levels <- cbind(pgs_levels, pgs_auto)

  if (!is.na(model_perf[model == "auto" & paramset == "mean", AUC])) {
    pgs_levels[, pgs_auto_final := final_auto_pgs]
  }

	fwrite(pgs_levels, sep="\t", quote=FALSE, compress="gzip", file=sprintf("%s/traininset_pgs_levels.txt.gz", outdir))

  #######################
  # Cleanup
  #######################
  system(sprintf("rm -rf %s", tmpdir), wait=TRUE)
  system(sprintf("rm -rf %s", ramdir), wait=TRUE)
}

