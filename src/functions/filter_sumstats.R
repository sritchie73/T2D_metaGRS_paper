library(data.table)
library(foreach)
source("src/functions/flip_strand.R")

# Note n_eff for case control studies is computed as 4/(1/cases + 1/controls) and sample size 
# for continuous traits. Case/control numbers are obtained directly from the study/GWAS Catalog,
# then also weighted by total sample size (if provided and varying by SNP in the summary stats)
# or total number of contributing studies (if provided for each SNP in the absence of per-SNP 
# sample size). 
n_eff <- function(cases, controls, per_snp_N=1, total_N=1) {
  cases <- cases / total_N * per_snp_N
  controls <- controls / total_N * per_snp_N
  4 / (1/cases + 1/controls)
}

# Some GWAS are missing the standard error, for additive effects in linear models we can 
# back-calculate the T-statistic using the sample size and P-value
lm_se <- function(beta, neg_log10_p, samples, n_covar=0) {
  df <- samples - 2 - n_covar # in practice for any sufficiently large N for GWAS, subtracting right number of covariates is irrelevant
  t_stat <- qt(10^(-neg_log10_p)/2, df = df)
  abs(beta)/abs(t_stat)
}

# Or from 95% Confidence intervals if those are reported
se_from_ci <- function(beta, L95, U95) {
  Zconst <- qnorm(p=0.05/2, lower.tail=F)
  se_from_U95 = (U95 - beta)/Zconst # Rearrangement of 95% CI = beta +/- 1.96 * se
  se_from_L95 = (beta - L95)/Zconst
  (se_from_L95 + se_from_U95)/2 # as either limit can give slightly different se
}

# Correction of BOLT-LMM beta and se for case %
# https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html#x1-470008
bolt_lmm_fix <- function(x, cases, controls) {
  u <- cases / (cases + controls)
  x  / (u * (1-u))
}

# Stop bad fread when running out of space in /tmp
fread <- function(...) {
  data.table::fread(..., tmpdir="tmp")
}

# Load varset table, if not already in global environment
if (!exists("varset")) {
  varset <- fread("data/filtered_sumstats/filtered_oriented_SNPs.txt", na.strings=c("NA", ""))
}

# Filter and format sumstats orienting to AoU alleles 
filter_sumstats <- function(gwas_ss, type, total_samples, total_cases, total_controls) {
  # Get chromosome and position if missing
  if ("rsid" %in% names(gwas_ss)) {
    gwas_ss[varset, on = .(rsid=rsid_1000G), c("chr", "pos_b38") := .(i.chr, i.pos_b38)]
    gwas_ss <- gwas_ss[!is.na(chr) & !is.na(pos_b38)]
    gwas_ss[, rsid := NULL]
  }

  # Filter to autosomal chromosomes
  if (typeof(gwas_ss$chr) != "integer") {
    gwas_ss <- gwas_ss[chr %in% 1:22]
    gwas_ss[, chr := as.integer(chr)]
  }
  gwas_ss <- gwas_ss[chr >= 1 & chr <= 22]

  # Get positions on build 38 if necessary
  if ("pos_b37" %in% names(gwas_ss)) {
    gwas_ss[varset, on = .(chr, pos_b37), pos_b38 := i.pos_b38]
    gwas_ss <- gwas_ss[!is.na(pos_b38)]
  } else if ("pos_b36" %in% names(gwas_ss)) {
    gwas_ss[varset, on = .(chr, pos_b36), pos_b38 := i.pos_b38]
    gwas_ss <- gwas_ss[!is.na(pos_b38)]
  }

  # Add in EAF column if missing to simplify downstream code
  if (!("EAF" %in% names(gwas_ss))) {
    gwas_ss[, EAF := NA]
  }

  # Require finite and non-missing weights and standard errors
  gwas_ss <- gwas_ss[is.finite(beta)]
  gwas_ss <- gwas_ss[is.finite(beta_se)]

  # Filter to variants in the candidate variant set by chromosome and position
  gwas_ss <- gwas_ss[varset[,.(chr, pos_b38)], on = .(chr, pos_b38), nomatch=0]

  # Add in All of Us variant identifier 
  gwas_ss[varset, on = .(chr, pos_b38), AoU_varID := i.AoU_varID]

  # Make sure we can also match by allele
  gwas_ss[, allele_match := FALSE]
  gwas_ss[varset, on = .(chr, pos_b38, EA=effect_allele, OA=other_allele), c("allele_match", "oriented", "flipped") := .(TRUE, TRUE, FALSE)]
  gwas_ss[varset, on = .(chr, pos_b38, EA=other_allele, OA=effect_allele), c("allele_match", "oriented", "flipped") := .(TRUE, FALSE, FALSE)]

  # Some SNPs may be on the opposite strand in the GWAS, which we detect and fix here:
  gwas_ss[!(allele_match), c("EA", "OA", "flipped") := .(flip_strand(EA), flip_strand(OA), TRUE)]
  gwas_ss[varset, on = .(chr, pos_b38, EA=effect_allele, OA=other_allele), c("allele_match", "oriented") := .(TRUE, TRUE)]
  gwas_ss[varset, on = .(chr, pos_b38, EA=other_allele, OA=effect_allele), c("allele_match", "oriented") := .(TRUE, FALSE)]

  # Remove any variants which matched on chromosome and position (or rsid) but not on alleles even after checking for strand mismatch
  gwas_ss <- gwas_ss[!gwas_ss[!(allele_match)], on = .(chr, pos_b38)] # anti-join in case multi-allelic variants with match for some alleles

  # Drop variants that were multi-allelic in the gwas
  mult <- gwas_ss[,.N,by=.(chr, pos_b38)][N > 1]
  gwas_ss <- gwas_ss[!mult, on =.(chr, pos_b38)]

  # For strand ambiguous alleles (A/T or G/C SNPs) we will assume same strand orientation as (majority) rest of sumstats 
  pct_flipped <- gwas_ss[EA != flip_strand(OA), sum(flipped)/.N]
  if (pct_flipped < 0.5) {
    gwas_ss[(flipped) & EA == flip_strand(OA), c("EA", "OA", "EAF", "beta", "oriented", "flipped") := .(OA, EA, 1 - EAF, -beta, !oriented, FALSE)]
  } else {
    gwas_ss[!(flipped) & EA == flip_strand(OA), c("EA", "OA", "EAF", "beta", "oriented", "flipped") := .(OA, EA, 1 - EAF, -beta, !oriented, TRUE)]
  }

  # Now fix orientation of effect alleles to match the candidate varset
  gwas_ss[!(oriented), c("EA", "OA", "EAF", "beta", "oriented") := .(OA, EA, 1-EAF, -beta, TRUE)]

  # If the GWAS has EAF information, we can also cross-check strand ambiguous SNPs with EAF in All of Us
  bad <- gwas_ss[EA == flip_strand(OA) & EAF > 0.42 & EAF < 0.58]
  gwas_ss <- gwas_ss[!bad, on = .(chr, pos_b38)]

  # Check cases where EAF does not agree with AoU EAF for strand ambiguous alleles
  # In these cases we need to (1) flip the strand of the alleles, and (2) reverse the orientation
  # of effect/other alleles. This (1) basically means we don't have to change the 'EA' and 'OA'
  # columns because the two operations cancel out, but (2) we have to flip the 'beta' and 'EAF'
  gwas_ss[varset, on=.(chr, pos_b38), EAF_AoU := EUR_EAF_AoU] # Note, variants filtered so all ancestries have same MAF direction
  gwas_ss[EA == flip_strand(OA) & ( (EAF < 0.5 & EAF_AoU > 0.5) | (EAF > 0.5 & EAF_AoU < 0.5) ),
          c("beta", "EAF", "flipped") := .(-beta, 1-EAF, !flipped)]

  # Compute effective sample size if not provided
  if (!missing(type) && type == "continuous") {
    if ("samples" %in% names(gwas_ss)) {
      gwas_ss[, n_eff := samples]
    } else {
      gwas_ss[, n_eff := total_samples]
    }
  } else if (!missing(type) && type == "case/control") {
    if (!("cases" %in% names(gwas_ss))) {
      gwas_ss[, cases := total_cases]
    }
    if (!("controls" %in% names(gwas_ss))) {
      gwas_ss[, controls := total_controls]
    }
    if ("samples" %in% names(gwas_ss)) {
      if (missing(total_samples)) {
        total_samples <- gwas_ss[, max(samples)]
      }
      gwas_ss[, n_eff := n_eff(cases, controls, samples, total_samples)]
    } else {
      gwas_ss[, n_eff := n_eff(cases, controls)]
    }
  }

  # Extract columns of interest and return
  gwas_ss[, .(AoU_varID, chr, pos=pos_b38, a1=EA, a0=OA, a1freq=EAF, beta, beta_se, neg_log10_p, n_eff)]
}

