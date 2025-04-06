library(data.table)
library(foreach)
source("src/functions/flip_strand.R")

# Load candidate variant set
varset <- fread("data/filtered_sumstats/candidate_varset.txt")

# Load variant information extracted from All of Us Research Platform by Carles
#
# ACAF Threshold subset of All of Us:
# - Variants that are frequent in the AoU computed ancestry subpopulations.
# - The cutoff we use is population-specific allele frequency (AF) > 1% OR population-specific allele count (AC) > 100, in any computed ancestry subpopulations
#
# https://support.researchallofus.org/hc/en-us/articles/14929793660948-Smaller-Callsets-for-Analyzing-Short-Read-WGS-SNP-Indel-Data-with-Hail-MT-VCF-and-PLINK
aou <- foreach(chrIdx = 1:22, .combine=rbind) %do% {
  fread(sprintf("data/All_of_Us/candidate_bim_files/candidate_varset_chr%s.bim", chrIdx))
}
setnames(aou, c("chr", "varid", "cm", "pos", "alt", "ref"))
aou[, chr := as.integer(gsub("chr", "", chr))]

# Filter to variants present in All of Us by chr and position (29,303 variants not present in ACAF subset)
varset <- varset[unique(aou[,.(chr, pos)]), on = .(chr, pos_b38=pos), nomatch=0]

# Make sure we can also match by allele to All of Us, orienting effect alleles to the ALT allele in All of Us
varset[, allele_match := FALSE]
varset[aou, on = .(chr, pos_b38=pos, effect_allele=alt, other_allele=ref), c("allele_match", "oriented", "flipped") := .(TRUE, TRUE, FALSE)]
varset[aou, on = .(chr, pos_b38=pos, effect_allele=ref, other_allele=alt), c("allele_match", "oriented", "flipped") := .(TRUE, FALSE, FALSE)]

# Some SNPs may be on the opposite strand in 1000G compared to All of Us, which we detect and fix here (n=1,237 SNPs)
varset[!(allele_match), c("effect_allele", "other_allele", "flipped") := .(flip_strand(effect_allele), flip_strand(other_allele), TRUE)]
varset[aou, on = .(chr, pos_b38=pos, effect_allele=alt, other_allele=ref), c("allele_match", "oriented") := .(TRUE, TRUE)]
varset[aou, on = .(chr, pos_b38=pos, effect_allele=ref, other_allele=alt), c("allele_match", "oriented") := .(TRUE, FALSE)]

# Remove any variants which matched on chromosome and position 
# but not on alleles even after checking for strand mismatch 
# (N=88 SNPs)
varset <- varset[(allele_match)]

# Orient effect alleles to the ALT allele in All of Us - LDpred2 loads this allele as the
# dosage, so we ultimately want to align gwas summary stat effect alleles to this
varset[!(oriented),
  c("effect_allele", "other_allele",  "EAS_EAF_1000G", "AMR_EAF_1000G", "AFR_EAF_1000G", "EUR_EAF_1000G", "SAS_EAF_1000G", "oriented") :=
  .(other_allele, effect_allele, 1-EAS_EAF_1000G, 1-AMR_EAF_1000G, 1-AFR_EAF_1000G, 1-EUR_EAF_1000G, 1-SAS_EAF_1000G, TRUE)]

# Add in SNP ID from All of Us
varset[aou, on = .(chr, pos_b38=pos, effect_allele=alt, other_allele=ref), AoU_varID := varid]

# Reorganise columns
varset <- varset[,.(chr, pos_b36, pos_b37, pos_b38, rsid_HapMap3, rsid_1000G, AoU_varID, effect_allele, other_allele, 
  EAF_1000G, EUR_EAF_1000G, AFR_EAF_1000G, AMR_EAF_1000G, EAS_EAF_1000G, SAS_EAF_1000G, 
  ASW_HapMap3, CEU_HapMap3, CHB_HapMap3, CHD_HapMap3, GIH_HapMap3, JPT_HapMap3, 
  LWK_HapMap3, MEX_HapMap3, MKK_HapMap3, TSI_HapMap3, YRI_HapMap3)]

# Write out
fwrite(varset[order(pos_b38)][order(chr)], sep="\t", quote=FALSE, "data/filtered_sumstats/filtered_oriented_SNPs.txt")

