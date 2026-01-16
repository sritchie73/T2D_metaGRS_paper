library(data.table)
library(foreach)
source("src/functions/flip_strand.R")
source("src/functions/eaf_to_maf.R")

# Load candidate variant set
varset <- fread("data/filtered_sumstats/candidate_varset.txt")

# Load variant information extracted from All of Us Research Platform by Carles
#
# ACAF Threshold subset of All of Us:
# - Variants that are frequent in the AoU computed ancestry subpopulations.
# - The cutoff we use is population-specific allele frequency (AF) > 1% OR population-specific allele count (AC) > 100, in any computed ancestry subpopulations
#
# https://support.researchallofus.org/hc/en-us/articles/14929793660948-Smaller-Callsets-for-Analyzing-Short-Read-WGS-SNP-Indel-Data-with-Hail-MT-VCF-and-PLINK
aou <- fread("data/All_of_Us/candidate_varset_AF.txt.gz")
aou[, chr := as.integer(gsub("chr", "", contig))]

# Filter to variants present in All of Us by chr and position (11,601 variants not present in ACAF subset)
varset <- varset[unique(aou[,.(chr, position)]), on = .(chr, pos_b38=position), nomatch=0]

# Make sure we can also match by allele to All of Us, orienting effect alleles to the ALT allele in All of Us
varset[, allele_match := FALSE]
varset[aou, on = .(chr, pos_b38=position, effect_allele=alt_allele, other_allele=ref_allele), c("allele_match", "oriented", "flipped") := .(TRUE, TRUE, FALSE)]
varset[aou, on = .(chr, pos_b38=position, effect_allele=ref_allele, other_allele=alt_allele), c("allele_match", "oriented", "flipped") := .(TRUE, FALSE, FALSE)]

# Some SNPs may be on the opposite strand in the HM3+ variant set compared to All of Us, which we detect and fix here (n=1,790 SNPs)
varset[!(allele_match), c("effect_allele", "other_allele", "flipped") := .(flip_strand(effect_allele), flip_strand(other_allele), TRUE)]
varset[aou, on = .(chr, pos_b38=position, effect_allele=alt_allele, other_allele=ref_allele), c("allele_match", "oriented") := .(TRUE, TRUE)]
varset[aou, on = .(chr, pos_b38=position, effect_allele=ref_allele, other_allele=alt_allele), c("allele_match", "oriented") := .(TRUE, FALSE)]

# Remove any variants which matched on chromosome and position 
# but not on alleles even after checking for strand mismatch 
# (N=658 SNPs, mostly indels at those positions in All of Us instead of SNPs)
varset <- varset[(allele_match)]

# Orient effect alleles to the ALT allele in All of Us - LDpred2 loads this allele as the
# dosage, so we ultimately want to align gwas summary stat effect alleles to this
varset[!(oriented),  c("effect_allele", "other_allele", "oriented") := .(other_allele, effect_allele, TRUE)] # N=2,743 SNPs

# Add in SNP ID from All of Us
varset[aou, on = .(chr, pos_b38=position, effect_allele=alt_allele, other_allele=ref_allele), AoU_varID := vid]

# Add in allele frequencies from ACAF callset
varset[aou, on = .(AoU_varID=vid), 
  c("rsid_AoU", "EUR_EAF_AoU", "AFR_EAF_AoU", "AMR_EAF_AoU", "EAS_EAF_AoU", "SAS_EAF_AoU", "MID_EAF_AoU", "OTH_EAF_AoU") :=
  .(dbsnp_rsid, gvs_eur_af, gvs_afr_af, gvs_amr_af, gvs_eas_af, gvs_sas_af, gvs_mid_af, gvs_oth_af)]

# SNPs with allele counts < 20 have their frequency set to NA for that ancestry in All of Us. To simplify filtering, set these to 0 in our
# varset table
varset[is.na(EUR_EAF_AoU), EUR_EAF_AoU := 0]
varset[is.na(AFR_EAF_AoU), AFR_EAF_AoU := 0]
varset[is.na(AMR_EAF_AoU), AMR_EAF_AoU := 0]
varset[is.na(EAS_EAF_AoU), EAS_EAF_AoU := 0]
varset[is.na(SAS_EAF_AoU), SAS_EAF_AoU := 0]
varset[is.na(MID_EAF_AoU), MID_EAF_AoU := 0]
varset[is.na(OTH_EAF_AoU), OTH_EAF_AoU := 0]

# Filter to variants that have >1% frequency in any All of Us ancestry - this filter is applied later on by the
# per-ancestry LDpred2 quality control when calculating the LD matrix (assuming 10K samples)
varset <- varset[
  maf(EUR_EAF_AoU) >= 0.01 |   # 1,326,176 SNPs pass, 105,761 below 1% frequency
  maf(AFR_EAF_AoU) >= 0.01 |   # 1,177,351 SNPs pass, 254,586 below 1% frequency
  maf(AMR_EAF_AoU) >= 0.01 |   # 1,268,464 SNPs pass, 163,473 below 1% frequency
  maf(EAS_EAF_AoU) >= 0.01 |   # 1,036,658 SNPs pass, 395,279 below 1% frequency
  maf(SAS_EAF_AoU) >= 0.01 |   # 1,217,613 SNPs pass, 214,324 below 1% frequency
  maf(MID_EAF_AoU) >= 0.01 |   # 1,255,463 SNPs pass, 176,474 below 1% frequency 
  maf(OTH_EAF_AoU) >= 0.01     # 1,308,505 SNPs pass, 123,432 below 1% frequency
  # 1,382,696 SNPs pass overall, 49,241 below 1% frequency in all ancestries
]

# No need to handle ambiguous SNPs, excluded by design in HM3+ set
varset[, ambig := effect_allele == flip_strand(other_allele)] # Sanity check, 0 SNPs

# Add in variant ID column used by AoU .bim files
varset[, AoU_varID2 := AoU_varID]
varset[, AoU_varID := paste0("chr", gsub("-", ":", AoU_varID2))]

# Rename some columns for consistency
setnames(varset, c("EUR_EAF_AoU", "AFR_EAF_AoU", "AMR_EAF_AoU", "EAS_EAF_AoU", "SAS_EAF_AoU", "MID_EAF_AoU", "OTH_EAF_AoU"),
                 c("EAF_EUR_AoU", "EAF_AFR_AoU", "EAF_AMR_AoU", "EAF_EAS_AoU", "EAF_SAS_AoU", "EAF_MID_AoU", "EAF_OTH_AoU"))

# Reorganise columns
varset <- varset[,.(chr, pos_b36, pos_b37, pos_b38, rsid_HapMap3, rsid_1000G, rsid_AoU, AoU_varID, AoU_varID2, effect_allele, 
  other_allele, EAF_EUR_AoU, EAF_AFR_AoU, EAF_AMR_AoU, EAF_EAS_AoU, EAF_SAS_AoU, EAF_MID_AoU, EAF_OTH_AoU, EAF_UKBB,
  EAF_1000G, EAF_AFR_1000G, EAF_AMR_1000G, EAF_EAS_1000G, EAF_EUR_1000G, EAF_SAS_1000G, ld, block_id)]

# Write out
fwrite(varset[order(pos_b38)][order(chr)], sep="\t", quote=FALSE, "data/filtered_sumstats/filtered_oriented_SNPs.txt")

