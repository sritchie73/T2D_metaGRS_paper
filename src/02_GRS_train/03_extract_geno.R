library(bigsnpr)
library(foreach)
library(data.table)

# Set up parallelism
source("src/functions/par_setup.R")

# Make output directory
system("mkdir -p data/ldpred2", wait=TRUE)

# Get list of training samples to keep
pheno <- fread("data/ukb/collated_curated_data.txt")
train <- pheno[(ldpred2_samples)]
fwrite(train[visit_index == 0,.(eid, eid)], col.names=FALSE, quote=FALSE, sep=" ", file="data/ldpred2/training_samples.txt")

# Get list of variants to keep
varset <- fread("data/filtered_sumstats/filtered_oriented_SNPs.txt")

# Create unique IDs per variant - we'll also do this in UKB.
# This means that we can extract the specific alleles of interest in cases
# where SNPs are multi-allelic in UKB (e.g. due to rare variants not in 1000G).
varset[, rn := .I]
varset[,snp_id := sprintf("%s:%s:%s", chr, pos_b37, paste(sort(c(effect_allele, other_allele)), collapse=":")), by=.(rn)]
varset[, rn := NULL]
fwrite(varset[,.(snp_id)], col.names=FALSE, quote=FALSE, file="data/ldpred2/filtered_SNPs.txt")

# Now do the same for UKB, creating new bim files with recoded IDs and symlinking to the bed/fam files
for (this_chr in 1:22) { 
  system(sprintf("ln -s $(realpath data/ukb/genetics/imputed_bed/ukb_imp_v3_dedup_chr%s.bed) data/ldpred2/full_ukb_chr%s.bed", this_chr, this_chr), wait=TRUE)
  system(sprintf("ln -s $(realpath data/ukb/genetics/imputed_bed/ukb_imp_v3_dedup_chr%s.fam) data/ldpred2/full_ukb_chr%s.fam", this_chr, this_chr), wait=TRUE)
  bim <- fread(sprintf("data/ukb/genetics/imputed_bed/ukb_imp_v3_dedup_chr%s.bim", this_chr))
  setnames(bim, c("chr", "rsid", "cm", "pos", "A1", "A2"))
  bim[, rn := .I]
  bim[varset, on = .(chr, pos=pos_b37), rsid := sprintf("%s:%s:%s", chr, pos, paste(sort(c(A1, A2)), collapse=":")), by=.(rn)]
  bim[, rn := NULL]
  fwrite(bim, col.names=FALSE, quote=FALSE, sep=" ", file=sprintf("data/ldpred2/full_ukb_chr%s.bim", this_chr))
}

# Extract the subset of genotype data we want
foreach(this_chr = 1:22) %dopar% {
  cmd <- "plink2 --bfile"
  cmd <- paste(cmd, sprintf("data/ldpred2/full_ukb_chr%s", this_chr))
  cmd <- paste(cmd, "--keep data/ldpred2/training_samples.txt")
  cmd <- paste(cmd, "--extract data/ldpred2/filtered_SNPs.txt")
  cmd <- paste(cmd, sprintf("--make-bed --out  data/ldpred2/filtered_ukb_chr%s", this_chr))
  cmd <- paste(cmd, "--threads 1 --memory 6000")
  system(cmd, wait=TRUE)
  return(NULL)
}

# Clean up files
system("rm data/ldpred2/full_*", wait=TRUE)
system("rm data/ldpred2/*.log", wait=TRUE)
system("rm data/ldpred2/training_samples.txt data/ldpred2/filtered_SNPs.txt", wait=TRUE)

# Re-assign rsIDs to filtered bim files
for (this_chr in 1:22) {
  bim <- fread(sprintf("data/ldpred2/filtered_ukb_chr%s.bim", this_chr))
  setnames(bim, c("chr", "rsid", "cm", "pos", "A1", "A2"))
  bim[varset, on = .(rsid=snp_id), rsid := rsid_1000G]
  fwrite(bim, sep="\t", quote=FALSE, file=sprintf("data/ldpred2/filtered_ukb_chr%s.bim", this_chr))
}

# Convert to LDpred2 bigsnpr backing file format
for (this_chr in 1:22) {
  invisible(snp_readBed(sprintf("data/ldpred2/filtered_ukb_chr%s.bed", this_chr)))
} 

# Remove bed/bim/fam files
system("rm data/ldpred2/*.{bed,bim,fam}", wait=TRUE)
