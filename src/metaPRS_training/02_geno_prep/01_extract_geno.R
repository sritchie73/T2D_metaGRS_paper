library(bigsnpr)
library(foreach)
library(data.table)

# What ancestry are we currently extracting the genotype data for? - need to run script once per ancestry (or modify to loop accordingly)
ancestry <- "EUR"  

# Make output directory
system(sprintf("mkdir -p data/ldpred2/%s", ancestry), wait=TRUE) # Hopefully AoU uses a bash environment?

# Get list of samples we will be extracting the genotype data for
# Key point here is we need to write out a file with FID and IID matching whats in the genotype data
# for UKB, this was just the eid column repeated twice.
pheno <- fread("data/ldpred2/collated_curated_data.txt")
train <- pheno[(ldpred2_samples)]
fwrite(train[visit_index == 0,.(eid, eid)], col.names=FALSE, quote=FALSE, sep=" ", file=sprintf("data/ldpred2/%s/training_samples.txt", ancestry))

# Determine the MAF cutoff for this ancestry (1/sqrt(N))
maf_cutoff <- 1/sqrt(train[,.N])

# Get list of variants to extract - these need to match what's in the genotype data
varset <- fread("data/filtered_sumstats/filtered_oriented_SNPs.txt")
fwrite(varset[,.(AoU_varID)], col.names=FALSE, quote=FALSE, sep=" ", file=sprintf("data/ldpred2/%s/filtered_SNPs.txt", ancestry))

# Use plink to create a new set of .bed files containing just the subsetted data
# Here, I've parallelised over the 22 chromosomes, but you could instead parallise plink itself instead
# by controlling the --threads argument
foreach(this_chr = 1:22) %dopar% {  # Requires a parallel backend to be registered already, e.g. via doMC::registerDoMC()
  cmd <- "plink2 --bfile" 
  cmd <- paste(cmd, sprintf("data/genotype_chr%s", this_chr)) # location to full genotype data
  cmd <- paste(cmd, sprintf("--keep data/ldpred2/%s/training_samples.txt", ancestry))
  cmd <- paste(cmd, "--extract data/ldpred2/filtered_SNPs.txt")
  cmd <- paste(cmd, "--maf", maf_cutoff)
  cmd <- paste(cmd, sprintf("--make-bed --out  data/ldpred2/%s/filtered_chr%s", ancestry, this_chr))
  cmd <- paste(cmd, "--threads 1 --memory 6000")
  system(cmd, wait=TRUE)
  return(NULL)
}

# Convert to LDpred2 bigsnpr backing file format
for (this_chr in 1:22) {
  invisible(snp_readBed(sprintf("data/ldpred2/%s/filtered_chr%s.bed", ancestry, this_chr)))
} 

# Remove files we no longer need
system(sprintf("rm data/ldpred2/%s/training_samples.txt data/ldpred2/%s/filtered_SNPs.txt", ancestry, ancestry), wait=TRUE)
system(sprinf("rm data/ldpred2/%s/*.{bed,bim,fam}", ancestry), wait=TRUE)

