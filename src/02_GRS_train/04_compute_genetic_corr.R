library(data.table)
library(bigsnpr)

# Set up parallelisation
source("src/functions/par_setup.R")

# Load set of variants 
varset <- fread("data/filtered_sumstats/filtered_oriented_SNPs.txt")

# Genetic correlation is computed for each chromosome separately
for (this_chr in 1:22) { 
  # Load filtered SNP data at this chromosome
  geno <- snp_attach(sprintf("data/ldpred2/filtered_ukb_chr%s.rds", this_chr)) 

  # Get genetic distance from 1000G
  gendist <- snp_asGeneticPos(
    infos.chr = geno$map$chromosome, 
    infos.pos = geno$map$physical.pos,
    dir = "data/1000G/genetic_maps/interpolated_OMNI/"
  )

  # Compute LD
  gencorr <- snp_cor(
    Gna = geno$genotypes, 
    infos.pos = gendist, 
    size = 3 / 1000,
    ncores = nCores # detected in src/functions/par_setup.R
  ) 

  # Save out
  saveRDS(gencorr, file=sprintf("data/ldpred2/filtered_ukb_ldcorr_chr%s.rds", this_chr))
}

