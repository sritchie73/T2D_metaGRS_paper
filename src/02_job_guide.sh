# Suggested SLURM commands to run for each of the scripts in src/02_GRS_train/, 
# documenting expected resource usage and requirements

echo "This should not be run as a script"
exit 1

# Obtains the candidate set of 1.04 million variants. Briefly, it:
#
# (1) Obtains the set of bi-allelic SNPs in HapMap3 and exome GWAS
# (2) Filters to SNPs present in 1000 Genomes phase 3
# (4) Orients alleles to the strand present in the UKB data
# (5) Sets the effect allele as the first column in the UKB bim file for
#     consistency with LDpred2.
#
# This process is also applied to variants present in the Exome GWAS 
# considered in this study to increase coverage of exome variants 
# (The filtering to HapMap3 variants is primarily to produce a small 
#  enough set of variants that LDpred2 can work with).
#
# The script also obtains genomic locations for all candidate SNPs on
# genome builds GRCh36/hg18, GRCh37/hg19, and GRCh38/hg38.
sbatch --mem 50000 --time 2:0:0 --wrap "Rscript --vanilla src/02_GRS_train/01_candidate_varset.R"

# Filters each set of downloaded GWAS summary statistics to the 
# candidate SNP set generated by the previous script. It:
#
# (1) Aligns the effect allele and strand to the one in the candidate SNP set
# (2) Converts positions to GRCh37/hg19 if needed
# (3) Estimate the effective sample size if not provided already as 
#     4/(1/case + 1/controls) for GWAS of binary traits, or as the total sample 
#     size otherwise.
#
# Filter GWAS summary stats generated by this script are ouput to data/filtered_sumstats/
sbatch --mem 20000 --time 2:0:0 --wrap "Rscript --vanilla src/02_GRS_train/02_filter_sumstats.R"

# Extracts the genotype data for the candidate variant set above and the training samples
# above from the imputed hard call genotype data and converts them to the format required
# by LDpred2. Data are output to data/ldpred2/
sbatch -c 22 --time 2:0:0 --wrap "Rscript --vanilla ./src/02_GRS_train/03_extract_geno.R"

# Compute SNP-SNP correlation matrices on each chromosome to be used by LDpred2.
sbatch --mem 80000 --time 2:0:0 --wrap "Rscript --vanilla ./src/02_GRS_train/04_compute_genetic_corr.R"

# Runs the LDpred2 hyperparameter tuning for each GWAS in the training participants.
#
# Output is saved to output/ldpred2/train/
#
# Must be run on a slurm partition whose nodes have a /ramdisks/ filesystem partition
# for storing the file-backed shared memory matrix storing the genetic correlations.
#
# Array indices should be updated if the number of GWAS change
sbatch --array=1-44 --mem 160000 --time 12:0:0 --partition skylake,skylake-himem \
  --wrap "Rscript --vanilla ./src/02_GRS_train/05_run_ldpred.R"

# note you may have to rerun a few of the above - they can take longer depending
# on how many array jobs end up being allocated per node



