The goal of the scripts in this folder are to extract the necessary genotype and correlation data that
LDpred2 will use in its training procedure. 

Before running these scripts you will need to have:

 (1) Downloaded data/filtered_sumstats/filtered_oriented_SNPs.txt to the All of Us TRE
 (2) Randomly assigned All of Us participants into three groups (LDpred2 training, metaPRS training, testing)
     as previously discussed.
 (3) Downloaded the 1000 Genomes genetic maps from data/1000G/genetic_maps/interpolated_OMNI/*

Reminder on the design of the random assignment of All of Us participants into the three groups:

 - Needs to be done in each ancestry separately
 - Random sampling should aim to make sure case/control numbers are representative (i.e. within an ancestry,
   all three groups should have the same % of cases and % of controls)
 - For EUR, AMR, and AFR ancestries, the LDpred2 training and metaPRS training groups should have 10k samples
   each, with the remaining samples allocated to the testing group
 - For EAS, SAS, and OTH ancestries we will do a 50/50 split between LDpred2 training and metaPRS training,
   with 0 samples withheld for testing due to smaller sample sizes.
 - The MID ancestry is excluded from analyses due to small sample size.

There are two scripts in this folder that need to be run in each ancestry separately:

 (1) 01_extract_geno.R: extracts the subset of the genotype data (LDpred2 training samples + filtered SNPs) in
                        the format LDpred2 requires to run.

 (2) 02_compute_genetic_corr.R: computes and saves the pairwise correlations between the filtered SNPs in the 
                                LDpred2 training samples, these are used by LDpred2 later for PRS training.

In terms of expected resource requirements:

 (1) 01_extract_geno.R: ran within 2 hours on UKB (11k samples) when parallelised across 22 cores on CSD3.
 (2) 02_compute_genetic_corr.R: ran within 2 hours but required 80 GB of RAM, with the number of cores set
                                in the script based on the node hardware (e.g. 12-15 cores)
 
Outputs that would be useful to have back from these scripts:

 (1) plink log files for each ancestry

