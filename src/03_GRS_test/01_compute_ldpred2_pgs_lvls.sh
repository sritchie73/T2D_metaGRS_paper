#!/bin/bash

# Make output directory:
mkdir -p output/ldpred2/all_hyperparam_grs_lvls

# Loop through GWAS and compute PGS for all LDpred2 hyperparameters
for gwas_dir in output/ldpred2/train/*; do
  gwas=$(basename $gwas_dir)

  rm -rf output/ldpred2/all_hyperparam_grs_lvls/$gwas

  # Use existing pipeline to compute scores
	./src/PGS_resources/calc_PS_lvls.sh \
    --score-file output/ldpred2/train/$gwas/ldpred2_pgs_varweights.txt.gz \
    --type 's' --score-weight 'm' \
		--genotype-prefix 'data/UKB/genetics/imputed_pgen/ukb_imp_v3_dedup_chr' \
		--keep-ambiguous \
		--work output/ldpred2/all_hyperparam_grs_lvls/$gwas \
		--out output/ldpred2/all_hyperparam_grs_lvls/$gwas \
    --partition icelake \
    --time 8:0:0 \
    --single-out collated_scores

  # SLURM doesn't like it when you make many simultaneous calls to sbatch, so sleep
  # for 10 seconds between sbatch calls
  sleep 10
done
