#!/bin/bash

# Need to run separately for each ancestry and gwas
ancestry="EUR"
gwas="T2D_MVP_Multi"

# Make output directory:
mkdir -p output/ldpred2/all_hyperparam_grs_lvls/$ancestry/$gwas

# The below uses my pilot script predating the PGS Catalog Calculator to
# simultaneously calculate all possible PRS arising from the LDpred2 training
# for each gwas and ancestry. This pipeline is not readily deployable to 
# other compute environments - you will need to use the PGS Catalog Calculator here
./src/PGS_resources/calc_PS_lvls.sh \
	--score-file output/ldpred2/train/$ancestry/$gwas/ldpred2_pgs_varweights.txt.gz \
	--type 's' --score-weight 'm' \
	--genotype-prefix 'data/UKB/genetics/imputed_pgen/ukb_imp_v3_dedup_chr' \
	--keep-ambiguous \
	--work output/ldpred2/all_hyperparam_grs_lvls/$ancestry/$gwas \
	--out output/ldpred2/all_hyperparam_grs_lvls/$ancestry/$gwas \
	--partition icelake \
	--time 8:0:0 \
	--single-out collated_scores

