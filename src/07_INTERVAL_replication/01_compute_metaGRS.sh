#!/bin/bash

mkdir -p data/INTERVAL/T2D_metaGRS

# Compute PGS
./src/PGS_resources/calc_PS_lvls.sh \
	--score-file output/metaGRS/hyperparam_selection/T2D_metaGRS.txt.gz \
	--genotype-prefix 'data/INTERVAL/genetics/imputed_pgen/impute_dedup_' \
  --genotype-suffix '_interval' \
	--keep-ambiguous \
  --out data/INTERVAL/T2D_metaGRS/ \
	--work data/INTERVAL/T2D_metaGRS/ \
	--time 8:0:0 \
	--single-out T2D_metaGRS

