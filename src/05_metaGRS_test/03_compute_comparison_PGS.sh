#!/bin/bash

# Build file listing PGS to compute
uuid=$(uuidgen)
mkdir -p tmp/$uuid
touch tmp/$uuid/paths_to_scores.txt
for ff in data/T2D_PGS/UKB_valid/*/*.txt.gz; do
  echo $ff >> tmp/$uuid/paths_to_scores.txt
done

# Compute PGS
./src/PGS_resources/calc_PS_lvls.sh \
	--score-file tmp/$uuid/paths_to_scores.txt \
	--type 'l' \
	--genotype-prefix 'data/UKB/genetics/imputed_pgen/ukb_imp_v3_dedup_chr' \
	--keep-ambiguous \
  --out data/UKB/T2D_PGS/ \
	--work data/UKB/T2D_PGS/ \
	--time 8:0:0 \
	--single-out collated_scores
