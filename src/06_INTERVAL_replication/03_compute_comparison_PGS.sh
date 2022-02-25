#!/bin/bash

mkdir -p data/INTERVAL/T2D_PGS 

# Build file listing PGS to compute
uuid=$(uuidgen)
mkdir -p tmp/$uuid
touch tmp/$uuid/paths_to_scores.txt
for ff in data/T2D_PGS/*/*/*.txt.gz; do
  echo $ff >> tmp/$uuid/paths_to_scores.txt
done

# Compute PGS
./src/PGS_resources/calc_PS_lvls.sh \
	--score-file tmp/$uuid/paths_to_scores.txt \
	--type 'l' \
	--genotype-prefix 'data/INTERVAL/genetics/imputed_pgen/impute_dedup_' \
  --genotype-suffix '_interval' \
	--keep-ambiguous \
  --out data/INTERVAL/T2D_PGS/ \
	--work data/INTERVAL/T2D_PGS/ \
	--time 8:0:0 \
	--single-out collated_scores

