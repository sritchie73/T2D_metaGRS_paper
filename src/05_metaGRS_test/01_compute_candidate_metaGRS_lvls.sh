#!/bin/bash

# Make output directory:
mkdir -p output/metaGRS/all_candidate_metaGRS_lvls

for sfile in output/metaGRS/train/*.txt.gz; do
  out_dir=$(echo $(basename $sfile) | sed 's/.txt.gz//')
  out_dir=output/metaGRS/all_candidate_metaGRS_lvls/$out_dir
  rm -rf $out_dir

  # Use existing pipeline to compute scores
	./src/PGS_resources/calc_PS_lvls.sh \
    --score-file $sfile \
    --type 's' --score-weight 'm' \
		--genotype-prefix 'data/ukb/genetics/imputed_pgen/ukb_imp_v3_dedup_chr' \
		--keep-ambiguous \
		--work $out_dir \
		--out $out_dir \
    --time 8:0:0 \
    --single-out collated_scores

  # SLURM doesn't like it when you make many simultaneous calls to sbatch, so sleep
  # for 10 seconds between sbatch calls
  sleep 10
done
