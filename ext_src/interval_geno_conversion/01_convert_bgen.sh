#!/bin/bash

ref_dir=/rds/project/jmmh2/rds-jmmh2-post_qc_data/interval/imputed/uk10k_1000g_b37/imputed
out_dir=$ref_dir/plink_format/pgen

chr=$SLURM_ARRAY_TASK_ID

# Chromosomes 1-22
plink2 --bgen $ref_dir/impute_${chr}_interval.bgen \
       --sample $ref_dir/interval.samples \
       --threads $SLURM_CPUS_ON_NODE \
       --memory $SLURM_MEM_PER_NODE \
       --silent \
       --make-pgen \
       --out $out_dir/impute_${chr}_interval

