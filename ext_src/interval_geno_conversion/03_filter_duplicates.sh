#!/bin/bash

out_dir=/rds/project/jmmh2/rds-jmmh2-post_qc_data/interval/imputed/uk10k_1000g_b37/imputed/plink_format/pgen
chr=$SLURM_ARRAY_TASK_ID

if [ -f $out_dir/chr${chr}_duplicates.txt ]; then
  # Exclude duplicates and create new pgen files
  plink2 --pfile $out_dir/impute_${chr}_interval \
         --exclude $out_dir/chr${chr}_duplicates.txt \
         --threads $SLURM_CPUS_ON_NODE \
         --memory $SLURM_MEM_PER_NODE \
         --silent \
         --make-pgen \
         --out $out_dir/impute_dedup_${chr}_interval

  # Remove old pgen files.
  rm $out_dir/impute_${chr}_interval.{pgen,pvar,psam}

  # Remove temporary exclusion file
  rm $out_dir/chr${chr}_duplicates.txt
else
  # No duplicates to remove, just mv the files.
  mv $out_dir/impute_${chr}_interval.pgen $out_dir/impute_dedup_${chr}_interval.pgen
  mv $out_dir/impute_${chr}_interval.pvar $out_dir/impute_dedup_${chr}_interval.pvar
  mv $out_dir/impute_${chr}_interval.psam $out_dir/impute_dedup_${chr}_interval.psam
  echo "There were no duplicate variants to remove." > $out_dir/impute_dedup_${chr}_interval.log
fi

