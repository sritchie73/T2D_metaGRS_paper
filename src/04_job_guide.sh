# Suggested SLURM commands to run for each of the scripts in src/03_metaGRS_train/,
# documenting expected resource usage and requirements

echo "This should not be run as a script"
exit 1

# Meta-GRS training - array job here covers 30x case-control definitions multiplied by 3 pre-filters
sbatch --array=1-90 --time 3:0:0 --mem 80000 \
  --wrap "Rscript --vanilla src/04_metaGRS_train/01_metaGRS_train.R"

# Collate scores into 3 files for PGS level calculation
sbatch --mem 80000 --time 3:0:0 \
  --wrap "Rscript --vanilla src/04_metaGRS_train/02_collate_scores.R"
