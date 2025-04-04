# Suggested SLURM commands to run for each of the scripts in src/03_GRS_test/,
# documenting expected resource usage and requirements

echo "This should not be run as a script"
exit 1

# Dispatch script to submit jobs to compute GRS levels in all UK Biobank for
# all LDpred2 hyperparameters for each of the GWAS. Takes a few minutes to 
# run due to sleep() command so that the repeated calls to sbatch (1 per GWAS)
# don't crash slurm.
./src/03_GRS_test/01_compute_ldpred2_pgs_lvls.sh

# Test the performance of all the LDpred2 models.
# Array indices should be updated if the number of GWAS change
sbatch --array=1-44 --mem 40000 --time 12:0:0 --partition skylake,cclake \
  --wrap "Rscript --vanilla src/03_GRS_test/02_ldpred2_test.R" 

# Extract best hyperparameter for each model/cohort/timepoint
Rscript src/03_GRS_test/03_hyperparam_select.R

