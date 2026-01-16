The scripts in this folder are designed to do the hyperparameter selection for the
LDpred2 training done on each GWAS in each ancestry, i.e. to select the optimal PRS
in terms of T2D prediction to carry forward for the metaPRS training

There are three scripts in this folder:

 (1) 01_compute_ldpred2_pgs_lvls.sh: for each GWAS and ancestry, compute all possible LDpred2
                                     PRS in the metaPRS training samples for the respective 
							   ancestry. You will need to adapt this to use the PGS 
                                     Catalog Calculator rather than my CSD3-specific pilot
                                     script.

 (2) 02_ldpred2_test.R: for each GWAS and ancestry, compute the associations between all
                        possible LDpred2 PRS and T2D

 (3) 03_hyperparam_select.R: for each ancestry, select the optimal LDpred2 model for each
                             GWAS aggregating them into a single file per ancestry to take
                             forward for metaPRS training.

In terms of resource requirements:

 (1) 01_compute_ldpred2_pgs_lvls.sh: unsure, have not tried the PGS Catalog Calculator, but assuming
                                     it can simultaneously calculate multiple PRSs it should be 
                                     relatively quick for each GWAS and ancestry if filtering the
                                     genotypes to the metaPRS training samples

 (2) 02_ldpred2_test.R: should be relatively quick and low memory - just some basic logistic regressions
                        and plots for each GWAS and ancestry

 (3) 03_hyperparam_select.R: should be relatively quick and low memory - just selecting, filtering, and
                             aggregating data already computed by previous scripts

Outputs that would be useful to have back from these scripts:

 (1) PGS Catalog Calculator log files
 (2) Association test plots comparing LDpred2 models and hyperparameters for each GWAS and ancestry
 (3) Association test summary statistics 'all_model_performance.txt' for each GWAS and ancestry
 (4) Summarised information on optimal PRSs for each ancestry in 'optimal_prs_info.txt'


