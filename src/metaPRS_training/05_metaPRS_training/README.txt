The scripts in this folder train the ancestry-specific metaPRSs and combine
them into a single multi-ancestry metaPRS

There are two scripts in this folder:

 (1) 01_metaPRS_training.R: trains a metaPRS in each ancestry

 (2) 02_metaPRS_aggregation.R: combines the ancestry-specific metaPRS into the
                               final multi-ancestry metaPRS

Outputs that would be useful to have back from these scripts:

 (1) For each ancestry, the enfit.png plot summarising the elasticnet model fit

 (2) For each ancestry, the table of all cross-validation AUCs 'enfit_model_fit.txt'

 (3) For each ancestry, the table of the cross-validation AUC for the best model  
     fit 'enfit_best_model_fit.txt'
 
 (4) For each ancestry, the set of log odds ratios for all PRS contributing to the
     ancestry-specific metaPRS 'best_enfit_coefs.txt'

 (5) For each ancestry, the ancestry-specific metaPRS variant weights e.g. 
     'T2D_EUR_metaPRS.txt.gz'

 (6) The final multi-ancestry metaPRS 'T2D_multiancestry_metaPRS.txt.gz'

