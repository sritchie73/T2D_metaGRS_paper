The script 01_run_ldpred.R in this folder runs the LDpred2 procedure. This needs to be run in each 
ancestry and for each of the 473 input GWASs separately.

Before running thise scripts you will need to have:

 (1) Downloaded the 473 harmonized and filtered GWASs from data/filtered_sumstats/filtered_gwas/
 (2) Downloaded the table of curated GWAS information data/gwas_summary_stats/curated_gwas_table.txt

In terms of resource requirements, the batch jobs I used for UKB were allocated 12 hours and 160 GB
of RAM, but I'm not sure if they ended up using this much time/memory in practice, and I've stripped
out a lot of extra code I had for computing the candidate PRSs in the training samples and testing
for T2D association as I did not end up using that step to inform model training.

The second script, 02_collate_qc.R, aggregates information about the LDpred2 SNP QC failure rate
across all GWASs, and needs to be run in each ancestry separately.

Outputs that would be useful to have back from these scripts:

 (1) auto_model_chain_convergence.png for each GWAS and ancestry
 (2) aggregated_ldpred2_gwasqc_fail_rate.txt for each ancestry

