# Suggested SLURM commands to run for each of the scripts in src/05_metaGRS_test/,
# documenting expected resource usage and requirements

echo "This should not be run as a script"
exit 1

# First set of scripts either dispatch to slurm themselves, or run with trivial resources

# Compute levels of all ~1300 candidate metaGRSs in UKB
./src/05_metaGRS_test/01_compute_candidate_metaGRS_lvls.sh

# Download from PGS Catalog all T2D PGS that can be tested in UKB
./src/05_metaGRS_test/02_download_comparison_PGS.sh

# Compute those external T2D PGS in UKB
./src/05_metaGRS_test/03_compute_comparison_PGS.sh

