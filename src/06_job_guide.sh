# Suggested SLURM commands to run for each of the scripts in src/06_INTERVAL_replication/,
# documenting expected resource usage and requirements

echo "This should not be run as a script"
exit 1

# First set of scripts either dispatch to slurm themselves, or run with trivial resources


# Download from PGS Catalog all T2D PGS that can be tested in INTERVAL (that haven't already
# been downloaded previously for calculation in UKB)
./src/06_INTERVAL_replication/02_download_comparison_PGS.sh

# Compute those external T2D PGS in UKB
./src/06_INTERVAL_replication/03_compute_comparison_PGS.sh

