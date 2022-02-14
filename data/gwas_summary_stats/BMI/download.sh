#!/bin/bash

# Get the version that has been harmonized to build 38 - the advantage here is that the original sumstats
# are missing chromosome and position information, and I don't have rsIDs for all candidate variants on
# build 36. Variants from HapMap3 (build36) that did not map to build 37 or build 38 have already been
# excluded from the candidate variant set.
wget -P $(dirname $0) ftp.ebi.ac.uk:/pub/databases/gwas/summary_statistics/GCST002001-GCST003000/GCST002783/harmonised/25673413-GCST002783-EFO_0004340.h.tsv.gz
