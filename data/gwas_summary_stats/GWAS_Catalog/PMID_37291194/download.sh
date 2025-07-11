#!/bin/bash

# N.b., while all of these have 'harmonised' folders, the harmonisation was incorrect; it only succesfully
# harmonised a fraction of the SNPs, due to incorrect reporting/assumption of genome build (assumed GRCh38,
# actually GRCh37)
wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90267001-GCST90268000/GCST90267574/
wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90267001-GCST90268000/GCST90267578/
wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90267001-GCST90268000/GCST90267573/
wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90267001-GCST90268000/GCST90267577/
wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90267001-GCST90268000/GCST90267568/
wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90267001-GCST90268000/GCST90267572/
wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90267001-GCST90268000/GCST90267567/
wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90267001-GCST90268000/GCST90267571/

mv ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90267001-GCST90268000/* .
rm -rf ftp.ebi.ac.uk
