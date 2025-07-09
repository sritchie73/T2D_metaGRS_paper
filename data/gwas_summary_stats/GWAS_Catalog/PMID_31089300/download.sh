#!/bin/bash

wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009403/
wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009402/
wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009399/
wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009405/
wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009406/
wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009407/

mv ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/* .
rm -rf ftp.ebi.ac.uk

