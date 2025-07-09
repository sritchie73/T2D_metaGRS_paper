#!/bin/bash

wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007140/
wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007141/
wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007142/
wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007143/

mv ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/* .
rm -rf ftp.ebi.ac.uk

