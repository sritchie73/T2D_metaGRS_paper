#!/bin/bash

wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008124
wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008125
wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008126
wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008127

mv ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/* .
rm -rf ftp.ebi.ac.uk

