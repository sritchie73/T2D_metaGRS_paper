#!/bin/bash

wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009151/
wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009150/
wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009145/
wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009152/

mv ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/* .
rm -rf ftp.ebi.ac.uk

