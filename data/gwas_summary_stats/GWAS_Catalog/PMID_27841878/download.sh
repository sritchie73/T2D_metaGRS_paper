#!/bin/bash

wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007098/
wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007095/
wget -m ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007097/

mv ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/* .
rm -rf ftp.ebi.ac.uk

