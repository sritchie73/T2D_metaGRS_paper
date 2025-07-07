#!/bin/bash

wget -P $(dirname $0) ftp.ebi.ac.uk:/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007847/BBJ_BetaBased1.MAF_001.AtLeast2studies.AllChr.txt.gz
wget -P $(dirname $0) ftp.ebi.ac.uk:/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007847/Suzuki_T2D_README.txt

