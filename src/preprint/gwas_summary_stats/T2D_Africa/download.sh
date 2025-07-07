#!/bin/bash

wget -P $(dirname $0) ftp.ebi.ac.uk:/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008114/ChenJ_31049640
wget -P $(dirname $0) ftp.ebi.ac.uk:/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008114/ChenJ_31049640.readme
gzip $(dirname $0)/ChenJ_31049640
