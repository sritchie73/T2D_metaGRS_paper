#!/bin/bash

wget -P $(dirname $0) ftp.ebi.ac.uk:/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010118/SpracklenCN_prePMID_T2D_ALL_Primary.txt.gz
wget -P $(dirname $0) ftp.ebi.ac.uk:/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010118/SpracklenCN_prePMID_T2D_ALL_Primary_README.txt

