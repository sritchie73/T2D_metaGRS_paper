#!/bin/bash

wget -P $(dirname $0) ftp.ebi.ac.uk:/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004773/METAANALYSIS_DIAGRAM_SE1.txt
wget -P $(dirname $0) ftp.ebi.ac.uk:/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004773/DIAGRAM_1000G_GWAS.pdf
