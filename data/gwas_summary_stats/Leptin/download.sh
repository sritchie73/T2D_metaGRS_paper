#!/bin/bash

wget -P $(dirname $0) ftp.ebi.ac.uk:/pub/databases/gwas/summary_statistics/GCST90007001-GCST90008000/GCST90007310/GCST90007310_buildGRCh37.tsv
gzip $(dirname $0)/GCST90007310_buildGRCh37.tsv
