#!/bin/bash

wget -P $(dirname $0) https://humandbs.biosciencedbc.jp/files/hum0014/hum0014_v20_gwas_v1.zip
unzip $(dirname $0)/hum0014_v20_gwas_v1.zip
rm $(dirname $0)/hum0014_v20_gwas_v1.zip
gzip $(dirname $0)/BBJCAD_2020.sumstats

wget -P $(dirname $0) https://humandbs.biosciencedbc.jp/files/hum0014/hum0014_v20_README_CAD_GWAS.txt
