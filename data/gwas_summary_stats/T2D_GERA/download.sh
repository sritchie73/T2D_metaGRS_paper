#!/bin/bash

echo "Warning: downloading results for all traits analysed in GERA (25GB)"
echo "This will take several hours (or more) depending on your connection speed"

wget -P $(dirname $0) https://personal.broadinstitute.org/ryank/Results_GERA.zip 
unzip $(dirname $0)/Results_GERA.zip Results_DIA2_GERA_chr_1_to_22.txt
gzip $(dirname $0)/Results_DIA2_GERA_chr_1_to_22.txt
unzip -a -p $(dirname $0)/Results_GERA.zip README.txt > Results_GERA_README.txt
rm $(dirname $0)/Results_GERA.zip

