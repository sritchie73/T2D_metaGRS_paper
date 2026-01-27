#!/bin/bash

# Intended to be run on the local machine where the app is built

wget https://github.com/dnanexus/dxfuse/releases/download/v1.6.1/dxfuse-linux
chmod +x dxfuse-linux
mv dxfuse-linux ../bin/dxfuse

wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250819.zip
unzip plink_linux_x86_64_20250819.zip
mv plink prettify ../bin/
rm toy.map toy.ped LICENSE plink_linux_x86_64_20250819.zip

wget https://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_20260110.zip
unzip plink2_linux_avx2_20260110.zip
mv plink2 vcf_subset ../bin/
rm intel-simplified-software-license.txt plink2_linux_avx2_20260110.zip
