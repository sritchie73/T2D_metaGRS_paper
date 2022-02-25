#!/bin/bash

# Have a separate directory later for T2D PGS we can't compare in UKB due to being built
# from UKB GWAS summary statistics, but can compare in other cohorts (e.g. INTERVAL).
mkdir -p data/T2D_PGS/ukb_valid

./src/PGS_resources/download_PGS_from_catalog.sh \
  --pgsid PGS000014 --name Khera2018 --out data/T2D_PGS/ukb_valid/Khera2018

./src/PGS_resources/download_PGS_from_catalog.sh \
  --pgsid PGS000036 --name Mahajan2018 --out data/T2D_PGS/ukb_valid/Mahajan2018

./src/PGS_resources/download_PGS_from_catalog.sh \
  --pgsid PGS000864 --name Aly2021 --out data/T2D_PGS/ukb_valid/Aly2021

./src/PGS_resources/download_PGS_from_catalog.sh \
  --pgsid PGS001357 --name Ye2021 --out data/T2D_PGS/ukb_valid/Ye2021

./src/PGS_resources/download_PGS_from_catalog.sh \
  --pgsid PGS001781 --name Tamlander2022 --out data/T2D_PGS/ukb_valid/Tamlander2022
