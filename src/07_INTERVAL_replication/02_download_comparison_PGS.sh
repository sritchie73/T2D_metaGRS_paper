#!/bin/bash

# T2D PGS that are more recent, but can't be evaluated in UK Biobank due to being built from 
# GWAS summary statistics that include UK Biobank participants
mkdir -p data/T2D_PGS/UKB_invalid

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS000864 --name Aly2021 --out data/T2D_PGS/UKB_valid/Aly2021

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS000330 --name Mars2020 --out data/T2D_PGS/UKB_invalid/Mars2020

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS000868 --name Aksit2020 --out data/T2D_PGS/UKB_invalid/Aksit2020

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS000713 --name SinnottArmstrong2021 --out data/T2D_PGS/UKB_invalid/SinnottArmstrong2021

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS000729 --name Ritchie2021 --out data/T2D_PGS/UKB_invalid/Ritchie2021

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS000805 --name Polfus2021 --out data/T2D_PGS/UKB_invalid/Polfus2021

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS001818 --name Prive2022_PLR --out data/T2D_PGS/UKB_invalid/Prive2022_PLR

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS002026 --name Prive2022_LDpred2 --out data/T2D_PGS/UKB_invalid/Prive2022_LDpred2

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS003103 --name Ma2022 --out data/T2D_PGS/UKB_invalid/Ma2022

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS002780 --name Wong2022 --out data/T2D_PGS/UKB_invalid/Wong2022

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS002671 --name Weissbrod2022 --out data/T2D_PGS/UKB_invalid/Weissbrod2022

