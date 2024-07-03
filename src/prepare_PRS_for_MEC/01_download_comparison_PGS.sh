#!/bin/bash

mkdir -p data/T2D_PGS/MEC_GRCh38

# From the UKB_valid download script
./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS000014 --harmonized GRCh38 --name Khera2018 --out data/T2D_PGS/MEC_GRCh38/Khera2018

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS000036 --harmonized GRCh38 --name Mahajan2018 --out data/T2D_PGS/MEC_GRCh38/Mahajan2018

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS001357 --harmonized GRCh38 --name Ye2021 --out data/T2D_PGS/MEC_GRCh38/Ye2021

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS001781 --harmonized GRCh38 --name Tamlander2022 --out data/T2D_PGS/MEC_GRCh38/Tamlander2022

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS002771 --harmonized GRCh38 --name Mars2022_AJHG --out data/T2D_PGS/MEC_GRCh38/Mars2022_AJHG

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS002243 --harmonized GRCh38 --name Mars2022_CellGenomics --out data/T2D_PGS/MEC_GRCh38/Mars2022_CellGenomics

# From the UKB_invalid download script 
./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS000864 --harmonized GRCh38 --name MonsourAly2021 --out data/T2D_PGS/MEC_GRCh38/MonsourAly2021

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS000330 --harmonized GRCh38 --name Mars2020 --out data/T2D_PGS/MEC_GRCh38/Mars2020

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS000868 --harmonized GRCh38 --name Aksit2020 --out data/T2D_PGS/MEC_GRCh38/Aksit2020

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS000713 --harmonized GRCh38 --name SinnottArmstrong2021 --out data/T2D_PGS/MEC_GRCh38/SinnottArmstrong2021

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS000729 --harmonized GRCh38 --name Ritchie2021 --out data/T2D_PGS/MEC_GRCh38/Ritchie2021

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS000805 --harmonized GRCh38 --name Polfus2021 --out data/T2D_PGS/MEC_GRCh38/Polfus2021

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS001818 --harmonized GRCh38 --name Prive2022_PLR --out data/T2D_PGS/MEC_GRCh38/Prive2022_PLR

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS002026 --harmonized GRCh38 --name Prive2022_LDpred2 --out data/T2D_PGS/MEC_GRCh38/Prive2022_LDpred2

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS003103 --harmonized GRCh38 --name Ma2022 --out data/T2D_PGS/MEC_GRCh38/Ma2022

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS002780 --harmonized GRCh38 --name Wong2022 --out data/T2D_PGS/MEC_GRCh38/Wong2022

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS002671 --harmonized GRCh38 --name Weissbrod2022 --out data/T2D_PGS/MEC_GRCh38/Weissbrod2022

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS003867 --harmonized GRCh38 --name Shim2023 --out data/T2D_PGS/MEC_GRCh38/Shim2023

# Three PRS that need to be combined in R script downstream:
./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS003443 --harmonized GRCh38 --name HuertaChagoya2023_EUR --out data/T2D_PGS/MEC_GRCh38/HuertaChagoya2023_EUR

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS003444 --harmonized GRCh38 --name HuertaChagoya2023_EAS --out data/T2D_PGS/MEC_GRCh38/HuertaChagoya2023_EAS

./src/PGS_resources/download_PGS_from_catalog.sh --nouuid --overwrite \
  --pgsid PGS003445 --harmonized GRCh38 --name HuertaChagoya2023_LAT --out data/T2D_PGS/MEC_GRCh38/HuertaChagoya2023_LAT

