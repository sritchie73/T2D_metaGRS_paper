# Integrated clinical risk prediction of type 2 diabetes with a multifactorial polygenic risk score 

This repository houses and documents the code used to generate the results in the study Ritchie SC *et al.* Integrated clinical risk prediction of type 2 diabetes with a multifactorial polygenic risk score. medRxiv, doi: [10.1101/2024.08.22.24312440](https://www.medrxiv.org/content/10.1101/2024.08.22.24312440v1)

## Repository information

The purpose of this repository is to provide a public record of the source code (i.e. methods) used for the entitled manuscript. The source code provided in this repository has not been designed to regenerate the results as-is for third-parties. The scripts herein contain numerous hard-coded filepaths and rely on data that cannot be made publicly available through this repository. 

The analysis scripts used to generate the results in our paper are contained in the `src/` folder. The `ext_src/` folder contains additional scripts primarily used for extraction and curation of UK Biobank and other data used for this project. This delineation has been made for pragmatic purposes: the scripts stored in `ext_src/` are copies of scripts located elsewhere on our HPC cluster as they have been written for cross-project purposes (see `ext_src/README.txt` for details).

## Software and versions used

The following software and versions were used to run these scripts:

- Rocky Linux release 8.10 (Green Obsidian) (HPC operating system)
- Slurm version 23.02.07 (HPC queue manager and job submission system)
- GNU bash version 4.4.20(1) (shell environment)
- docopts v0.6.41 commit ccd24b6 (part of some of the commandline tools in `ext_src/`) 
- PLINK v2.00a5.7LM AVX2 Intel (30 Oct 2023) (www.cog-genomics.org/plink/2.0/), aliased as plink2 (used to curate UKB genotype data and calculate PRS)
- UKBiobank ukbconv_lx (c) CTSU. Compiled Mar 14 2018 (Used for converting UK Biobank data to csv)
- R version 4.3.1 (2023-06-16), along with the R packages:
  - Data wrangling:
    - data.table version 1.14.8
    - bit64 version 4.0.5
    - lubridate version 1.9.3
    - openxlsx version 4.2.5.2
    - R.utils version 2.12.2
    - Rmpfr version 0.9-5
  - Programming:
    - docopt version 0.7.1
    - foreach version 1.5.2
    - doMC version 1.3.8
  - Statistics:
    - bigsnpr version 1.12.2
    - MASS version 7.3-60
    - RNOmni version 1.0.1.2
    - survival version 3.5-7
    - boot version 1.3-28.1
    - glmnet version 4.1-8
    - caret version 6.0-94
    - nricens version 1.6
    - pROC version 1.18.5
    - QDiabetes version 1.0-2
  - Visualisation:
    - ggplot2 version 3.4.4
    - ggstance version 0.3.6
    - ggthemes version 4.2.4
    - ggbeeswarm version 0.7.2
    - ggforce version 0.4.2 
    - gghalves version 0.1.4
    - scales version 1.2.1
    - RColorBrewer version 1.1-3
    - cowplot version 1.1.1
    - patchwork version 1.2.0
    - hexbin version 1.28.3

Inkscape version 1.2 was used to layout and annotate figures from the figure components generated within the R scripts. Microsoft Office 365 was used to draft the manuscript (Microsoft Word) and curate supplemental tables (Microsoft Excel) on MacOS Ventura 13.6.9
