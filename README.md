# Integrated clinical risk prediction of type 2 diabetes with a multifactorial polygenic risk score 

This repository houses and documents the code used to generate the results in the study Ritchie SC *et al.* Integrated clinical risk prediction of type 2 diabetes with a multifactorial polygenic risk score. medRxiv, doi: [10.1101/2024.08.22.24312440](https://www.medrxiv.org/content/10.1101/2024.08.22.24312440v1)

## Repository information

The purpose of this repository is to provide a public record of the source code (i.e. methods) used for the entitled manuscript. Importantly, the source code provided in this repository has not been designed to regenerate the results as-is for third-parties. The scripts herein contain numerous hard-coded filepaths and rely on data that cannot be made publicly available through this repository and also could not be analysed in a single location, as each dataset had to be analysed within its own dedicated trusted research environments, each with their own organisational quirks and dramatically different compute configurations and by different analysts with their own coding practices and levels of experience. 

## Compute systems used

The following compute systems were used for this project:

 (1) The Cambridge Service for Data Driven Discovery (CSD3) high performance computing (HPC) cluster: https://docs.hpc.cam.ac.uk/hpc/index.html
 (2) The All of Us (AoU) Research Platform : https://www.researchallofus.org/data-tools/workbench/
 (3) The UK Biobank Research Analysis Platform (UKB RAP): https://dnanexus.gitbook.io/uk-biobank-rap
 (4) The National University of Singapore High Performance Computing Facility (NUS HPC): https://research.nus.edu.sg/research-facilities/project/central-high-performance-computing-facility/
 (5) Local compute
 
CSD3 was used to prepare GWAS summary statistics for T2D multiancestry metaPRS training and for running analysis on the INTERVAL cohort. In the original preprint, CSD3 was also used to analyse the UK Biobank cohort (before use of the RAP was mandated) and to generate tables and figures.
 
The AoU Research Platform was used for all analyses of the All of Us cohort, which included multi-ancestry metaPRS training and initial validation and comparison to other published PRSs in independent samples not used for metaPRS training.

The UKB RAP was used for all analyses of the UK Biobank cohort, which included independent replication of multi-ancestry metaPRS performance, comparison to published PRSs, and assessing benefits for screening and 10-year risk prediction in comparison to the QDiabetes risk score.

Local compute was used to collate summary statistics and prepare tables and figures for publication, as well as manuscript drafting.

## Repository Organisation

Code in this repository lacks a coherent centralized organisation in part due to the different architectures of compute systems involved, long running and evolving nature of the project, independent analysts involved, and presence of dataset QC and processing pipelines across many independent projects.

Broadly speaking, the logical ordering of the source code (in terms of analysis sequence) is as follows:

 (1) Code for downloading GWAS summary statistics used as inputs for the multiancestry metaPRS training can be found under `data/gwas_summary_statistics/`. Note some GWAS summary statistics required manual rather than automated download; see Table S14 in the manuscript for further details.
 
 (2) Code for filtering and harmonizing the GWAS summary statistics for the multiancestry metaPRS training can be found under `src/prepare_sumstats/`
 
 (3) Code for metaPRS training can be found in `src/metaPRS_training/` 

 (4) Code for computing metaPRS and other T2D PRS in UK Biobank can be found in `src/calc_PRS_UKB/`
 
 (5) Code for curating common phenotypes (e.g. cross-project) in UK Biobank on the RAP can be found in `common/`
 
 (6) Code for testing associations between PRS and T2D case status in UK Biobank can be found in `src/validation_UKB/`
 
 (7) Code for computing metaPRS and other T2D PRS in INTERVAL can be found in `src/calc_PRS_INTERVAL/`
 
 (8) Code for testing associations between PRS and T2D case status in UK Biobank can be found in `src/validation_INTERVAL/`
 

 
 

## Software and versions used

The following software and versions were used to run these scripts:

### Local compute

 - OSX (Tahoe 26.2)
 - homebrew for managing software installation
 - RStudio
 - GNU bash version 5.3.9(1) (shell environment)
 - golang (for compiling docopts)
 - docopts
 - python 3.14 along with libraries:
   - dxpy
   - pgscatalog-core
 - Inkscape was used to layout and annotate figures from the figure components generated within the R scripts.
 - Microsoft Office 365 was used to draft the manuscript (Microsoft Word) and curate supplemental tables (Microsoft Excel)
 
### The All of Us research workbench


 
### UKB Research Analysis Platform

- The run_script applet for running scripts as jobs (https://github.com/sritchie73/dxapplet-run_script)
- The dxutils R package providing extensions that enhance the DNAnexus command line utilities (https://github.com/sritchie73/dxutils)

### CSD3

- Rocky Linux release 8.10 (Green Obsidian) (HPC operating system)
- Slurm version 23.02.07 (HPC queue manager and job submission system)
- GNU bash version 4.4.20(1) (shell environment)
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

