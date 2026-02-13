#!/bin/bash

bash src/PGS_resources/UKB_RAP/calc_PS_lvls.sh \
  --score-file 'PRS_score_files' \
  --type 'd' \
  --out 'PRS_levels' \
  --single-out 'all_PRSs' \
  --work 'PRS_levels/'
