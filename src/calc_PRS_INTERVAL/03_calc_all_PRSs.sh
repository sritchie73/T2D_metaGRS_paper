#!/bin/bash

bash src/PGS_resources/CSD3/calc_PS_lvls.sh \
  --score-file 'data/PRS_score_files/' \
  --type 'd' \
  --out 'data/INTERVAL/PRS_levels' \
  --single-out 'all_PRSs' \
  --work 'data/INTERVAL/PRS_levels'
  