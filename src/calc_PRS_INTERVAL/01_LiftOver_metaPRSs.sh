#!/opt/homebrew/bin/bash

# Download chain file (and wget via homebrea)
mkdir -p chain_files/
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz -P chain_files/

# Needs pygscatalog, on OSX:
# pip3 install --user --break-system-packages pgscatalog-core
# export PATH=$HOME/Library/Python/3.14/bin:$PATH
mkdir -p output/AoU_training/scorefiles_grch37/
pgscatalog-format -s output/AoU_training/scorefile_formatted/*.txt.gz \
  --chain_dir chain_files/ \
  --liftover --target_build GRCh37 \
  --outfile output/AoU_training/scorefiles_grch37/ \
  --threads 4 --verbose

# Upload to CSD3
scp output/AoU_training/scorefiles_grch37/* csd3:projects/T2D_metaGRS/data/PRS_score_files/

  
