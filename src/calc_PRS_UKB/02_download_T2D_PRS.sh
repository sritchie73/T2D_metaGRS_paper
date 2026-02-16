#!/opt/homebrew/bin/bash

# Needs pygscatalog, on OSX:
# pip3 install --user --break-system-packages pgscatalog-core
# export PATH=$HOME/Library/Python/3.14/bin:$PATH

mkdir -p PGS_Catalog/MONDO_0005148
pgscatalog-download --efo MONDO_0005148 --build GRCh38 --outdir PGS_Catalog/MONDO_0005148/ --verbose

# And DNAnexus command line tools, on OSX (with homebrew):
# brew install dxpy
dx upload PGS_Catalog/MONDO_0005148/* --destination PRS_score_files/

# No need to keep local copy
rm -rf PGS_Catalog
