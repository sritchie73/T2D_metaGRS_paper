#!/usr/bin/env bash
src_dir=/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/GRS_resources

eval "$($src_dir/software_dependencies/docopts -h - : "$@" <<EOF
Downloads a PGS from the catalog 

Given a PGS ID from the PGS Catalog, downloads and saves the score file

Usage:
  download_PGS_from_catalog.sh --pgsid <PGS_ID> [options]
  download_PGS_from_catalog.sh -h | --help

Options:
  -h --help                   Show this screen.
  --pgsid <PGS_ID>            Identifier in the PGS Catalog.
  --name <name>               Optional additional name to give the PGS folder and score file.
                              If not provided, the reported trait in the PGS score file header 
                              will be used. [default: NULL]
EOF
)"

# Download score
tmp_dir=$HOME/rds/hpc-work/GRS_resources_v2_tmpwd/PGS_catalog_downloads
mkdir -p $tmp_dir
wget -P $tmp_dir ftp://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/$pgsid/ScoringFiles/$pgsid.txt.gz

# Get trait name from header if --name not provided
if [ $name == "NULL" ]; then
  trait=$(zcat $tmp_dir/$pgsid.txt.gz | grep -w -m 1 'trait_reported' | sed 's/^.*=//') # grab the line starting '# Reported Trait = ' then remove that part of the line
  name=$(echo $trait | sed "s/'//g" | sed 's/\b\(.\)/\u\1/g' | sed 's/[^[:alnum:]]//g') # convert first letter of each word to uppercase, then strip out non-alphanumeric characters
fi

# Set unique identifier for the PGS
uuid=$(uuidgen | sed 's/-.*//')

# Build the name of the score
pgs_name="${name}_${pgsid}_${uuid}"

# Save the score
mkdir -p $src_dir/$pgs_name
mv $tmp_dir/$pgsid.txt.gz $pgs_name/$pgs_name.txt.gz
chmod -w $src_dir/$pgs_name/$pgs_name.txt.gz

