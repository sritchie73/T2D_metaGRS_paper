#!/bin/bash

# Script designed to be run via 'dx run swiss-army-knife' to compute all 
# ancestry-specific metaPRS and multi-ancestry metaPRS alongside all T2D PRS
# currently in the PGS Catalog on the UK Biobank RAP
# Following  https://pgsc-calc.readthedocs.io/en/v3-alpha.1/tutorial/ukbrap

# Skip the dx-fuse installation step, we will use -imount_inputs=true in dx run
# to automatically mount the project to /mnt/project

# Install nextflow
sudo apt-get update && sudo apt-get install -y default-jre
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin
nextflow run hello

# Install the PGS Catalog Calculator
git clone -b v3-alpha.1 https://github.com/PGScatalog/pgsc_calc.git
cd pgsc_calc
nextflow run main.nf -profile test_full,docker

# Create configuration file
echo "process {
    withName: 'PGSC_CALC_FORMAT' {
        cpus = 2
        memory = { 10.GB * task.attempt }
    }
    withName: 'PGSC_CALC_LOAD' {
        cpus = 1
        memory = { 16.GB * task.attempt }
    }
    withName: 'PGSC_CALC_SCORE' {
        cpus = 4
        memory = { 60.GB * task.attempt }
    }
}
params {
    variant_batch_size = 1000
}
" > ${HOME}/ukb_config.config

# Prepare the sample sheet
find '/mnt/project/Bulk/Imputation/Imputation from genotype (TOPmed)' -regex '.*c[0-9]+.*.bgen$' \
  | awk -v OFS="," '
BEGIN { ORS = ""; print " [ "}
{
    full = $0

    # extract number after "c"
    chrom = ""
    if (match(full, /c[0-9]+/)) {
        chrom = substr(full, RSTART+1, RLENGTH-1)
    }

    # get .sample path
    sample = full
    sub(/\.[^.\/]+$/, ".sample", sample)

    printf "%s{\"sampleset\": \"%s\", \"path\": \"%s\", \"chrom\": \"%s\",  \"file_format\": \"%s\",  \"genotyping_method\": \"%s\",  \"bgen_sample_file\": \"%s\"}",
        separator, "ukb_topmed", full, chrom, "bgen", "array", sample
        separator=", "
}
END { print " ] " }
' | jq > ${HOME}/samplesheet.json

# Run PGS calculator
nextflow run main.nf \
  -profile docker \
  --input $HOME/samplesheet.json \
  --target_build GRCh38 \
  -c $HOME/ukb_config.config \
  --outdir ${HOME} \
  --efo_id MONDO_0005148 \
  --scorefile /mnt/project/AoU_training/scorefile_formatted/*
  
# Clean up files we don't want to upload to project storage
rm -rf ${HOME}/pgsc_calc
rm ${HOME}/samplesheet.json
rm ${HOME}/ukb_config.config
  