library(data.table)
library(R.utils)
library(dxutils) # remotes::install_github("sritchie73/dxutils")

# Get list of score files to convert to PGS Catalog format
scores <- list.files(path="AoU_training/metaPRS_training", pattern="*.txt.gz", 
  recursive=TRUE, full.names=TRUE)

# Create output directory for upload
out_dir <- "AoU_training/scorefile_formatted"
dir.create(out_dir)

# Iterate through scores and create score-file formatted version
for (sf in scores) {
  # Create header information 
  # (see https://pgsc-calc.readthedocs.io/en/v3-alpha.1/howto/custom)
  out_file <- file.path(out_dir, gsub(".gz", "", basename(sf)))
  score_name <- gsub(".txt.gz", "", basename(sf))
  
  cat("#pgs_name=", score_name, "\n", file=out_file)
  cat("#pgs_id=", score_name, "\n", file=out_file, append=TRUE)
  cat("#trait_reported=type 2 diabetes mellitus\n", file=out_file, append=TRUE)
  cat("#genome_build=GRCh38\n", file=out_file, append=TRUE)
      
  # Reformat and add weights
  weights <- fread(sf)
  weights <- weights[!is.na(weight) & weight != 0]
  weights <- weights[, .(chr_name=chr, chr_position=pos, effect_allele, 
    other_allele, effect_weight=weight)]
  write.table(weights, sep="\t", quote=FALSE, file=out_file, append=TRUE, row.names=FALSE)
  
  # gzip
  gzip(out_file)
}

# Upload to DNAnexus project storage 
dx_upload("AoU_training/metaPRS_training/", "PRS_score_files/")
