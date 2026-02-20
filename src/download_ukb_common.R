library(data.table)

download_files <- function(pattern) {
  files <- system(sprintf("dx find data --name '%s' --path common", pattern), intern=TRUE) 
  files <- tstrsplit(files, "(bytes */)|(KB */)|( \\(file-)")[[2]]
  for (ff in files) {
    if (!file.exists(ff)) {
      system(sprintf("mkdir -p '%s'", dirname(ff)))
      system(sprintf("dx download '%s' -o '%s/'", ff, dirname(ff)))
      cat(sprintf("Downloaded '%s' to '%s/'\n", ff, dirname(ff))) 
    }
  }
}

download_files("*.R")
download_files("README.txt")
