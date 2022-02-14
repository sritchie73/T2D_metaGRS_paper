library(data.table)
library(foreach)
library(doMC)

# Auto detect if we're on compute nodes and determine available CPUs for parallization
nCores <- Sys.getenv("SLURM_CPUS_ON_NODE")
if (nCores == "") {
  registerDoSEQ()
  setDTthreads(1L)
} else {
  nCores <- as.integer(nCores)
  registerDoMC(nCores)
  setDTthreads(nCores)
}

# Determine mem per CPU available
partition <- Sys.getenv("SLURM_JOB_PARTITION")
if (partition %in% c("cclake", "skylake")) {
  memPerCore <- 6*1024
} else if (partition %in% c("cclake-himem", "skylake-himem")) {
  memPerCore <- 11*1024
} else if (partition == "cardio") {
  memPerCore <- 4*1024
} else {
  memPerCore <- NA
}
