########################################################################################
# Load R package dependencies
########################################################################################

suppressMessages(library("data.table"))
suppressMessages(library("foreach"))
suppressMessages(library("docopt"))
suppressMessages(library("bit64"))
suppressMessages(library("dxutils"))

########################################################################################
# Parse user input and sanity check
########################################################################################

# Parse input arguments
"Calculate the levels of a polygenic score in a group of samples

Note this script is intended to be run within a compute node on as part of an array job
with one job per chromosome set by the environment variable SLURM_ARRAY_TASK_ID

Usage:
  calc_PS_lvls.R --score-file <file> --work <directory> [options]
  calc_PS_lvls.R -h | --help

Options:
  -h --help                   Show this screen.
  --score-file <file>         Path to polygenic score file, directory, or file containing list
                              of score files (see --type) on the RAP project storage.
  --work <directory>          Working directory on RAP project storage to store intermediate files 
                              and logs shared across all scores for the duration of the run. 
  --type <type>               Type of filepath given to --score-file. 's': filepath points to
                              a single polygenic score. 'd': filepath is a directory on the RAP
                              project storage containing multiple score files. 'l': filepath points 
                              to a file on RAP project storage containing a list of paths, one per 
                              line. Each line may optionally be followed by three to ten column 
                              position or column name arguments to specify for each file the 
                              --score-X options (described below) to handle different locations, 
                              names, or presence/absence of the respective fields in the score files. 
                              If not provided, the columns for that particular score will be loaded 
                              as provided to the respective program arguments. Reasonable combinations 
                              can be inferred when listing only a subset of columns (assuming they are 
                              in the same order as documented below) unless you are using the ambiguous 
                              MAF threshold option, or there are columns indicating that some 
                              or all effect weights should be taken as either dominant or recessive. 
                              The score_summary.txt output will always report which columns have been used 
                              for which field. Note that files from the PGS Catalog are automatically
                              detected and handled, so you do not need to specify these manually. [default: s]
  --score-rsid <col>          Name or number of the column in the polygenic score file corresponding to
                              the variant marker id to pass to plink. If column is not present, set to
                              'NULL'. [default: rsid]
  --score-chr <col>           Name or number of the column in the polygenic score file corresponding to
                              the variant's chromosome. If column is not present, set to 'NULL'.
                              If the score contains non-autosomal variants, then the chromosome field
                              *must* contain one of X, Y, XY, or MT. [default: chr_name]
  --score-pos <col>           Name or number of the column in the polygenic score file corresponding to
                              the variant's position. If column is not present, set to 'NULL'. [default: chr_position]
  --score-EA <col>            Name or number of the column in the polygenic score file corresponding to
                              the variant's effect allele. [default: effect_allele]
  --score-EAF <col>           Name or number of the column in the polygenic score file corresponding to
                              the effect allele frequency (i.e. if using --keep-ambiguous and --ambiguous-thresh).
                              If column is not present, set to 'NULL'. [default: NULL]
  --score-OA <col>            Name or number of the column in the polygenic score file corresponding to
                              the variant's non-effect allele. If column is not present, set to 'NULL'.
                              [default: other_allele]
  --score-weight <col>        Name or number of the column in the polygenic score file corresponding to
                              the effect allele's weight in the polygenic score. If the score file has
                              multiple weight columns for multiple scores, set this to 'm' to calculate
                              the levels of all these scores. In this case, the program will assume that
                              all columns not listed in the arguments above are weights columns. [default: effect_weight]
  --score-dominant <col>      Column of TRUE/FALSE values indicating whether the effect for each variant
                              should be considered dominant (i.e. the weight is multiplied by the effect
                              allele presence/absence rather than by the number of copies). [default: NULL]
  --score-recessive <col>     Column of TRUE/FALSE values indicating whether the effect for each variant
                              should be considered recessive (i.e. the weight is counted only when there
                              are two copies of the effect allele). [default: NULL]
  --match-by-rsid             Flag to indicate that the score file should be matched to the 
                              genotype data by the rsid column instead of by chromosome and
                              position when all three columns are provided.
  --cohort-name <name>        Name of the group of samples, used to name the output file name
                              and folder (if --out not provided). [default: UKB_TopMed]
  --out <directory>           Directory on RAP project storage to save the final results to. By default, 
                              the results are stored in a folder named <--cohort-name>_sample_levels/ 
                              in the same directory as each input --score-file. If multiple score files 
                              are detected (files ending in .txt.gz) in a score file's directory then a 
                              further sub-folder is created <--cohort-name>_sample_levels/<score_name>/ 
                              for the score being calculated by this program. [default: NULL]
  --single-out <name>         If provided, saves only a single file on RAP project storage containing the 
                              levels of all scores being calculated in a file with the given name (i.e. 
                              <--out>/<--single-out>.sscore.gz. The default behaviour is to otherwise save 
                              the levels of each score into separate files. [default: NULL]
  --genotype-prefix <prefix>  Path and prefix occurring before the chromosome number for the genotype
                              data on RAP project storage you want to use for polygenic scoring.
                              Defaults to the plink2 binary files extracted for the TopMed imputation
                              in the CEU_overarching project.
                              [default: common/Imputed Genotypes/TopMed\ (GRCh38)/ukb21007_c]
  --genotype-suffix <suffix>  Suffix for the filename occurring after the chromosome number but before the
                              .pgen/.pvar/.psam extension for the genotype data. [default: _b0_v1]
  --single-geno               Flag to indicate that the genotype data is stored as a single file, not split across
                              multiple chromosomes. In this case, you can ignore the --genotype-suffix argument.
  --genotype-format <format>  Format the genotype data is stored in, must correspond to one of the arguments to
                              plink, e.g. the default, 'pfile' is passed directly to plink as '--pfile'. To use
                              plink version 1 binary data (bed/bim/fam) set this as 'bfile'. [default: pfile]
  --keep <file>               Optional, path to file on RAP project storage to pass to plink2 --keep to subset 
                              to a given set of samples when calculating the polygenic score levels. [default: NULL]
  --keep-ambiguous            Flag to force the program to keep variants with ambiguous alleles,
                              (e.g. A/T and G/C SNPs), which are normally excluded. In this case
                              the program proceeds assuming that the genotype data is on the
                              same strand as the GWAS whose summary statistics were used to 
                              construct the score.
  --ambiguous-thresh <maf>    When provided, then variants with ambiguous alleles are kept only when their
                              minor allele frequency is below this threshold. Requires each score file to
                              have a column giving the effect allele's frequency. If multiple scores, then
                              any scores missing this column will exclude all variants with ambiguous alleles.
                              If using, we would suggest a threshold of 0.42. [default: NULL]
  --freqx-prefix <prefix>     Optional, path on RAP project storage and prefix occurring before the chromosome
                              number in the filepaths for the plink1.9 --freqx reports containing the allele 
                              frequencies. These may be used instead of those directly estimated from the 
                              genotype data, e.g. when matching ambiguous alleles, or when imputing missing 
                              alleles when using plink1 binary data. If --single-geno has been set, a single 
                              freqx file is also assumed. [default: NULL]
  --freqx-suffix <suffix>     Optional, suffix for the above plink1.9 --freqx reports. [default: NULL]
  --remove-multiallelic       Flag that controls whether multi-allelic variants are kept or not. If given,
                              variants that have > 2 alleles in either the score file or genotype data are
                              discarded.
" -> doc

# Parse arguments and make additional sanity checks so errors are
# caught before long-running code fails
args <- docopt(doc)

# Set any "NULL" arguments to NULL
for (an in names(args)) {
  if (is.null(args[[an]]) || args[[an]] == "NULL") {
    args[[an]] <- NULL
  }
}

# Make sure we're running as part of an array job
if (Sys.getenv("SLURM_ARRAY_TASK_ID") == "") {
  stop("This script should not be run directly, use calc_PS_lvls.sh")
}

# Make local working directories for carrying out the work
dir.create("input_data")
dir.create("output")
dir.create("work")
dir.create("checkpointing")
dir.create("errors")
dir.create("finished")

# Make sure the work directory on project storage is interpreted as a directory
args[["--work"]] <- paste0(args[["--work"]], "/")
args[["--work"]] <- gsub("//", "/", args[["--work"]])
args[["work"]] <- args[["--work"]]

# Save arguments for easier debugging
saveRDS(args, file="checkpointing/args.rds")
dx_upload("checkpointing/args.rds", args[["--work"]], exists = "skip") # only 1 copy needed

# Check for --keep file if provided and download local copy
if (!is.null(args[["--keep"]])) assert_dx_exists(args[["--keep"]])

# Check type
stopifnot(args[["type"]] %in% c('s', 'd', 'l'))

# Check thresholds
if (!is.null(args[["--ambiguous-thresh"]])) {
	tryCatch({ args[["--ambiguous-thresh"]] <- as.numeric(args[["--ambiguous-thresh"]]) }, warning=function(w) { stop("--ambiguous-thresh must be numeric") })
	if ((args[["--ambiguous-thresh"]] < 0) || (args[["--ambiguous-thresh"]] > 0.5)) stop("--ambiguous-thresh must be between 0 and 0.5")
}

# Determine which chromosome we're running on:
task_to_chr <- function(num) {
  switch(num, "23" = "X", "24" = "Y", "25" = "XY", "26" = "MT", num)
}
taskIdx <- Sys.getenv("SLURM_ARRAY_TASK_ID")
chrIdx <- task_to_chr(taskIdx)
taskIdx <- as.integer(taskIdx)
taskMax <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_MAX"))

# Determine chromosome specific file paths.
if (args[["--genotype-format"]] == "pfile") {
	varfile <- paste0(args[["--genotype-prefix"]], ifelse(args[["--single-geno"]], "", chrIdx), args[["--genotype-suffix"]], ".pvar")
	genofile <- paste0(args[["--genotype-prefix"]], ifelse(args[["--single-geno"]], "", chrIdx), args[["--genotype-suffix"]], ".pgen")
	samplefile <- paste0(args[["--genotype-prefix"]], ifelse(args[["--single-geno"]], "", chrIdx), args[["--genotype-suffix"]], ".psam")
} else if (args[["--genotype-format"]] == "bfile") {
	varfile <- paste0(args[["--genotype-prefix"]], ifelse(args[["--single-geno"]], "", chrIdx), args[["--genotype-suffix"]], ".bim")
	genofile <- paste0(args[["--genotype-prefix"]], ifelse(args[["--single-geno"]], "", chrIdx), args[["--genotype-suffix"]], ".bed")
	samplefile <- paste0(args[["--genotype-prefix"]], ifelse(args[["--single-geno"]], "", chrIdx), args[["--genotype-suffix"]], ".fam")
} else {
  stop("Unsuported genotype format: ", args[["--genotype-format"]])
}
freqxfile <- paste0(args[["--freqx-prefix"]], ifelse(args[["--single-geno"]], "", chrIdx), args[["--freqx-suffix"]])
if (freqxfile == "" || freqxfile == chrIdx) {
  freqxfile <- NULL
}

# Check files exist (but don't download yet)
assert_dx_exists(varfile)
assert_dx_exists(genofile)
assert_dx_exists(samplefile)
if (!is.null(freqxfile)) assert_dx_exists(freqxfile)

# Check that another chromosome hasn't thrown an error that means we shouldn't
# proceed:
if (dx_exists(sprintf("%s/errors/score_summary.txt", args[["--work"]]))) {
  stop("All score files had errors, see ", args[["--work"]], "errors/score_summary.txt")
}

# Work out which checkpoint we're up to, if any
if (dx_exists(sprintf("%s/finished/score_summary_%s.txt", args[["--work"]], chrIdx))) {
  checkpoint <- 3
} else if (dx_exists(sprintf("%s/checkpoint2/score_summary_%s.txt", args[["--work"]], chrIdx))) {
  checkpoint <- 2
} else if (dx_exists(sprintf("%s/checkpoint1/score_info_chr%s.txt", args[["--work"]], chrIdx))) {
  checkpoint <- 1
} else {
  checkpoint <- 0
}

# Get information about all the score files
if (checkpoint < 1) {
  dir.create("checkpointing/checkpoint1")
  
  ########################################################################################
  # Now we need to construct a single score file for all input scores by matching to the 
  # designated genotype data. Each chromosome does this, this reduces memory load when
  # loading the score variants (rather than loading them all for each chromosome).
  ########################################################################################
  
  # Make sure the score file (or directory) exists
  assert_dx_exists(args[["--score-file"]])
  
  # construct table of score information - first step is to determine paths and
  # column positions / names for each relevant column
  is.integer <- function(value) {
    suppressWarnings(!is.na(as.integer(value)))
  }
  parsecolarg <- function(value) {
    if (is.null(value)) {
      return(NA)
    } else {
      return(value)
    }
  }
  if (args[["--type"]] == 's') {
  	score_info <- data.table(path = args[["--score-file"]], 
  	                        local_path = dx_download(args[["--score-file"]], "input_data/"),
  													rsid = parsecolarg(args[["--score-rsid"]]),
  													chr = parsecolarg(args[["--score-chr"]]),
  													pos = parsecolarg(args[["--score-pos"]]),
  													EA = parsecolarg(args[["--score-EA"]]),
  													EAF = parsecolarg(args[["--score-EAF"]]),
  													OA = parsecolarg(args[["--score-OA"]]),
  													weight = parsecolarg(args[["--score-weight"]]),
  													is_dom = parsecolarg(args[["--score-dominant"]]),
  													is_rec = parsecolarg(args[["--score-recessive"]]),
  													error = NA_character_)
  } else if (args[["--type"]] == 'd') {
    dx_download(paste0(args[["--score-file"]], "/"), "input_data/scores/")
    score_files <- setdiff(list.files(path="input_data/scores/"), list.dirs(path="input_data/scores/", full.names=FALSE, recursive=FALSE))
    score_info <- data.table(path = sprintf("%s/%s", args[["--score-file"]], score_files),
                            local_path = sprintf("input_data/scores/%s", score_files),
  													rsid = parsecolarg(args[["--score-rsid"]]),
  													chr = parsecolarg(args[["--score-chr"]]),
  													pos = parsecolarg(args[["--score-pos"]]),
  													EA = parsecolarg(args[["--score-EA"]]),
  													EAF = parsecolarg(args[["--score-EAF"]]),
  													OA = parsecolarg(args[["--score-OA"]]),
  													weight = parsecolarg(args[["--score-weight"]]),
  													is_dom = parsecolarg(args[["--score-dominant"]]),
  													is_rec = parsecolarg(args[["--score-recessive"]]),
  													error = NA_character_)
  	if (nrow(score_info) == 1 && is.na(score_info$path)) {
  		stop("Directory is empty: ", args[["--score-file"]])
  	}
  } else if (args[["--type"]] == 'l') {
    score_file_list <- dx_download(args[["--score-file"]], "input_data/")
    
  	tryCatch({ l <- readLines(score_file_list) }, error=function(e) { 
  		stop("Unable to read from file: ", args[["--score-file"]]) })
    
  	score_info <- foreach(line = l, .combine=rbind) %do% {
  		# Each line corresponds to a single score, optionally followed by 
  		# 3-7 columns indicating the --score-<X> columns. Where this number
  		# is < 7 reasonable defaults are inferred.
  		fields <- strsplit(line, "\\s+")[[1]] # split on all whitespace.
  		if (length(fields) == 1) { # Just a path to a score file, fill with defaults
  			return(data.table(path = fields[1],
  			                  local_path = dx_download(fields[1], "input_data/scores/", missing = "skip"),
  												rsid = parsecolarg(args[["--score-rsid"]]),
  												chr = parsecolarg(args[["--score-chr"]]),
  												pos = parsecolarg(args[["--score-pos"]]),
  												EA = parsecolarg(args[["--score-EA"]]),
  												EAF = parsecolarg(args[["--score-EAF"]]),
  												OA = parsecolarg(args[["--score-OA"]]),
  												weight = parsecolarg(args[["--score-weight"]]), 
  												is_dom = parsecolarg(args[["--score-dominant"]]),
  												is_rec = parsecolarg(args[["--score-recessive"]]),
  												error = NA_character_))
  		} else if (length(fields) == 10) { # All fields including dominant and recessive flags provided.
  			return(data.table(path = fields[1], rsid = fields[2], chr = fields[3],
  							pos = fields[4], EA = fields[5], EAF = fields[6],
  							OA = fields[7], weight = fields[8], is_dom = fields[9], is_rec=fields[10],
  							error = NA_character_))
  		} else if (length(fields) == 9) { 
  			return(data.table(path = fields[1], rsid = NA, chr = NA, pos = NA, EA = NA, EAF = NA, OA = NA, weight = NA, is_dom = NA, is_rec = NA,
  							error = "Unable to infer use case when 9 fields provided."))
  		} else if (length(fields) == 8) { # All fields provided, excluding dominant and recessive flags
  			return(data.table(path = fields[1], rsid = fields[2], chr = fields[3], 
  												pos = fields[4], EA = fields[5], EAF = fields[6],
  												OA = fields[7], weight = fields[8], is_dom = NA, is_rec = NA, error = NA_character_))
  		} else if (length(fields) == 7) {
  			return(data.table(path = fields[1], rsid = fields[2], chr = fields[3],
  							 pos = fields[4], EA = fields[5], EAF = NA, OA = fields[6], 
  							 weight = fields[7], is_dom = NA, is_rec = NA, error = NA_character_))
  		} else if (length(fields) == 6) {
  			return(data.table(path = fields[1], rsid = NA, chr = fields[2],
  							 pos = fields[3], EA = fields[4], EAF = NA, OA = fields[5], 
  							 weight = fields[6], is_dom = NA, is_rec = NA, error = NA_character_))
  		} else if (length(fields) == 5) {
  			return(data.table(path = fields[1], rsid = NA, chr = fields[2],
  							 pos = fields[3], EA = fields[4], EAF = NA, OA = NA, 
  							 weight = fields[4], is_dom = NA, is_rec = NA, error = NA_character_))
  		} else if (length(fields) == 4) {
  			return(data.table(path = fields[1], rsid = fields[2], chr = NA,
  							 pos = NA, EA = fields[3], EAF = NA, OA = NA, 
  							 weight = fields[4], is_dom = NA, is_rec = NA, error = NA_character_))
  		} else {
  			return(data.table(path = fields[1], rsid = NA, chr = NA, pos = NA, EA = NA, EAF = NA, OA = NA, weight = NA, is_dom = NA, is_rec = NA,
  												error = "At least three field to column mappings must be provided."))
  		} 
  	}
  } else {
  	stop("Internal error: allowed unknown --type")
  }
  
  # check all fields unique
  field_check <- suppressWarnings(melt(score_info[!is.na(error)], id.vars="path"))
  bad <- field_check[!is.na(value), .(fail=(length(unique(value)) == .N)), by=path][(fail)]
  score_info[field_check, on = .(path), error := "The same column cannot be used for multiple fields."]
  rm(bad, field_check)
  invisible(gc())
  
  # Check minimum set of fields present:
  score_info[!is.na(error) & (is.na(EA) | is.na(weight) | (is.na(rsid) & (is.na(chr) | is.na(pos)))),
    error := "Must provide effect allele, weight, and either the rsid, or chromosome and position columns"]
  
  # Give each score a unique number to uniquely identify it during computation
  score_info[, compName := as.character(.I)]
  
  # Reorganise
  score_info = score_info[, .(path, local_path, compName, rsid, chr, pos, EA, EAF, OA, weight, is_dom, is_rec, error)]
  
  # Check if each file exists
  score_info[is.na(error) & is.null(local_path), error := "score file does not exist"]
  
  # At this point, if all scores have errored, we can stop all tasks.
  if (all(!is.na(score_info$error))) {
    fwrite(score_info, sep="\t", quote=FALSE, file="errors/score_summary.txt")
    dx_upload("errors", args[["--work"]], exists="skip")
  	stop("All score files had errors, see ", args[["--work"]], "errors/score_summary.txt")
  }
  
  #################################################################################################
  # Now, attempt to load each score file and map the columns. At this point we also detect whether
  # they're PGS Catalog files and load appropriately
  #################################################################################################
  
  # Check that another chromosome hasn't thrown an error that means we shouldn't
  # proceed:
  if (dx_exists(sprintf("%s/errors/score_summary.txt", args[["--work"]]))) {
    stop("All score files had errors, see ", args[["--work"]], "errors/score_summary.txt")
  }
  
  rbindf <- function(...) rbind(..., fill=TRUE)
  scores <- foreach(idx = score_info[,.I], .combine=rbindf) %do% {
  	# previously failed scores are skipped
  	if (!is.na(score_info[idx, error])) {
  		return(NULL)
  	}
  
  	# Read header line so we can determine if its a PGS Catalog file (and check we can read the file)
  	tryCatch({ 
  		line1 <- readLines(score_info[idx, local_path], 1)
  	}, error = function(e) {
  		score_info[idx, error := "Could not open score file"]
  	})
  	if (score_info[idx, !is.na(error)]) {
  		return(NULL)
  	}
  
  	# Load the score with fread (and record an error if this fails)
  	tryCatch({
      score <- fread(cmd=sprintf("zgrep -v '^#' %s | zgrep -v '^[[:space:]]*$'", score_info[idx, local_path]))
  	}, error = function(e) {
  		score_info[idx, error := "Error when reading score with fread"]
  	})
  	if (score_info[idx, !is.na(error)]) {
  		return(NULL)
  	}
  
  	pgs_catalog_file <- grepl("### ?PGS CATALOG SCORING FILE", line1)
  	if (!pgs_catalog_file) {
  		# Set the column names based on the information in the score_info field
  		namecol <- function(score, old, new) {
  			if (!is.na(old)) {
  				if (is.integer(old)) {
            old <- as.integer(old)
  					if (names(score)[old] != new && new %in% names(score)) {
  						warning("Column named ", new, " found in score file (column ", 
  						paste(which(names(score) == new), collapse=", "), " but using column ", 
  						old, " [", names(score)[old], "] as --score-", new, " column instead")
  						# rename extra columns so we don't pick the wrong one later
  						which.new <- setdiff(which(names(score) == new), as.integer(old))
  						names(score)[which.new] <- paste0(new, ".", seq_along(which.new))
  					}
            score_info[idx, c(new) := names(score)[old]]
  					setnames(score, names(score)[old], new)
  				} else {
  					if (sum(names(score) == old) > 1) {
  						warning("Multiple columns named ", old, "in score file ", score_info[idx, path], " using first occurence as --score-", new, " column")
  						# rename extra columns so we don't pick the wrong one later
  						which.dup = which(names(score) == old)[-1]
  						names(score)[which.dup] <- paste0(old, ".", seq_along(which.dup))
  					} 
  					if (new != old && new %in% names(score)) {
  						warning("Column named ", new, " found in score file (column ", 
  						paste(which(names(score) == new), collapse=", "), " but using column ", 
  						old, " (column ", paste(which(names(score) == old), collapse=", "), 
  						") as --score-", new, " column instead")
  						# rename extra columns so we don't pick the wrong one later
  						which.new <- which(names(score) == new)
  						names(score)[which.new] <- paste0(new, ".", seq_along(which.new))
  					}
   
  					setnames(score, old, new, skip_absent=TRUE)
  				}
  			}
  		}
  		tryCatch({
  			namecol(score, score_info[idx, rsid], "rsid")
  			namecol(score, score_info[idx, chr], "chr")
  			namecol(score, score_info[idx, pos], "pos")
  			namecol(score, score_info[idx, EA], "EA")
  			namecol(score, score_info[idx, EAF], "EAF")
  			namecol(score, score_info[idx, OA], "OA")
  			namecol(score, score_info[idx, is_dom], "is_dom")
  			namecol(score, score_info[idx, is_rec], "is_rec")
  			if (score_info[idx, weight != "m"]) { 
  				namecol(score, score_info[idx, weight], "weight")
  			}
  		}, error=function(e) {
  			score_info[idx, error := "Could not match provided --score-<X> columns to columns in score file"]
  		})
      if (score_info[idx, !is.na(error)]) {
        return(NULL)
      }
   
      # Record how many variants are in the score file, and how many are on this chromosome
      score_info[idx, n_var := score[,.N]]
      if ("chr" %in% names(score)) {
  			score = score[chr == chrIdx] # filter to this chromosome. 
      }
  
  		# If multiple scores, transform to long format
  		if (score_info[idx, weight == "m"]) {
  			# first, make sure score name columns are unique, generate warning if necessary
  			sn = data.table(name = names(score))
  			dups = sn[,.N,by=name][N > 1, name]
  			sn[name %chin% dups, new := paste0(name, ".", seq_len(.N)), by=name]
  			sn[is.na(new), new := name]
  			if (sn[, any(name != new)]) {
  				warning("Non-unique score weight column names found in multi-score file", score_info[idx, path], 
  								". Numbers have been appended to make each score's name unique")
  			}
   
        # Give each sub-score a 'compName' 
        sn[!(name %chin% c("rsid", "chr", "pos", "EA", "EAF", "OA", "is_dom", "is_rec")),
             compName := gsub(" ", "0", paste0(score_info[idx, compName], ".", format(.I)))]
        sn[is.na(compName), compName := name]
        setnames(score, sn[,compName])
  			
        # Melt to long format
  			score <- melt(score, id.vars=intersect(names(score), c("rsid", "chr", "pos", "EA", "EAF", "OA", "is_dom", "is_rec")), 
  										variable.name="compName", value.name="weight")
  
        # create score sub-table to append to score-info
        sub_info <- score_info[idx]
        sub_info[, compName := NULL]
        sub_info <- cbind(sn[!(name %chin% c("rsid", "chr", "pos", "EA", "EAF", "OA", "is_dom", "is_rec")), .(compName)], sub_info)
        sub_info[sn, on = .(compName), weight := new]
        score_info <- rbindf(score_info, sub_info)
        
        # Temporarily drop zero-weights to save on memory before we later cast back to wide format.
        if (score[, any(is.na(weight))]) {
          warning("Missing values found in one or more weight columns and converted to 0 for multi-score file ", score_info[idx, path])
          score[is.na(weight), weight := 0]
        }
        score <- score[weight != 0] 
  		}
  
  		# If we're working with a single score drop extra columns that might be in the file and add the compName
  		if (score_info[idx, weight != "m"]) {
  			score <- score[, which(names(score) %chin% c("rsid", "chr", "pos", "EA", "EAF", "OA", "is_dom", "is_rec", "weight")), with=FALSE]
        if (score[, any(is.na(weight))]) {
          warning("Missing values found in weight column and converted to 0 for score file ", score_info[idx, path])
          score[is.na(weight), weight := 0]
        }
        score[, compName := score_info[idx, compName]]
  		}
      return(score)
  	} else {
  		# For PGS catalog files, we need to map score names to ones provided standard in the catalog.
      mapname <- function(old, new) {
        if (old %in% names(score)) {
          setnames(score, old, new)
          if (!is.character(score_info[[new]])) {
            score_info[, c(new) := as.character(score_info[[new]])]
          }
          score_info[idx, c(new) := old]
        } else {
          score_info[idx, c(new) := NA]
        }
      }
      tryCatch({
  			mapname("rsID", "rsid")
  			mapname("chr_name", "chr")
  			mapname("chr_position", "pos")
  			mapname("effect_allele", "EA")
  			mapname("allelefrequency_effect", "EAF")
  			mapname("effect_weight", "weight")
  			mapname("is_dominant", "is_dom")
  			mapname("is_recessive", "is_rec")
  
        # Old vs. new score file format has different names
        if ("reference_allele" %in% names(score)) {
  				mapname("reference_allele", "OA")
        } else {
  				mapname("other_allele", "OA")
        }
      }, error=function(e) {
        score_info[idx, error := "PGS Catalog file detected, but could not map column names"]
      })
      if (score_info[idx, !is.na(error)]) {
        return(NULL)
      }
  
      # Record how many variants are in the score file
      score_info[idx, n_var := score[,.N]]
      if ("chr" %in% names(score)) {
  			score = score[chr == chrIdx] # filter to this chromosome. 
      }
  
  		# Drop interaction terms
  		if ("is_interaction" %in% names(score)) {
  			int_var <- score[(is_interaction), .N,]
  			score <- score[!(is_interaction)]
  			score[, is_interaction := NULL]
  			score_info[idx, n_interaction_skipped := int_var]
  		}
  
      # Error if any variants have the same effect weight for both dominant and recessive effects
      if ("is_dom" %in% names(score) && "is_rec" %in% names(score)) {
        if (nrow(score[(is_dom) & (is_rec)]) > 0) {
          score_info[idx, error := "Entries with TRUE in both --score-dominant and --score-recessive columns"]
          return(NULL)
        }
      }
  
      # Add compName
      score[, compName := score_info[idx, compName]]
  
  		# Drop any columns we won't use
  	  score <- score[, which(names(score) %chin% c("rsid", "chr", "pos", "EA", "EAF", "OA", "is_dom", "is_rec", "weight", "compName")), with=FALSE]
  
      # Detect and convert missing values
  		if (score[, any(is.na(weight))]) {
  			warning("Missing values found in weight column and converted to 0 for score file ", score_info[idx, path])
  			score[is.na(weight), weight := 0]
  		}
      return(score)
  	}
  }
  rm(score)
  invisible(gc())
  
  # At this point, if all scores have errored, we can stop all tasks.
  if (all(!is.na(score_info$error))) {
    fwrite(score_info, sep="\t", quote=FALSE, file="errors/score_summary.txt")
    dx_upload("errors", args[["--work"]], exists="skip")
    stop("All score files had errors, see ", args[["--work"]], "errors/score_summary.txt")
  }
  
  # make sure chromosome column is character:
  if ("chr" %in% names(scores) && !is.character(scores$chr)) {
    scores[, chr := as.character(chr)]
  }
  
  # Fix any scores that use lower case alleles (WHY?)
  scores[, EA := toupper(EA)]
  if ("OA" %in% names(scores)) {
    scores[!is.na(OA), OA := toupper(OA)]
  }
  
  # Reorganise rows
  score_info = score_info[order(as.numeric(compName))]
      
  # There may be no variants on this chromosome 
  if (nrow(scores) == 0) {
    fwrite(score_info, sep="\t", quote=FALSE, file=sprintf("finished/score_summary_%s.txt", chrIdx))
    dx_upload("finished", args[["--work"]])
    quit(save="no")
  }
  
  # Checkpoint - first point where it makes sense to save intermediate work, i.e.
  # if restarted after this point we can skip downloading and processing the score
  # files
  fwrite(scores, sep="\t", file=sprintf("checkpointing/checkpoint1/scores_chr%s.txt", chrIdx))
  fwrite(score_info, sep="\t", file=sprintf("checkpointing/checkpoint1/score_info_chr%s.txt", chrIdx))
  dx_upload("checkpointing/checkpoint1", args[["--work"]])
} 

# If we have already computed all the information needed about the score files
# in an early job that was restarted, load that information so we can skip 
# recomputing it.
if (checkpoint == 1) {
  dir.create("checkpointing/checkpoint1")
  dx_download(sprintf("%s/checkpoint1/score_info_chr%s.txt", args[["--work"]], chrIdx), "checkpointing/checkpoint1/")
  dx_download(sprintf("%s/checkpoint1/scores_chr%s.txt", args[["--work"]], chrIdx), "checkpointing/checkpoint1/")
  score_info <- fread(sprintf("checkpointing/checkpoint1/score_info_chr%s.txt", chrIdx), na.strings = c("", "NA"), colClasses = c("compName"="character", "error"="character"))
  scores <- fread(sprintf("checkpointing/checkpoint1/scores_chr%s.txt", chrIdx), na.strings = c("", "NA"), colClasses = c("compName"="character", "chr"="character"))
}

# Match the score files to the variant information in the genotype data
if (checkpoint < 2) {
  dir.create("checkpointing/checkpoint2")
  
  ################################################################################################
  # Now each chromosome attempts to load the genotype variant information and orient variants 
  # to the same allele in UK Biobank. Also at this point we want to filter out ambiguous alleles
  # unless otherwise specified. When matching variants, default to matching by chromosome, 
  # position, and alleles unless otherwise asked (i.e. --match-by-rsid, or score missing
  # chr and pos columns).
  ################################################################################################

  # Function for flipping the strand of an allele.
  # Uses a series of gsub calls to replace A's with T's,
  # G's with C's, and vice-versa. Also works for alleles
  # with more than one nucleotide (e.g. indels).
  flip_strand <- function(x) {
    # Swap each letter for a dummy, we need this intermediate
    # step so we can distinguish between alleles when swapping.
    # E.g if we did A -> T then T -> A we'd end up with all A's
    # and no T's. instead we do A -> V -> T and T -> X -> A.
    x <- gsub("A", "V", x)
    x <- gsub("T", "X", x)
    x <- gsub("C", "Y", x)
    x <- gsub("G", "Z", x)
    x <- gsub("V", "T", x)
    x <- gsub("X", "A", x)
    x <- gsub("Y", "G", x)
    x <- gsub("Z", "C", x)
    return(x)
  }
  
  # Download the varfile (we have already checked it exists)
  varfile <- dx_download(varfile, sprintf("input_data/chr%s.%s", chrIdx, ifelse(args[["--genotype-format"]] == "pfile", "pvar", "bim")))
  
  # Attempt to load the variant information for this chromosome, exiting if not found
  tryCatch({
    # Detect if pgen from vcf
    line1 <- readLines(varfile, 1)
    if (grepl("^##", line1)) {
      varinfo <- fread(cmd=sprintf("grep -v '^#' %s", varfile), sep="\t")
    } else {
  		varinfo <- fread(varfile)
    }
  	if(args[["--genotype-format"]] == "pfile") {
      varinfo <- varinfo[,1:5] # drop extra columns from VCF pvar files
  		setnames(varinfo, c("chr", "pos", "rsid", "ref", "alt"))
    } else if (args[["--genotype-format"]] == "bfile") {
  		setnames(varinfo, c("chr", "rsid", "cm", "pos", "alt", "ref"))
    }
    if (!is.character(varinfo$chr)) {
      varinfo[,chr := as.character(chr)]
    }
  }, error=function(e) {
    score_info[is.na(error) & compName %in% unique(scores[!is.na(chr), compName]), c("n_match", "chr_fail") := .(0, TRUE)]
    fwrite(score_info, sep="\t", quote=FALSE, file=sprintf("errors/score_summary_%s.txt", chrIdx))
    fwrite(score_info, sep="\t", quote=FALSE, file=sprintf("finished/score_summary_%s.txt", chrIdx))
    dx_upload("errors", args[["--work"]])
    dx_upload("finished", args[["--work"]])
    quit(save="no")
  })
  
  # Make sure each variant has a unique identifier, and save for downstream
  # processing
  if (args[["--genotype-format"]] == "pfile") {
    fwrite(varinfo[, .(`#CHROM`=chr, POS=pos, ID=paste(chr, pos, alt, ref, sep=":"), REF=ref, ALT=alt)],
           sep="\t", quote=FALSE, file=varfile)
  } else if (args[["--genotype-format"]] == "bfile") {
    fwrite(varinfo[, .(chr, rsid=paste(chr, pos, alt, ref, sep=":"), cm, pos, alt, ref)],
           sep="\t", quote=FALSE, col.names=FALSE, file=varfile)
  }
  system(sprintf("cp %s checkpointing/checkpoint2/", varfile))
  
  if (args[["--single-geno"]]) {
  	varinfo <- varinfo[chr == chrIdx]
  	if (nrow(varinfo[chr == chrIdx]) == 0) {
      score_info[compName %in% unique(scores$compName), c("n_match", "chr_fail") := .(0, TRUE)]
  	  fwrite(score_info, sep="\t", quote=FALSE, file=sprintf("errors/score_summary_%s.txt", chrIdx))
  	  fwrite(score_info, sep="\t", quote=FALSE, file=sprintf("finished/score_summary_%s.txt", chrIdx))
  	  dx_upload("errors", args[["--work"]])
  	  dx_upload("finished", args[["--work"]])
  	  quit(save="no")
  	}
  }

  # How many variants on this chromosome for each score?
  if ("chr" %in% names(scores)) {
  	score_info[is.na(error) & !is.na(chr), n_chr := 0]
  	score_info[scores[!is.na(chr), .N, by=compName], on = .(compName), n_chr := i.N]
  	if (!("chr" %in% names(scores)) && !("pos" %in% names(scores))) {
  		score_info[scores[!is.na(chr), .N, by=.(compName=as.character(as.integer(as.vector(compName))), rsid, EA)][, .N, by=compName], on = .(compName), n_chr := i.N] # overall summary for multi-score file
  	} else {
  		score_info[scores[!is.na(chr), .N, by=.(compName=as.character(as.integer(as.vector(compName))), chr, pos, EA)][, .N, by=compName], on = .(compName), n_chr := i.N] # overall summary for multi-score file
  	}
  } else {
    score_info[, chr := NA_integer_]
  }

  # Obtain chromosome and position for variants missing this information,
  # (or override if explictly asked to match by rsid)
  # then subset again to just variants on this task's chromosome
  if ("rsid" %in% names(scores)) {
    scores[varinfo[!is.na(rsid)], on = .(rsid), c("chr", "pos") :=
           .(ifelse((args[["--match-by-rsid"]] & !is.na(rsid)) | !("chr" %in% names(scores)) | is.na(chr), i.chr, chr),
             ifelse((args[["--match-by-rsid"]] & !is.na(rsid)) | !("pos" %in% names(scores)) | is.na(pos), i.pos, pos))]
  }
  scores <- scores[chr == chrIdx]

  # can exit if no variants on this chromosome
  if (nrow(scores) == 0) {
    fwrite(score_info, sep="\t", quote=FALSE, file=sprintf("finished/score_summary_%s.txt", chrIdx))
    dx_upload("finished", args[["--work"]])
    quit(save="no")
  }

  # Before matching variants by alleles, check for field separators in
  # the other allele column. The matching process won't work with these,
  # so we can just set those other allele fields to NA to match by effect
  # allele only
  if ("OA" %in% names(scores)) {
  	scores[!grepl("^(A|G|C|T)+$", OA), OA := NA]
  }

  # Check if we can match by chromosome, position, and alleles
  scores[, match := FALSE]
  if ("OA" %in% names(scores)) {
    # Where both effect and other allele provided, must match entry by both
    scores[varinfo, on = .(chr, pos, EA=alt, OA=ref), match := !is.na(OA)]
    scores[varinfo, on = .(chr, pos, EA=ref, OA=alt), match := !is.na(OA)]
    # If no other allele provided, match just by effect allele
    scores[varinfo, on = .(chr, pos, EA=alt), match := match | is.na(OA)]
    scores[varinfo, on = .(chr, pos, EA=ref), match := match | is.na(OA)]
  } else {
    scores[varinfo, on = .(chr, pos, EA=alt), match := TRUE]
    scores[varinfo, on = .(chr, pos, EA=ref), match := TRUE]
  }

  # For variants that did not match above, but can match by chromosome and position,
  # flip the strand of the alleles and attempt match again
  if ("OA" %in% names(scores)) {
    # Flip alleles
    scores[varinfo, on = .(chr, pos), c("EA", "OA") := 
           .(ifelse(match, EA, flip_strand(EA)), ifelse(match, OA, flip_strand(OA)))]
    # Where both effect and other allele provided, must match entry by both
    scores[varinfo, on = .(chr, pos, EA=alt, OA=ref), match := ifelse(match, match, !is.na(OA))]
    scores[varinfo, on = .(chr, pos, EA=ref, OA=alt), match := ifelse(match, match, !is.na(OA))]
    # If no other allele provided, match just by effect allele
    scores[varinfo, on = .(chr, pos, EA=alt), match := ifelse(match, match, match | is.na(OA))]
    scores[varinfo, on = .(chr, pos, EA=ref), match := ifelse(match, match, match | is.na(OA))]
  } else {
    scores[varinfo, on = .(chr, pos), c("EA") := .(ifelse(match, EA, flip_strand(EA)))]
    scores[varinfo, on = .(chr, pos, EA=alt), match := ifelse(match, match, TRUE)]
    scores[varinfo, on = .(chr, pos, EA=ref), match := ifelse(match, match, TRUE)]
  }

  # Drop score variants that could not be matched to the genotype data
  scores <- scores[(match)]
  scores[, match := NULL]
  invisible(gc())

  # How many variants could we match to the genotype data?
  score_info[is.na(error), n_match := 0]
  score_info[scores[, .N, by=compName], on = .(compName), n_match := i.N]
  score_info[scores[, .N, by=.(compName=as.character(as.integer(as.vector(compName))), chr, pos)][, .N, by=compName], on = .(compName), n_match := i.N] # overall summary for multi-score file

  # can exit if no variants matched
  if (nrow(scores) == 0) {
    fwrite(score_info, sep="\t", quote=FALSE, file=sprintf("finished/score_summary_%s.txt", chrIdx))
    dx_upload("finished", args[["--work"]])
    quit(save="no")
  }

  # Identify and flag variants that are multi-allelic in either the score file or in the genotype data
  genomulti <- varinfo[,.N,by=pos][N > 1, .(pos)]
  multi <- unique(rbind(
    scores[genomulti, on=.(pos), nomatch=0, .(pos, compName)],
    scores[,.N,by=.(pos,compName)][N > 1, .(pos, compName)]))

  # Collate number of multi-allelic variants
  score_info[is.na(error), n_multiallele := 0]
  score_info[multi[,.N, by=compName], on = .(compName), n_multiallele := i.N]
  score_info[multi[, .N, by=.(compName=as.character(as.integer(as.vector(compName))), pos)][, .N, by=compName], on = .(compName), n_multiallele := i.N] # overall summary for multi-score file

  # Remove them if requested
  if (args[["--remove-multiallelic"]]) {
    scores <- scores[!multi, on = .(pos)]
    score_info[, n_multiallele_removed := n_multiallele]
  } else {
    score_info[, n_multiallele_removed := ifelse(n_multiallele > 0, 0, NA)]
  }

  # Check at this point whether each score has only 1 effect weight per variant.
  bad <- unique(scores[,.N,by=.(pos, compName)][N > 1, .(compName)])
  score_info[bad, on = .(compName), error := "Score has multiple effect weights for the same variant/position"]
  scores <- scores[!bad, on = .(compName)]

  # can exit if all errors and no variants remain
  if (nrow(scores) == 0) {
    fwrite(score_info, sep="\t", quote=FALSE, file=sprintf("errors/score_summary_%s.txt", chrIdx))
    fwrite(score_info, sep="\t", quote=FALSE, file=sprintf("finished/score_summary_%s.txt", chrIdx))
    dx_upload("errors", args[["--work"]])
    dx_upload("finished", args[["--work"]])
    quit(save="no")
  }

  # Now that we've done the matching, give the variants in the score file
  # the same unique identifier they have in the variant information file.
  # For multi-allelic variants we want to avoid double counting alleles,
  # e.g. where the effect allele is the one in common across the multiple
  # rows.
  scores[varinfo, on = .(pos, EA=alt), rsid := paste(i.chr, i.pos, i.alt, i.ref, sep=":")]
  scores[varinfo, on = .(pos, EA=ref), rsid := paste(i.chr, i.pos, i.alt, i.ref, sep=":")]

  # Now identify variants whose effect allele is ambiguous (e.g. A/T or G/C SNPs). The 
  # following works for biallelic sites and some multi-allelic sites
  scores[varinfo, on = .(pos, EA=alt), ambig := EA == flip_strand(ref)]
  scores[varinfo, on = .(pos, EA=ref), ambig := EA == flip_strand(alt)]

  # For multi-allelic sites that otherwise weren't identified as ambiguous:
  # 1. if the OA column is provided, check if the effect allele is strand ambiguous
  #    with any other allele at that position.
  # 2. Remove the OA column.
  # 3. Add in all OA from the varinfo table, then check again
  mult_ambig <- scores[multi, on = .(pos, compName)][!(ambig)]
  if ("OA" %in% names(mult_ambig)) {
    mult_ambig[, ambig := any(EA == flip_strand(OA)), by=pos]
    scores[mult_ambig[(ambig)], on = .(pos, EA), ambig := TRUE]
    mult_ambig <- mult_ambig[!(ambig)]
    mult_ambig[, OA := NULL]
  }
  mult_ambig <- rbind(
    varinfo[, .(pos, EA=alt, OA=ref)][mult_ambig, on = .(pos, EA), nomatch=0],
    varinfo[, .(pos, EA=ref, OA=alt)][mult_ambig, on = .(pos, EA), nomatch=0])
  mult_ambig[, ambig := any(EA == flip_strand(OA)), by=pos]
  mult_ambig <- unique(mult_ambig[,.(pos, EA, ambig)])
  scores[mult_ambig[(ambig)], on = .(pos, EA), ambig := TRUE]
  rm(mult_ambig)

  # Extract ambiguous variants so we can do further checks by allele frequency if requested
  ambig <- scores[(ambig)]
  scores <- scores[!(ambig)]
  scores[, ambig := NULL]
  ambig[, ambig := NULL]
  invisible(gc())

  # How many ambiguous alleles in each score?
  score_info[is.na(error), n_ambig := 0]
  if (nrow(ambig) > 0) {
    score_info[ambig[,.N,by=compName], on = .(compName), n_ambig := i.N]
    score_info[ambig[,.N, by=.(compName=as.character(as.integer(as.vector(compName))), chr, pos)][, .N, by=compName], on = .(compName), n_ambig := i.N] # overall summary for multi-score file
  }

  # If the freqx file has been provided (either for allele frequency checking of ambiguous SNPs,
  # or for providing directly to plink2 to skip allele frequency calculations) we need to add 
  # the missing position information and give it the same unique variant identifiers as the
  # variant info file.
  tryCatch({
    if (!is.null(freqxfile)) {
      freqxfile <- dx_download(freqxfile, "input_data/")
  	  freqx <- fread(freqxfile)
    }
  }, error=function(e) {
    warning("Error when attempting to read freqx file: ", freqxfile, " computing from genotype data instead")
  })
  
  # Parse the freqx file. Depends on the version of plink and format requested:
  if (exists("freqx")) {
    if ("C(HOM A1)" %in% names(freqx)) { # Plink 1.9 --freqx
      if (inherits(freqx$CHR, "integer")) freqx[, CHR := as.character(CHR)]
      freqx[varinfo, on = .(SNP=rsid, CHR=chr, A2=ref, A1=alt), POS := i.pos]
      freqx[, SNP := paste(CHR, POS, A1, A2, sep=":")]
      freqx[, POS := NULL]
      fwrite(freqx, sep="\t", quote=FALSE, file=sprintf("checkpointing/checkpoint2/chr%s.frqx", chrIdx))
      freqx_ext <- "frqx"
    } else if ("MAF" %in% names(freqx)) { # Plink 1.9 --freq
      if (inherits(freqx$CHR, "integer")) freqx[, CHR := as.character(CHR)]
      freqx[varinfo, on = .(SNP=rsid, CHR=chr, A2=ref, A1=alt), POS := i.pos]
      freqx[, SNP := paste(CHR, POS, A1, A2, sep=":")]
      freqx[, POS := NULL]
      fwrite(freqx, sep="\t", quote=FALSE, file=sprintf("checkpointing/checkpoint2/chr%s.frq", chrIdx))
      freqx_ext <- "frq"
    } else if ("ALT" %in% names(freqx)) { # Plink 2 --freq
      setnames(freqx, "#CHROM", "CHROM")
      if (inherits(freqx$CHROM, "integer")) freqx[, CHROM := as.character(CHROM)]
      freqx[varinfo, on = .(ID=rsid, CHROM=chr, ALT=alt, REF=ref), POS := i.pos]
      freqx[, ID := paste(CHROM, POS, ALT, REF, sep=":")]
      freqx[, POS := NULL]
      fwrite(freqx, sep="\t", quote=FALSE, file=sprintf("checkpointing/checkpoint2/chr%s.afreq", chrIdx))
      freqx_ext <- "afreq"
    } else {
      rm(freqx)
      warning("Unrecognised freqx format, computing from genotype data instead")
    }
  }

  # Now we can give the loaded varinfo table the appropriate unique identifiers
  # and subset to just the score variants
  varinfo[, rsid := paste(chr, pos, alt, ref, sep=":")]
  varinfo <- varinfo[pos %in% unique(c(scores$pos, ambig$pos))] # we want to preserve row order, so filter instead of join.

  # Handle ambiguous variants depending of various input options and collate
  if (nrow(ambig) > 0) {
  	if (!args[["--keep-ambiguous"]]) {
  		ambig[, keep := FALSE]
  	} else if (args[["--keep-ambiguous"]] && is.null(args[["--ambiguous-thresh"]])) {
  		ambig[, keep := TRUE]
  	} else if (!("EAF" %in% names(scores))) {
  		# Requested to match by effect allele frequency, but no scores provided this.
  		ambig[, keep := FALSE] 
  	} else { # keep if can match by effect allele frequency
  		if (!exists("freqx")) {
  			# Need to make sure we give variants unique identifiers when writing out
  			fwrite(ambig[,.(unique(rsid))], col.names=FALSE, quote=FALSE,
  						 file=sprintf("work/ambig_freqx_extract_chr%s", chrIdx))
  		  
  		  # Download genotype data
  		  dx_download(genofile, sprintf("input_data/chr%s.%s", chrIdx, ifelse(args[["--genotype-format"]] == "pfile", "pgen", "bed")))
  		  dx_download(samplefile, sprintf("input_data/chr%s.%s", chrIdx, ifelse(args[["--genotype-format"]] == "pfile", "psam", "fam")))
  
        # Do frequency calculation for ambiguous variants
  			cmd <- normalizePath(sprintf("%s/plink2", SoftwareDir))
  			cmd[2] <- sprintf("--%s input_data/chr%s", args[["--genotype-format"]], chrIdx)
  			cmd[3] <- sprintf("--extract work/ambig_freqx_extract_chr%s", chrIdx)
  			cmd[4] <- sprintf("--freq --out checkpointing/checkpoint2/ambig_freqx_extract_chr%s", chrIdx)
  			if (!is.null(args[["--keep"]])) cmd[6] <- paste("--keep", args[["--keep"]])
  			errorcode <- system(paste(cmd, collapse=" "), wait=TRUE)
  
  			if (errcode != 0) {
  				score_info[, error := sprintf("%sError from plink2 when computing frequencies of ambiguous alleles, see log file: %s", 
  					ifelse(is.na(error) | error == "", "", paste0(error, ".")), 
  					sprintf("%s/errors/ambig_freqx_extract_chr%s.log", args[["--work"]], chrIdx))]
  			  system("mv checkpointing/checkpoint2/*.log errors/")
  				fwrite(score_info, sep="\t", quote=FALSE, file=sprintf("errors/score_summary_%s.txt", chrIdx))
  				fwrite(score_info, sep="\t", quote=FALSE, file=sprintf("finished/score_summary_%s.txt", chrIdx))
  				dx_upload("errors/", args[["--work"]])
  				dx_upload("finished/", args[["--work"]])
  				quit(save="no")
  			}
  
        # load allele frequencies
  			freqx <- fread(sprintf("checkpointing/checkpoint2/ambig_freqx_extract_chr%s.afreq", work_dir, chrIdx))
  		  freqx <- freqx[, .(rsid=ID, EA=ALT, EAF=ALT_FREQS)]
  		} else {
        # Need to extract EAF from loaded freqx file
        if (freqx_ext == "frqx") {
          freqx <- freqx[, .(rsid=SNP, EA=A1, EAF=(`C(HOM A1)`*2+`C(HET)`)/(`C(HOM A1)`*2+`C(HET)`+`C(HOM A2)`*2))]
        } else if (freqx_ext == "frq") {
          freqx <- freqx[, .(rsid=SNP, EA=A1, EAF=MAF)]
        } else if (freqx_ext == "afreq") {
          freqx <- freqx[, .(rsid=ID, EA=ALT, EAF=ALT_FREQS)]
        }
      }
  
      # Flag multi-allelic SNPs in freqx file - for these we can only match the the allele whose frequency has been calculated
      freqx[, pos := as.numeric(gsub(":.*", "", gsub("^.*?:", "", rsid)))]
      freqx[varinfo[,.N,by=pos], on=.(pos), multi := N > 1]
  
  		# Where EAF column not provided, set keep to FALSE
  		ambig[, keep := TRUE]
  		ambig[is.na(EAF), keep := FALSE]
  
  		# Attempt to match ambiguous variants by EAF below threshold:
      ambig[freqx, on = .(rsid, EA), # effect allele is the one we have frequency for in freqx.
            keep := ifelse(keep &    # ignore score variants with no EAF information
              ((EAF < args[["--ambiguous-thresh"]] & i.EAF < args[["--ambiguous-thresh"]]) |    # effect allele is minor allele, and frequency below threshold
               (EAF > (1 - args[["--ambiguous-thresh"]]) & i.EAF > (1 - args[["--ambiguous-thresh"]]))), # effect allele is major allele, and frequency below threshold
              TRUE, FALSE)]
      # The same as above, but instead the effect allele is the allele whose frequency was not calculated.
      # In this case we can take the effect allele's EAF in the genotype data as 1 - the calculated EAF,
      # provided we're dealing only with bi-allelic SNPs.
      ambig[freqx[!(multi)], on = .(rsid, OA=EA),
            keep := ifelse(keep &   
              ((EAF < args[["--ambiguous-thresh"]] & (1 - i.EAF) < args[["--ambiguous-thresh"]]) |    # effect allele is minor allele, and frequency below threshold
               (EAF > (1 - args[["--ambiguous-thresh"]]) & (1 - i.EAF) > (1 - args[["--ambiguous-thresh"]]))), # effect allele is major allele, and frequency below threshold
              TRUE, FALSE)]
      
      # Remove no longer needed freqx object
      rm(freqx)
  	}
    
    # Filter to just ambiguous variants we want to keep
    ambig <- ambig[(keep)]
    ambig[, keep := NULL]
    invisible(gc())
  }

  # How many ambiguous alleles are kept?
  score_info[is.na(error), n_ambig_kept := 0]
  if (nrow(ambig) > 0) {
  	score_info[ambig[,.N,by=compName], on = .(compName), n_ambig_kept := i.N]
  	score_info[ambig[,.N, by=.(compName=as.character(as.integer(as.vector(compName))), chr, pos)][, .N, by=compName], on = .(compName), n_ambig_kept := i.N] # overall summary for multi-score file
  }

  # Add kept variants back to score and collate
  scores <- rbind(scores, ambig)
  rm(ambig)
  invisible(gc())
  score_info[, n_used := n_match - n_ambig + n_ambig_kept]

  # can exit if no variants left (i.e. all matched were ambiguous)
  if (nrow(scores) == 0) {
    fwrite(score_info, sep="\t", quote=FALSE, file=sprintf("finished/score_summary_%s.txt", chrIdx))
    dx_upload("finished", args[["--work"]])
    quit(save="no")
  }

  # Split scores into those for linear, dominant, and recessive effects,
  # transform to wide format, write out, and remove.
  sdcast <- function(...) {
    tryCatch({
      dcast(...)
    }, warning=function(w) {
      score_info[!is.na(error), error := paste0("Critical error collating score weights on chromosome ", 
                                                chrIdx, ": encountered multiple weights for the same effect allele.")]
      fwrite(score_info, sep="\t", quote=FALSE, file=sprintf("errors/score_summary_%s.txt", chrIdx))
      fwrite(score_info, sep="\t", quote=FALSE, file=sprintf("finished/score_summary_%s.txt", chrIdx))
      dx_upload("errors", args[["--work"]])
      dx_upload("finished", args[["--work"]])
      quit(save="no")
    })
  }

  # If different scores use different effect alleles, we will need to create more than one file since
  # plink can handle only one entry per variant
  make_score_files <- function(dt, model_name, file_prefix) {
    var_ea <- unique(dt[,.(rsid, EA)])
    var_ea[, group := seq_len(.N), by=rsid]
    sfinfo <- data.table(model = model_name, group = unique(var_ea$group), ncol=0) 
    for (group_id in unique(var_ea$group)) {
      subdt <- dt[var_ea[group == group_id], on = .(rsid, EA)]
      wide <- sdcast(subdt, rsid + EA ~ compName, value.var="weight", fill=0)
      wide <- wide[varinfo[, .(rsid)], on = .(rsid), nomatch=0] # row order to match that of genotype data.
      fwrite(wide, sep="\t", quote=FALSE, file=sprintf("%s_group%s.txt", file_prefix, group_id))
      sfinfo[model == model_name & group == group_id, ncol := ncol(wide)]
    }
    return(sfinfo)
  }

  plink_input_info <- data.table()
  if ("is_dom" %in% names(scores)) {
    dominant <- scores[(is_dom)]
    scores <- scores[!(is_dom) | is.na(is_dom)]
  
    if (nrow(dominant) > 0) {
      plink_input_info <- rbind(fill=TRUE, plink_input_info, 
          make_score_files(dominant, "dominant", sprintf("checkpointing/checkpoint2/collated_scores_dominant_chr%s", chrIdx)))
    }
  
    rm(dominant)
  }

  if ("is_rec" %in% names(scores)) {
    recessive <- scores[(is_rec)]
    scores <- scores[!(is_rec) | is.na(is_rec)]
  
    if (nrow(recessive) > 0) {
      plink_input_info <- rbind(fill=TRUE, plink_input_info, 
          make_score_files(recessive, "recessive", sprintf("checkpointing/checkpoint2/collated_scores_recessive_chr%s", chrIdx)))
    }
  
    rm(recessive)
  }

  if (nrow(scores) > 0) {
    plink_input_info <- rbind(fill=TRUE, plink_input_info, 
        make_score_files(scores, "linear", sprintf("checkpointing/checkpoint2/collated_scores_linear_chr%s", chrIdx)))
    rm(scores)
  }

  # save score_info and remove unnecessary objects in memory while plink2 runs
  fwrite(plink_input_info, sep="\t", quote=FALSE, file=sprintf("checkpointing/checkpoint2/plink_input_info_%s.txt", chrIdx))
  fwrite(score_info, sep="\t", quote=FALSE, file=sprintf("checkpointing/checkpoint2/score_summary_%s.txt", chrIdx))
  dx_upload("checkpointing/checkpoint2", args[["--work"]])
  rm(varinfo)
  invisible(gc())
}

# If we have already done the score file matching to the genotype data download
# it so we can skip recomputing it
if (checkpoint == 2) {
  dir.create("checkpointing/checkpoint2")
  
  if (args[["--genotype-format"]] == "pfile") {
    dx_download(sprintf("%s/checkpoint2/chr%s.pvar", args[["--work"]], chrIdx), "input_data/")
  } else {
    dx_download(sprintf("%s/checkpoint2/chr%s.bim", args[["--work"]], chrIdx), "input_data/")
  }
  
  if (dx_exists(sprintf("%s/checkpoint2/chr%s.frqx", args[["--work"]], chrIdx))) {
    dx_download(sprintf("%s/checkpoint2/chr%s.frqx", args[["--work"]], chrIdx), "checkpointing/checkpoint2/")
    freqx_ext <- "frqx"
  } else if (dx_exists(sprintf("%s/checkpoint2/chr%s.frq", args[["--work"]], chrIdx))) {
    dx_download(sprintf("%s/checkpoint2/chr%s.frq", args[["--work"]], chrIdx), "checkpointing/checkpoint2/")
    freqx_ext <- "frq"
  } else if (dx_exists(sprintf("%s/checkpoint2/chr%s.afreq", args[["--work"]], chrIdx))) {
    dx_download(sprintf("%s/checkpoint2/chr%s.afreq", args[["--work"]], chrIdx), "checkpointing/checkpoint2/")
    freqx_ext <- "afreq"
  }
  
  dx_download(sprintf("%s/checkpoint2/score_summary_%s.txt", args[["--work"]], chrIdx), "checkpointing/checkpoint2/")
  score_info <- fread(sprintf("checkpointing/checkpoint2/score_summary_%s.txt", chrIdx), colClasses = c("compName"="character", "error"="character"))
  
  dx_download(sprintf("%s/checkpoint2/plink_input_info_%s.txt", args[["--work"]], chrIdx), "checkpointing/checkpoint2/")
  plink_input_info <- fread(sprintf("checkpointing/checkpoint2/plink_input_info_%s.txt", chrIdx))
  
  for (idx in plink_input_info[,.I]) {
    dx_download(sprintf("%s/checkpoint2/collated_scores_%s_chr%s_group%s.txt", 
      args[["--work"]], plink_input_info[idx, model], chrIdx, plink_input_info[idx, group]),
      "checkpointing/checkpoint2/")
  }
}

# Compute the scores
if (checkpoint < 3) {
  dir.create("checkpointing/checkpoint3")
  
  ##############################################################################################
  # Now that we have collated the score files, we can compute the score
  ##############################################################################################
  
  # For hard call genotype data, mean-imputation of missing alleles only works when > 50 samples
  no_mean_imp <- FALSE
  if (args[["--genotype-format"]] == "bfile") {
    n_samples <- as.numeric(system(sprintf("wc -l %s | cut -f 1 -d ' '", samplefile), intern=TRUE))
    if (n_samples < 50) {
  		warning("Less than 50 samples, no-mean-imputation flag to --score has been set.") 
  		no_mean_imp <- TRUE
    }
  }
  
  # Download genotype data to local storage
  if (args[["--genotype-format"]] == "pfile") {
    dx_download(genofile, sprintf("input_data/chr%s.pgen", chrIdx))
    dx_download(samplefile, sprintf("input_data/chr%s.psam", chrIdx))
  } else {
    dx_download(genofile, sprintf("input_data/chr%s.bed", chrIdx))
    dx_download(samplefile, sprintf("input_data/chr%s.fam", chrIdx))
  }
  
  # Check for --keep file if provided and download local copy
  if (!is.null(args[["--keep"]])) {
    args[["--keep"]] <- dx_download(args[["--keep"]], "input_data/")
  }

  # Calculate the score levels for each model. Also extract subcohort if asked.
  for (idx in plink_input_info[,.I]) {
    # Get set of variants to extract
    system(sprintf("tail -n +2 checkpointing/checkpoint2/collated_scores_%s_chr%s_group%s.txt | cut -f 1 > checkpointing/checkpoint3/collated_scores_%s_variants_chr%s_group%s.txt", 
           plink_input_info[idx, model], chrIdx, plink_input_info[idx, group],
           plink_input_info[idx, model], chrIdx, plink_input_info[idx, group]), wait=TRUE)
  
    # Calculate the score sums
    cmd <- "plink2 --silent"
    cmd[2] <- sprintf("--%s input_data/chr%s", args[["--genotype-format"]], chrIdx)
    cmd[3] <- sprintf("--out checkpointing/checkpoint3/collated_scores_%s_chr%s_group%s", plink_input_info[idx, model], chrIdx, plink_input_info[idx, group])
    cmd[4] <- sprintf("--extract checkpointing/checkpoint3/collated_scores_%s_variants_chr%s_group%s.txt",  plink_input_info[idx, model], chrIdx, plink_input_info[idx, group])
    cmd[5] <- sprintf("--score checkpointing/checkpoint2/collated_scores_%s_chr%s_group%s.txt",  plink_input_info[idx, model], chrIdx, plink_input_info[idx, group])
    cmd[7] <- sprintf("'header-read' 'ignore-dup-ids' 'cols=scoresums'")
    if (no_mean_imp) cmd[8] <- paste("no-mean-imputation")
    if (plink_input_info[idx, model] != "linear") cmd[9] <- sprintf("'%s'", plink_input_info[idx, model])
    if (plink_input_info[idx, ncol] > 3) cmd[10] <- sprintf("--score-col-nums 3-%s", plink_input_info[idx, ncol])
    if (exists("freqx_ext") && args[["--genotype-format"]] == "bfile") 
      cmd[11] <- sprintf("--read-freq %s/chr%s.%s", work_dir, chrIdx, freqx_ext)
    if (!is.null(args[["--keep"]])) cmd[12] <- paste("--keep", args[["--keep"]])
    errcode <- system(paste(na.omit(cmd), collapse=" "), wait=TRUE)
  
    if (errcode != 0) {
      logfile <- sprintf("collated_scores_%s_chr%s_group%s.log", plink_input_info[idx, model], chrIdx, plink_input_info[idx, group])
      score_info[, error := sprintf("%sError from plink2 when computing scores, see log file: %s", 
        ifelse(is.na(error) | error == "", "", paste0(error, ".")),
        sprintf("%s/errors/collated_scores_%s_chr%s_group%s.log", 
                args[["--work"]], plink_input_info[idx, model], chrIdx, plink_input_info[idx, group]))]
      system("mv checkpointing/checkpoint3/*.log errors/")
      fwrite(score_info, sep="\t", quote=FALSE, file=sprintf("errors/score_summary_%s.txt", chrIdx))
      fwrite(score_info, sep="\t", quote=FALSE, file=sprintf("finished/score_summary_%s.txt", chrIdx))
      dx_upload("errors/", args[["--work"]])
      dx_upload("finished/", args[["--work"]])
      quit(save="no")
    }
  }
  
  # Upload computed scores
  dx_upload("checkpointing/checkpoint3", args[["--work"]])

  # Create file to let others know this chromosome has finished processing
  system(sprintf("mv checkpointing/checkpoint2/score_summary_%s.txt finished/", chrIdx))
  dx_upload("finished", args[["--work"]])
}

########################################################################################
# Collate the results and move to target output directory(s)
########################################################################################

# Check if this is the last running chromosome, and only proceed in thats the case
for (ii in 1:taskMax) {
  ii <- task_to_chr(as.character(ii))
  if (!dx_exists(sprintf("%s/finished/score_summary_%s.txt", args[["--work"]], ii))) {
    quit(save="no") # other tasks still running, finish.
  }
}

# Download all the work in progress files from all chromosomes
dx_download(sprintf("%s/checkpoint2/", args[["--work"]]), "checkpointing/checkpoint2/")
dx_download(sprintf("%s/checkpoint3/", args[["--work"]]), "checkpointing/checkpoint3/")
dx_download(sprintf("%s/finished", args[["--work"]]))

# Collate the plink logs
system("touch output/collated_plink_logs.txt")

freqxlogs <- list.files(path="checkpointing/checkpoint2", pattern="ambig_freqx_extract_.*.log", full.names=TRUE)
for (ff in freqxlogs) {
  system(sprintf("cat %s >> output/collated_plink_logs.txt", ff), wait=TRUE)
}

scoringlogs <- list.files(path="checkpointing/checkpoint3", pattern="collated_scores_.*.log", full.names=TRUE)
for (ff in scoringlogs) {
  system(sprintf("cat %s >> output/collated_plink_logs.txt", ff), wait=TRUE)
}

# Collate a list of matched variants
rbindu <- function(...) { unique(rbind(...)) }
for(ii in 1:taskMax) {
  chr <- task_to_chr(as.character(ii))
  # Load subset of variants that were matched to any score
  varmatchfiles <- list.files(path="checkpointing/checkpoint3", pattern=sprintf("collated_scores_.*_variants_chr%s_.*.txt", chr), full.names=TRUE)
  if (length(varmatchfiles) == 0) next # nothing on this chromosome
  varmatch <- foreach(ff = varmatchfiles, .combine=rbindu) %do% {
	  fread(ff, header=FALSE)
  }
  setnames(varmatch, "match_id")

  # Load full variant information for this chromosome
  if (args[["--genotype-format"]] == "pfile") {
    varinfo <- fread(sprintf("checkpointing/checkpoint2/chr%s.pvar", chr))
    setnames(varinfo, c("chromosome", "position", "match_id", "ref_allele", "alt_allele"))
  } else if (args[["--genotype-format"]] == "bfile") {
    varinfo <- fread(sprintf("checkpointing/checkpoint2/chr%s.bim", chr))
    setnames(varinfo, c("chromosome", "match_id", "centimorgan", "position", "minor_allele", "major_allele"))
  }
  
  # Obtain the rsid in the original data
  origfile <- paste0(args[["--genotype-prefix"]], ifelse(args[["--single-geno"]], "", chr), args[["--genotype-suffix"]])  
  if (args[["--genotype-format"]] == "pfile") {
    origfile <- sprintf("%s.pvar", origfile)
    origfile <- dx_download(origfile, "input_data/")
    line1 <- readLines(origfile, 1)
		if (grepl("^##", line1)) {
			varinfo[, rsid := fread(cmd=sprintf("grep -v '^#' %s", origfile), sep="\t")[["V3"]]]
    } else {
			varinfo[, rsid := fread(cmd=sprintf("tail -n +2 %s | cut -f 3", origfile), header=FALSE)[[1]]]
    }
  } else if (args[["--genotype-format"]] == "bfile") {
    origfile <- sprintf("%s.bim", origfile)
    origfile <- dx_download(origfile, "input_data/")
    varinfo[, rsid := fread(cmd=sprintf("cut -f 2 %s.bim", origfile), header=FALSE)[[1]]]
  }

  # Obtain the effect allele frequency, if applicable
  if (file.exists(sprintf("checkpointing/checkpoint2/ambig_freqx_extract_chr%s.afreq", chr))) {
    freqx <- fread(sprintf("checkpointing/checkpoint2/ambig_freqx_extract_chr%s.afreq", chr))
    freqx <- freqx[, .(match_id=ID, alt_frequency=ALT_FREQS)]
    varinfo[freqx, on = .(match_id), alt_frequency := i.alt_frequency]
  } else if (file.exists(sprintf("checkpointing/checkpoint2/chr%s.frqx", chr))) {
    freqx <- fread(file.exists(sprintf("checkpointing/checkpoint2/chr%s.frqx", chr)))
		freqx <- freqx[, .(match_id=SNP, MAF=(`C(HOM A1)`*2+`C(HET)`)/(`C(HOM A1)`*2+`C(HET)`+`C(HOM A2)`*2))]
    varinfo[freqx, on = .(match_id), MAF := i.MAF]
  } else if (file.exists(sprintf("checkpointing/checkpoint2/chr%s.frq", chr))) {
    freqx <- fread(file.exists(sprintf("checkpointing/checkpoint2/chr%s.frq", chr))) 
		freqx <- freqx[, .(match_id=SNP, MAF)]
    varinfo[freqx, on = .(match_id), MAF := i.MAF]
  } else if (file.exists(sprintf("checkpointing/checkpoint2/chr%s.afreq", chr))) { 
    freqx <- fread(file.exists(sprintf("checkpointing/checkpoint2/chr%s.afreq", chr)))
		freqx <- freqx[, .(match_id=ID, alt_frequency=ALT_FREQS)]
    varinfo[freqx, on = .(match_id), alt_frequency := i.alt_frequency]
  }
  
  # Extract the subset of matched variants
  varinfo <- varinfo[varmatch, on = .(match_id)]
  
  # reorganise columns (not all will be present, hence intersect())
  varinfo <- varinfo[, intersect(c("match_id", "chromosome", "position", 
     "alt_allele", "ref_allele", "minor_allele", "major_allele", "rsid", 
     "alt_frequency", "MAF"), 
     names(varinfo)), with=FALSE]

  # write out
  if (!file.exists("output/matched_variants.txt")) {
    cat("# This file contains variants that matched the genotype data across all input score files\n", file="output/matched_variants.txt")
    cat("# which may include scores whose levels are not saved in this folder (see run_score_file_list.txt)\n", file="output/matched_variants.txt", append=TRUE)
    cat("# Note different input scores may use either pair of alleles as the effect allele,\n", file="output/matched_variants.txt", append=TRUE)
    cat("# and may be oriented to the opposite strand.\n", file="output/matched_variants.txt", append=TRUE)
    fwrite(varinfo[0], sep="\t", quote=FALSE, file="output/header.txt")
    system("cat output/header.txt >> output/matched_variants.txt", wait=TRUE)
    system("rm -f output/header.txt", wait=TRUE)
  }
  fwrite(varinfo, sep="\t", quote=FALSE, append=TRUE, file="output/matched_variants.txt")
}
system("gzip output/matched_variants.txt", wait=TRUE)

# Load all the score sum files and total across chromosomes and models. We need to 
# progressively load and sum to avoid having to load all scores across all chromosomes
# into memory at once. 
sscorefiles <- list.files(path="checkpointing/checkpoint3", pattern="*.sscore", full.names=TRUE)
sscores <- fread(sscorefiles[1])
setnames(sscores, names(sscores)[1], "IID")
setnames(sscores, names(sscores)[-1], gsub("_SUM", "", names(sscores)[-1]))
cn <- names(sscores)[-1] # keep track of scores loaded so far (may not be all if score has no variants on this chromosome)
integer_scores <- setdiff(which(sapply(sscores, class) == "integer"), 1L)
for (ic in integer_scores) sscores[[ic]] <- as.numeric(sscores[[ic]])
sscores <- melt(sscores, id.vars=c("IID"), variable.name="compName", value.name="score_sum")
if (length(sscorefiles) > 1) {
  for (sidx in 2:length(sscorefiles)) {
    sscores2 <- fread(sscorefiles[sidx])
    setnames(sscores2, names(sscores2)[1], "IID")
    setnames(sscores2, names(sscores2)[-1], gsub("_SUM", "", names(sscores2)[-1]))
    cn_new <- setdiff(names(sscores2)[-1], cn)
		integer_scores <- setdiff(which(sapply(sscores2, class) == "integer"), 1L)
		for (ic in integer_scores) sscores2[[ic]] <- as.numeric(sscores2[[ic]])
    sscores2 <- melt(sscores2, id.vars=c("IID"), variable.name="compName", value.name="score_sum")
    sscores[sscores2, on = .(compName, IID), score_sum := .(score_sum + i.score_sum)]
    sscores <- rbind(sscores, sscores2[compName %in% cn_new]) # need to add in scores that did not have variants on any previous chromosome
    cn <- c(cn, cn_new)
    rm(sscores2)
    invisible(gc())
  }
}

# Load all score_info files and collate
sinfofiles <- list.files(path="finished", pattern="score_summary_*")
score_info <- lapply(sprintf("finished/%s", sinfofiles), fread, colClasses=c("compName"="character", "error"="character"))
names(score_info) <- gsub("score_summary_", "", gsub(".txt", "", sinfofiles))
score_info <- rbindlist(score_info, idcol="chromosome", fill=TRUE, use.names=TRUE)
score_info[, local_path := NULL]

# Collate numbers that are always present
nasum <- function(...) { base::sum(..., na.rm=TRUE) }
collated_info <- score_info[, .(n_chr=nasum(n_chr), n_match=nasum(n_match), n_multiallele=nasum(n_multiallele),
                                n_multiallele_removed=nasum(n_multiallele_removed),
                                n_ambig=nasum(n_ambig), n_ambig_kept=nasum(n_ambig_kept), n_used=nasum(n_used)),
                            by=.(path, compName, rsid, chr, pos, EA, EAF, OA, weight, is_dom, is_rec, n_var)]

# Collate numbers that are sometimes present:
if ("n_interaction_skipped" %in% names(score_info)) {
  col2 <- score_info[, .(n_interaction_skipped=nasum(n_interaction_skipped)), by=compName]
  collated_info <- collated_info[col2, on = .(compName)]
  collated_info <- collated_info[, .(path, compName, rsid, chr, pos, EA, EAF, OA, weight, is_dom, is_rec,
                                     n_var, n_chr, n_interaction_skipped, n_match, n_ambig, n_ambig_kept,
                                     n_used)]
} 

# Collate errors: first those that are the same across all chromosomes (should always be the case,
# unless theres' been I/O blocking shenanigans (like the disk quota running out partway through a run).
all_errors <- score_info[,.(error=unique(error)), by=compName]
n_errors <- all_errors[,.(N=length(unique(error))), by=compName]
errors <- all_errors[n_errors[N == 1, .(compName)], on=.(compName)]

# Collate errors that are different on each chromosome into a single error message per score
diff_errors <- score_info[n_errors[N > 1, .(compName)], on = .(compName), .(compName, chromosome, error)][!is.na(error)]
diff_errors <- diff_errors[, .(msg = sprintf("For chromosomes %s:", paste(chromosome, collapse=", "))), by=.(compName, error)]
diff_errors <- diff_errors[, .(error = paste(sprintf("%s %s.", msg, error), collapse=" ")), by=.(compName)]
errors[diff_errors, on=.(compName), error := i.error] # Add to global error table

# Collate warnings generated when the score has variants on a chromosome for which there
# was no genotype data (or where that file could not be opened).
if ("chr_fail" %in% names(score_info)) {
  chr_fail <- score_info[(chr_fail), 
    .(error = sprintf("Variants on chromosome(s) %s not counted due to failure to locate or open genotype data for those chromosome(s).", 
                      paste(chromosome, collapse=", "))),
     by=compName]

  errors[chr_fail, on = .(compName), error := paste(error, i.error, sep=". ")]
  errors[, error := gsub("^NA. ", "", error)]
}

# Add errors to collated_info
collated_info[errors, on = .(compName), error := error]

# Give a name to each score that we may use in the output
collated_info[, name := gsub("\\.txt(\\..*)?", "", basename(path))]
collated_info[as.integer(compName) != as.numeric(compName), name := weight]

# We also want the integer version of the compName to handle multi-score files
collated_info[, compNameInt := as.integer(compName)]

# Determine path for score file output(s)
if (is.null(args[["--out"]]) && is.null(args[["--single-out"]])) {
  # The default. The levels of each score are stored in the same location in which the score-file are
  # stored.
  if (is.null(args[["--cohort-name"]])) {
    collated_info[, outpath := sprintf("%s/sample_levels", dirname(path))]
    fname <- collated_info[compName == compNameInt][,.(compNameInt, outpath = sprintf("%s/%s.sscore.gz", outpath, name))]
    collated_info[fname, on = .(compNameInt), outpath := i.outpath]
    rm(fname)
  } else {
    collated_info[, outpath := sprintf("%s/%s_sample_levels", dirname(path), args[["--cohort-name"]])]
    fname <- collated_info[compName == compNameInt][,.(compNameInt, outpath = sprintf("%s/%s_%s.sscore.gz", outpath, name, args[["--cohort-name"]]))]
    collated_info[fname, on = .(compNameInt), outpath := i.outpath]
    rm(fname)
  }

  # if there are multiple score files in that directory create new sub folders named for the score.
  n_files <- collated_info[compName == compNameInt][, .(N=length(list.files(path=dirname(path), pattern="*\\..*"))), by=.(compNameInt)]
  newpath <- collated_info[compName == compNameInt][n_files[N > 1], on=.(compNameInt), 
     .(compNameInt, outpath = sprintf("%s/%s/%s", dirname(outpath), name, basename(outpath)))]
  collated_info[newpath, on=.(compNameInt), outpath := i.outpath]
} else if (!is.null(args[["--out"]]) && is.null(args[["--single-out"]])) {
  # Single root output directory given, but still want to store each input score-file separately 
  if (is.null(args[["--cohort-name"]])) {
    collated_info[, outpath := sprintf("%s/sample_levels", args[["--out"]])]
    fname <- collated_info[compName == compNameInt][,.(compNameInt, outpath = sprintf("%s/%s/%s.sscore.gz", outpath, name, name))]
    collated_info[fname, on = .(compNameInt), outpath := i.outpath]
    rm(fname)
  } else {
    collated_info[, outpath := sprintf("%s/%s_sample_levels", args[["--out"]], args[["--cohort-name"]])]
    fname <- collated_info[compName == compNameInt][,.(compNameInt, outpath = sprintf("%s/%s/%s_%s.sscore.gz", outpath, name, name, args[["--cohort-name"]]))]
    collated_info[fname, on = .(compNameInt), outpath := i.outpath]
    rm(fname)
  }
} else if (is.null(args[["--out"]]) && !is.null(args[["--single-out"]])) {
  # Single output file requested, but output directory not given (default to directory --score-file is in)
  if (args[["--type"]] == "d") {
    root_dir <- args[["--score-file"]]
  } else {
    root_dir <- dirname(args[["--score-file"]])
  }
  if (is.null(args[["--cohort-name"]])) {
    collated_info[, outpath := sprintf("%s/sample_levels/%s.sscore.gz", root_dir, args[["--single-out"]])]
  } else {
    collated_info[, outpath := sprintf("%s/%s_sample_levels/%s.sscore.gz", root_dir, args[["--cohort-name"]], args[["--single-out"]])]
  }
} else {
  # both --out and --single-out given
  collated_info[, outpath := sprintf("%s/%s.sscore.gz", args[["--out"]], args[["--single-out"]])]
}

# if there are files with only one score, we will name these as "score_sum" in the score-file output.
n_scores <- collated_info[,.N,by=outpath]
collated_info[n_scores[N == 1], on = .(outpath), name := "score_sum"]

# conversely, if there are any duplicated names (for example when requesting a single output file),
# make them unique:
dup_names <- collated_info[weight != "m",.N,by=.(name, outpath)][N > 1]
collated_info[dup_names, on = .(name, outpath), name := paste0(name, ".", seq_len(.N)), by=.(name, outpath)]

# rename some of the columns and drop any unused ones:
setnames(collated_info, c("rsid", "chr", "pos", "EA", "EAF", "OA", "weight", "is_dom", "is_rec"),
  c("rsid_column", "chr_column", "pos_column", "effect_allele_column", "EAF_column", "other_allele_column",
    "weight_column", "is_dominant_column", "is_recessive_column"))
collated_info[, compNameInt := NULL]
setnames(collated_info, "name", "score_name")
setnames(collated_info, "path", "score_path")
collated_info <- collated_info[, which(!sapply(collated_info, function(x) { all(is.na(x)) })), with=FALSE] # drop all columns that are filled with all NA
if (!("error" %in% names(collated_info))) collated_info[, error := NA_character_] # need to keep this column for later errors
collated_info[weight_column == "m", c("weight_column", "score_name") := .(NA, "(Totals for input multi-score):")]
collated_info <- collated_info[, c("score_path", "score_name", setdiff(names(collated_info), c("score_path", "score_name"))), with=FALSE]

# Copy across information about scores run:
if (args[["--type"]] == 'l') {
  dx_download(args[["--score-file"]], "output/run_score_file_list.txt")
} else if (args[["--type"]] == 'd') {
  cat(sep="", gsub("/$", "", args[["--score-file"]]), "/*\n", file="output/run_score_file_list.txt")
} else {
  cat(sep="", args[["--score-file"]], "\n", file="output/run_score_file_list.txt")
}

# Copy this script across to freeze version
system("cp calc_PS_lvls.R output/", wait=TRUE)

# Get launcher scripts and command
dx_download(sprintf("%s/command_log.txt", args[["--work"]]), "output/")
dx_download(sprintf("%s/calc_PS_lvls.sh", args[["--work"]]), "output/")

# if there is only a single output requested, and its the same
# as the working directory then all we need to do is dcast the
# score file and clean up.
outfiles <- unique(collated_info$outpath)
if (length(outfiles) == 1) {
  outfile <- unique(outfiles)
  out_name <- basename(outfile)
  out_dir <- paste0(dirname(outfile), "/")
  # give the scores their names
  sscores <- sscores[collated_info[, .(compName, score_name)], on = .(compName), nomatch=0, .(IID, score_name, score_sum)]
  sscores <- dcast(sscores, IID ~ score_name, value.var="score_sum")
  sscores <- sscores[, intersect(c("IID", collated_info$score_name), names(sscores)), with=FALSE] # preserve score order
  tryCatch({
    fwrite(sscores, sep="\t", quote=FALSE, compress="gzip", file=sprintf("output/%s", out_name))
  }, error = function(e) {
    collated_info[, error := paste0(error, ". Error when trying to write wide-format collated scores to", outfile)]
  })

  # Clean up collated info
  dropnames <- c("compName", "outpath", "score_fail")
  if (collated_info[,all(is.na(error))]) dropnames <- c(dropnames, "error")
  if (collated_info[,all(score_name == "score_sum")]) dropnames <- c(dropnames, "score_name")
  fwrite(collated_info[, setdiff(names(collated_info), dropnames), with=FALSE],
         sep="\t", quote=FALSE, file="output/score_summary.txt")

	# Upload to project storage
	dx_upload("output/", out_dir)
	
	# Remove working files
	if (!dx_exists(sprintf("%s/errors/", args[["--work"]])) && args[["--work"]] != out_dir) {
	  dx_rm(args[["--work"]])
	} else {
	  dx_rm(sprintf("%s/checkpoint1", args[["--work"]]))
	  dx_rm(sprintf("%s/checkpoint2", args[["--work"]]))
	  dx_rm(sprintf("%s/checkpoint3", args[["--work"]]))
	  dx_rm(sprintf("%s/finished", args[["--work"]]))
	  dx_rm(sprintf("%s/args.rds", args[["--work"]]))
	}
	
	# Can now exit
  quit(save="no")
}

# For each output directory, attempt to:
#
# (1) save the sscore file
# (2) save the relevant subset of the information file
# (3) copy the collated plink logs
#
# We only do this if the file does not exist, otherwise we preserve
# the working directory and error out.
for (outIdx in unique(collated_info$outpath)) {
  out_dir <- paste0(dirname(outIdx), "/")
  if (dx_exists(out_dir)) {
    collated_info[outpath == outIdx, error := paste0(error, ". Output directory ", out_dir, " already exists!")] 
    next
  }
 
  # Extract subset of scores to write out and cast to wide format
  this_info <- collated_info[outpath == outIdx]
  this_info[, score_fail := FALSE]
  this_info[!(compName %in% unique(sscores$compName)) || is.na(weight_column), score_fail := TRUE]
  this_sscore <- sscores[this_info[, .(compName, score_name)], on = .(compName), nomatch=0, .(IID, score_name, score_sum)]
  this_sscore <- dcast(this_sscore, IID ~ score_name, value.var="score_sum")
  this_sscore <- this_sscore[, intersect(c("IID", collated_info$score_name), names(this_sscore)), with=FALSE] # preserve score order
  
  # Write out collated score to output folder for upload
  fwrite(this_sscore, sep="\t", quote=FALSE, compress="gzip", file=sprintf("output/%s", basename(outIdx)))

  # Write out subset of the score info
  this_info <- this_info[!(score_fail)]
  dropnames <- c("compName", "outpath", "score_fail")
  if (this_info[,all(is.na(error))]) dropnames <- c(dropnames, "error")
  if (this_info[,all(score_name == "score_sum")]) dropnames <- c(dropnames, "score_name")
  fwrite(this_info[, setdiff(names(this_info), dropnames), with=FALSE],
         sep="\t", quote=FALSE, file="output/score_summary.txt")
  
  # Upload to project storage
  dx_upload("output/", out_dir)
  
  # Remove score_summary and collated score file for next iteration
  system("rm output/score_summary.txt")
  system(sprintf("rm output/%s", basename(outIdx)))

  # Remove successfully transferred scores from the collated information
  collated_info <- collated_info[!this_info, on=.(compName)]
}
collated_info[, error := gsub("^NA\\. ", "", error)]
collated_info[, error := gsub("^\\. ", "", error)]
collated_info[, error := gsub("^ \\. ", "", error)]

# Save error output for individual scores in errors directory
if (nrow(collated_info) > 0) {
  fwrite(collated_info, sep="\t", quote=FALSE, file="errors/score_summary.txt")
  fwrite(sscores, sep="\t", quote=FALSE, compress="gzip", file="errors/collated_scores.sscore.gz")
  dx_upload("errors", args[["--work"]])
  
}

# Remove working files
dx_rm(sprintf("%s/checkpoint1", args[["--work"]]))
dx_rm(sprintf("%s/checkpoint2", args[["--work"]]))
dx_rm(sprintf("%s/checkpoint3", args[["--work"]]))
dx_rm(sprintf("%s/finished", args[["--work"]]))
dx_rm(sprintf("%s/args.rds", args[["--work"]]))

# Now finally throw error if needed
if (nrow(collated_info) > 0) {
  stop("Some scores had errors: see ", args[["--work"]], "errors/score_summary.txt")
}
