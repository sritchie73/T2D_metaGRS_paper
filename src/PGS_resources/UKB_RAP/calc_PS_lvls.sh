#!/usr/bin/env bash
#
# Requires docopts in your $PATH; download at https://github.com/docopt/docopts
# Also requires R with dxutils (https://github.com/sritchie73/dxutils) installed

# Location on RAP project storage where this and the calc_PS_lvls.R script are
# located
src_dir=src/PGS_resources/UKB_RAP

eval "$(docopts -h - : "$@" <<EOF
Calculate the levels of a polygenic score in UK Biobank on the RAP

By default, this program calculates polygenic score levels using 
the TopMed imputation (extracted into plink2 binary format) by
matching variants by chromosome, position, and alleles between
the genotype data (build GRCh38) and the score file (must
be the same genome build). Variants with ambiguous alleles (A/T or
G/C SNPs are excluded with warning). Please read through the list
of options to modify these default behaviours.

This script uses the run_script applet to run the calc_PGS_lvls_RAP.R script
on a cloud workstation batched by chromosome. Sensible defaults are provided.

NOTE: when calculating multiple scores, please submit one job using
the multiple score functionality.

Usage:
  calc_PS_lvls.sh --score-file <file> [options]
  calc_PS_lvls.sh -h | --help

Options:
  -h --help                   Show this screen.
  --score-file <file>         Path to polygenic score file, directory, or file containing list
                              of score files (see --type) on the RAP project storage.
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
  --work <directory>          Working directory on RAP project storage to store intermediate files 
                              and logs shared across all scores for the duration of the run. The 
                              default is to store this in a directory called 'work/' in the folder 
                              in which the input --score-file is scored, then all logs will be copied 
                              to each score output folder under 'checkpointing/' at the end of the 
                              run. [default: NULL]
  --genotype-prefix <prefix>  Path and prefix occurring before the chromosome number for the genotype
                              data on RAP project storage you want to use for polygenic scoring.
                              Defaults to the plink2 binary files extracted for the TopMed imputation
                              in the CEU_overarching project.
                              [default: 'common/Imputed Genotypes/TopMed (GRCh38)/ukb21007_c']
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
  --instance-type <name>      Instance type to use for each of the 26 parallel jobs (1 per chromosome) 
                              submitted via dx run run_script [default: mem3_ssd1_v2_x8]
  --priority <type>           Priority for the job submitted by dx run run_script [default: low]
EOF
)"

# Check type switch
if ! [[ $type = 's' || $type = 'd' || $type = 'l' ]]; then
  echo "--type must be one of 's', 'd', or 'l'. See --help." >&2
  exit 1
fi


# Check for logging/working directory
if [[ $work = "NULL" ]]; then
  indir=$(dirname $score_file)
  if [[ $indir = "." ]]; then
    work=checkpointing/
  else 
    work=$indir/checkpointing/
  fi
fi
if [[ $work != */ ]]; then
  work=$work/
fi
if [[ $(Rscript -e "dxutils::dx_exists('"$work"')") = "[1] TRUE" ]]; then
  echo "Working directory $work already exists on RAP project storage. Overwrite? (y/n)" 1>&2
  read ans
  while true; do
    if [[ $ans = "y" || $ans = "Y" || $ans = "Yes" || $ans = "YES" || $ans = "yes" ]]; then
      Rscript -e "dxutils::dx_rm('"$work"')"
      if [[ $? -ne 0 ]]; then
        exit 1
      fi
      break
    elif [[ $ans = "n" || $ans = "N" || $ans = "No" || $ans = "NO" || $ans = "no" ]]; then
      echo "Resume from existing checkpointing data in $work? (y/n)"
      read ans
      while true; do
        if [[ $ans = "y" || $ans = "Y" || $ans = "Yes" || $ans = "YES" || $ans = "yes" ]]; then
          if [[ $(Rscript -e "dxutils::dx_exists('"$work"/checkpoint1')") = "[1] FALSE" ]]; then
            echo "No checkpointing data found, aborting"
            exit 1
          fi
          tasks=()
          echo "Determining which chromosomes still need to be computed..."
          for ii in {1..22} "X"; do
            if [[ $(Rscript -e "dxutils::dx_exists('"$work"/finished/score_summary_"$ii".txt')") = "[1] FALSE" ]]; then
              if [[ $ii = "X" ]]; then
                tasks+=(23)
                echo "Chromosome 23 added to task list"
              else
                tasks+=($ii)
                echo "Chromosome $ii added to task list"
              fi
            fi
          done
          if [[ ${#tasks[@]} -eq 0 ]]; then
            echo "All chromosomes completed PRS computation, launching single task with chromosome 22 label to collate results"
            tasks=(22)
          fi
          break
        elif [[ $ans = "n" || $ans = "N" || $ans = "No" || $ans = "NO" || $ans = "no" ]]; then
          exit 1
        else
         echo "Unrecognised user input. Please answer 'y' or 'n'." 1>&2
         read ans
        fi
      done
    else
     echo "Unrecognised user input. Please answer 'y' or 'n'." 1>&2
     read ans
    fi
  done
else 
  echo "Working and temporary logging directory on RAP project storage is: $work" 1>&2
  
  # Log submitted command
  arg_string=$@
  echo "Batch command:" > command_log.txt
  echo "------------------------------------------------------------------------" >> command_log.txt
  echo "$src_dir/calc_PS_lvls.sh $arg_string" >> command_log.txt
  echo "" >> command_log.txt
  Rscript -e "dxutils::dx_upload('command_log.txt', '"$work"')"
  rm command_log.txt
  
  # Copy across this script file
  Rscript -e "dxutils::dx_upload('"$src_dir"/calc_PS_lvls.sh', '"$work"')"
  if [[ $? -ne 0 ]]; then
    exit 1
  fi
  
  # Create task list
  tasks=$(seq 1 23)
fi

# build command string
cmd[0]="Rscript calc_PS_lvls.R"
cmd[1]="--score-file $score_file"
cmd[3]="--work $work"
cmd[4]="--type $type"
cmd[5]="--score-rsid $score_rsid"
cmd[6]="--score-chr $score_chr"
cmd[7]="--score-pos $score_pos"
cmd[8]="--score-EA $score_EA"
if [[ $score_EAF != "NULL" ]]; then  cmd[9]="--score-EAF $score_EAF"; fi
cmd[10]="--score-OA $score_OA"
cmd[11]="--score-weight $score_weight"
if [[ $score_dominant != "NULL" ]]; then  cmd[12]="--score-dominant $score_dominant"; fi
if [[ $score_recessive != "NULL" ]]; then  cmd[13]="--score-recessive $score_recessive"; fi
if $match_by_rsid; then  cmd[14]="--match-by-rsid"; fi
cmd[15]="--cohort-name $cohort_name"
if [[ $out != "NULL" ]]; then  cmd[16]="--out $out"; fi
if [[ $single_out != "NULL" ]]; then  cmd[17]="--single-out $single_out"; fi
cmd[18]="--genotype-prefix $genotype_prefix"
if [[ $genotype_suffix != "NULL" ]]; then  cmd[19]="--genotype-suffix $genotype_suffix"; fi
cmd[20]="--genotype-format $genotype_format"
if $single_geno; then  cmd[21]="--single-geno"; fi
if [[ $keep != "NULL" ]]; then  cmd[22]="--keep $keep"; fi
if $keep_ambiguous; then  cmd[23]="--keep-ambiguous"; fi
if [[ $ambiguous_thresh != "NULL" ]]; then  cmd[24]="--ambiguous-thresh $ambiguous_thresh"; fi
if [[ $freqx_prefix != "NULL" ]]; then  cmd[25]="--freqx-prefix $freqx_prefix"; fi
if [[ $freqx_suffix != "NULL" ]]; then  cmd[26]="--freqx-suffix $freqx_suffix"; fi
if $remove_multiallelic; then  cmd[27]="--remove-multiallelic"; fi

cmd_string=${cmd[@]}

# Create array job
for task_id in ${tasks[@]}; do
  job_id=$(dx run run_script \
    --name "Calculate PGS, chromosome $task_id" \
    -iscript="$src_dir/calc_PS_lvls.R" \
    -icmd="$cmd_string" \
    -ienv="SLURM_ARRAY_TASK_MAX=23" \
    -ienv="SLURM_ARRAY_TASK_ID=$task_id" \
    --instance-type="$instance_type" \
    --priority="$priority" \
    --brief --yes)
  echo "Job to calculate PGS on chromosome $task_id submitted with DNAnexus job ID: $job_id"
done
