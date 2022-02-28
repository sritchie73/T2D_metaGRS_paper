library(data.table)
library(foreach)
source("src/functions/define_case_control_status.R")

# curate set of case control definitions and models used for metaGRS training 
#
# Here we drop QDiabetes models - these use the same definition as
# short_name == "incident_basic_adjudicated" but additionally drop samples with
# missing data in the given QDiabetes column - i.e. used later only for testing
# metaGRS
case_control_definitions <- case_control_definitions[!(short_name %like% "QDiabetes")]

# Collate all candidate metaGRS into three files: 
#
#  (1) for metaGRS trained on prevalent T2D
#  (2) for metaGRS trained on incident T2D
#  (3) for metaGRS trained on prevalent + incident T2D combined
#
# We do this to minimise the number of PGS level calculation programs we need to run.
# The PGS calculation script can calculate many PGS simultaneously, which is much faster
# than running them separately. However, we are also constrained by memory, so better to
# split into a few files than try to run all ~1200 scores at once.
#
for (this_model_type in unique(case_control_definitions$type)) {
  this_model <- case_control_definitions[type == this_model_type]
  scores <- foreach(mname = this_model$short_name, .combine=rbind) %do% {
    foreach(prefilter = c("none", "sig_auc", "no_auc_ci_overlap"), .combine=rbind) %do% {
      this_scores <- fread(sprintf("output/metaGRS/train/%s/prefilter_%s/candidate_metaGRS.txt.gz", mname, prefilter), tmpdir="tmp")
      setnames(this_scores, gsub("alpha", sprintf("%s_%s_alpha", mname, prefilter), gsub("\\.", "_", names(this_scores))))
      melt(this_scores, id.vars=c("rsid", "chr", "pos", "effect_allele", "other_allele"), variable.name="metaGRS", value.name="beta", na.rm=TRUE)
    } 
  }   
  scores <- dcast(scores, rsid + chr + pos + effect_allele + other_allele ~ metaGRS, value.var="beta", fill=0)
  scores <- scores[order(pos)][order(chr)]
  fwrite(scores, sep="\t", quote=FALSE, compress="gzip", file=sprintf("output/metaGRS/train/collated_candidate_metaGRSs_%s_T2D.txt.gz", gsub(" \\+ ", "_plus_", this_model_type)))
}
