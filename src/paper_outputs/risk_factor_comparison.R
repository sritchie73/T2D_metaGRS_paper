library(data.table)
library(ggplot2)
library(cowplot)

system("mkdir -p output/paper_outputs/{figures,tables}", wait=TRUE)

# Load assocations
uni <- fread("output/UKB_tests/incident_T2D_associations_QDiabetes2018C_subset.txt")
multi <- fread("output/UKB_tests/incident_T2D_multivariate_associations.txt")

# Curate C-indices
cinds <- unique(rbind(uni[,.(model_type, model, C.index, C.L95, C.U95)], multi[,.(model_type, model, C.index, C.L95, C.U95)]))
cinds <- cinds[model %in% c("age + sex + assessment_centre", "T2D_metaGRS", "QDiabetes2018A", "QDiabetes2018B_non_fasting", "QDiabetes2018C", "QDiabetes2013") |
  (model_type == "risk factor" & model != "fasting_glucose")]

# Create labels
cinds[uni[,.SD[1],by=.(model_type, model)], on = .(model_type, model), coefficient := i.coefficient]
cinds[coefficient == "Age", coefficient := "Reference (age + sex)"]
cinds[coefficient == "Body mass index", coefficient := "Body mass index (BMI)"]
cinds[coefficient == "Smoking status: Moderate vs. Non-smoker", coefficient := "Smoking status"]
cinds[coefficient == "QDiabetes 2018 model B (non-fasting glucose)", coefficient := "QDiabetes 2018 model B"]
cinds[model_type == "QDiabetes2013 + T2D_metaGRS", coefficient := "QDiabetes 2013 + T2D metaPRS"]
cinds[model_type == "QDiabetes2018A + T2D_metaGRS", coefficient := "QDiabetes 2018 model A + T2D metaPRS"]
cinds[model_type == "QDiabetes2018B + T2D_metaGRS", coefficient := "QDiabetes 2018 model B + T2D metaPRS"]
cinds[model_type == "QDiabetes2018C + T2D_metaGRS", coefficient := "QDiabetes 2018 model C + T2D metaPRS"]
cinds[, model_type := fcase(
  model_type == "reference", "Reference (age + sex)",
  model_type == "PGS", "Genetics",
  model_type == "risk factor", "Established risk factors",
  default = "Combined scores")]

# Order
cinds[, model_type := factor(model_type, levels=c("Reference (age + sex)", "Genetics", "Established risk factors", "Combined scores"))]
cinds <- cinds[order(C.index)][order(model_type)]
cinds[, coefficient := factor(coefficient, levels=rev(coefficient))]

# Plot 
cinds[, color_anno := model == "T2D_metaGRS"]
ref_lines <- cinds[coefficient %in% c("Reference (age + sex)", "T2D_metaGRS", "QDiabetes 2018 model C"), .(C.index)]

g <- ggplot(cinds) +
  aes(x=C.index, xmin=C.L95, xmax=C.U95, y=coefficient, color=color_anno) +
  facet_grid(model_type ~ ., space="free_y", scales="free_y") + 
  geom_vline(data=ref_lines, aes(xintercept=C.index), linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, fill="white") +
  scale_color_manual(values=c("TRUE"="#bd0026", "FALSE"="#08306b")) +
  xlab("C-index (95% CI)") +
  theme_bw() +
  theme(
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
    axis.text.y=element_text(size=6, color="black"), axis.title.y=element_blank(),
    strip.background=element_blank(), strip.text=element_blank(),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="none"
  )

ggsave(g, width=3.6, height=4, file="output/paper_outputs/figures/risk_factor_comparison.pdf")

# Curate supplementary table
dt <- rbind(
  uni[,.(model_type, model, Samples, Cases, coefficient, HR, L95, U95, Pvalue, C.index, C.L95, C.U95)],
  multi[,.(model_type, model, Samples, Cases, coefficient, HR, L95, U95, Pvalue, C.index, C.L95, C.U95)]
)
dt <- dt[model %in% c("age + sex + assessment_centre", "T2D_metaGRS", "QDiabetes2018A", "QDiabetes2018B_non_fasting", "QDiabetes2018C", "QDiabetes2013") |
  (model_type == "risk factor" & model != "fasting_glucose")]

dt <- dt[!(coefficient %like% "Assessment centre")]
dt <- dt[!(coefficient %in% c("Age", "Sex: Male vs. Female")) | model_type == "reference"]

dt[uni[,.SD[1],by=.(model_type, model)], on = .(model_type, model), model := i.coefficient]
dt[model == "Age", model := "Reference (age + sex)"]
dt[model == "Body mass index", model := "Body mass index (BMI)"]
dt[model == "Smoking status: Moderate vs. Non-smoker", model := "Smoking status"]
dt[model == "QDiabetes 2018 model B (non-fasting glucose)", model := "QDiabetes 2018 model B"]
dt[model == "QDiabetes 2018 model A", coefficient := "QDiabetes 2018 model A"]
dt[model == "QDiabetes 2018 model B", coefficient := "QDiabetes 2018 model B"]
dt[model == "QDiabetes 2018 model C", coefficient := "QDiabetes 2018 model C"]
dt[coefficient == "QDiabetes   2013", coefficient := "QDiabetes 2013"]
dt[coefficient == "QDiabetes 2018 model B_non_fasting", coefficient := "QDiabetes 2018 model B"]
dt[model_type == "QDiabetes2013 + T2D_metaGRS", model := "QDiabetes 2013 + T2D metaPRS"]
dt[model_type == "QDiabetes2018A + T2D_metaGRS", model := "QDiabetes 2018 model A + T2D metaPRS"]
dt[model_type == "QDiabetes2018B + T2D_metaGRS", model := "QDiabetes 2018 model B + T2D metaPRS"]
dt[model_type == "QDiabetes2018C + T2D_metaGRS", model := "QDiabetes 2018 model C + T2D metaPRS"]
dt[, model_type := fcase(
  model_type == "reference", "Reference (age + sex)",
  model_type == "PGS", "Genetics",
  model_type == "risk factor", "Established risk factors",
  default = "Combined scores")]
dt <- unique(dt)

# Order
dt[, model_type := factor(model_type, levels=c("Reference (age + sex)", "Genetics", "Established risk factors", "Combined scores"))]
dt <- dt[order(-HR)][order(-C.index)][order(model_type)]

# Write out
fwrite(dt, sep="\t", quote=FALSE, file="output/paper_outputs/tables/risk_factor_comparison.txt")

# Write out delta-C index
delta_C <- fread("output/UKB_tests/delta_cind_bootstraps.txt")
fwrite(delta_C, sep="\t", quote=FALSE, file="output/paper_outputs/tables/qdiabetes_deltaC_boot.txt")

