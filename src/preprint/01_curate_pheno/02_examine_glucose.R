library(data.table)
library(ggplot2)
library(cowplot)
source("src/functions/glm_test.R")
source("src/functions/cox_test.R")
source("src/functions/factor.R") # set reference group (first level) without specifying all levels
source("src/functions/factor_by_size.R") # set reference group (first level) as largest group

# Make output directories
system("mkdir -p output/slide_deck_plots", wait=TRUE)
system("mkdir -p output/UKB_tests", wait=TRUE)

# Load curated phenotype data
pheno <- fread("data/UKB/collated_curated_data.txt")

# Extract glucose
dt <- melt(pheno, id.vars=c("eid", "fasting_time", "sex"), 
  measure.vars=c("fasting_glucose", "non_fasting_glucose"),
  na.rm=TRUE, variable.name="glucose_type", 
  value.name="glucose_concentration")

# Process
dt[, glucose_type := gsub("_glucose", "", glucose_type)]
dt[, glucose_type := gsub("_", "-", glucose_type)]
dt[, glucose_type := factor(glucose_type, levels=c("non-fasting", "fasting"))]

# Bin fasting time
dt[, fast_gt_8 := fasting_time]
dt[fast_gt_8 > 8, fast_gt_8 := 8]
dt[, fast_gt_8 := paste(fast_gt_8, "hours")]
dt[fast_gt_8 == "1 hours", fast_gt_8 := "1 hour"]
dt[fast_gt_8 == "8 hours", fast_gt_8 := "≥ 8 hours"]
dt[, fast_gt_8 := factor(fast_gt_8, levels=c(
  "0 hours", "1 hour", "2 hours", "3 hours", "4 hours",
  "5 hours", "6 hours", "7 hours", "≥ 8 hours"
))]

dt[, fast_gt_3 := fasting_time]
dt[fast_gt_3 > 3, fast_gt_3 := 3]
dt[, fast_gt_3 := paste(fast_gt_3, "hours")]
dt[fast_gt_3 == "1 hours", fast_gt_3 := "1 hour"]
dt[fast_gt_3 == "3 hours", fast_gt_3 := "≥ 3 hours"]
dt[, fast_gt_3 := factor(fast_gt_3, levels=c(
  "0 hours", "1 hour", "2 hours", "≥ 3 hours"
))]

# Compare glucose concentrations before and after adjustment
g <- ggplot(dt, aes(x=glucose_concentration, color=fast_gt_3)) +
  geom_density(trim=TRUE) +
  facet_grid(glucose_type ~ sex) +
  scale_color_manual("Fasting time", values=structure(
    colorRampPalette(c("#9ecae1", "#08306b"))(4),
    name=levels(dt$fast_gt_3)
  )) +
  xlab("Glucose (mmol/L)") + 
  xlim(0, 10) +
  theme_bw()
ggsave(g, width=13, height=7.2, file="output/slide_deck_plots/UKB_glucose_adj_fast_lte_3.png")

g <- ggplot(dt, aes(x=glucose_concentration, color=fast_gt_8)) +
  geom_density(trim=TRUE) +
  facet_grid(glucose_type ~ sex) +
  scale_color_manual("Fasting time", values=structure(
    colorRampPalette(c("#9ecae1", "#08306b"))(9),
    name=levels(dt$fast_gt_8)
  )) +
  xlab("Glucose (mmol/L)") + 
  xlim(0, 10) +
  theme_bw()
ggsave(g, width=13, height=7.2, file="output/slide_deck_plots/UKB_glucose_adj_fast_lte_8.png")

# ---------------------------------
# Model fits for glucose vs. T2D
# ---------------------------------
pheno[, age := scale(age)]
pheno[, sex := factor(sex, reference="Female")]
pheno[, assessment_centre := factor_by_size(assessment_centre)]
pheno[, censor_hospital_nation := factor_by_size(censor_hospital_nation)]
pheno[is.na(earliest_hospital_nation) & !(type_2_diabetes), earliest_hospital_nation := assessment_nation]
pheno[, earliest_hospital_nation := factor_by_size(earliest_hospital_nation)]

# Drop participants free of T2D at baseline who have withdrawn consent for hospital record
# linkage and people with uncertain diabetes status
pheno <- pheno[!is.na(type_2_diabetes) & !is.na(earliest_hospital_nation)]

# Fit glucose on prevalent disease
mf <- "type_2_diabetes ~ %s + age + sex + assessment_centre + earliest_hospital_nation"
prev <- rbind(idcol="model",
  "Reference"=glm.test(type_2_diabetes ~ age + sex + assessment_centre + earliest_hospital_nation, "type_2_diabetes", pheno, ci.method="wald"),
  "Non-fasting glucose"=glm.test(sprintf(mf, "non_fasting_glucose"), "type_2_diabetes", pheno, ci.method="wald"),
  "Glucose (fasting time ≥ 3 hours)"=glm.test(sprintf(mf, "non_fasting_glucose"), "type_2_diabetes", pheno[fasting_time >= 3], ci.method="wald"),
  "Glucose (fasting time ≥ 8 hours)"=glm.test(sprintf(mf, "non_fasting_glucose"), "type_2_diabetes", pheno[fasting_time >= 8], ci.method="wald"),
  "Adjusted glucose"=glm.test(sprintf(mf, "fasting_glucose"), "type_2_diabetes", pheno, ci.method="wald")
)

# Fit glucose on incident disease
y <- "Surv(incident_censor_years, incident_type_2_diabetes)"
follow_diff <- "assessment_centre + earliest_hospital_nation + censor_hospital_nation"
mf <- paste(y, "~ %s + age + sex +", follow_diff)
inci <- rbind(idcol="model",
  "Reference"=cox.test(paste(y, "~ age + sex +", follow_diff), "incident_type_2_diabetes", pheno),
  "Non-fasting glucose"=cox.test(sprintf(mf, "non_fasting_glucose"), "incident_type_2_diabetes", pheno),
  "Glucose (fasting time ≥ 3 hours)"=cox.test(sprintf(mf, "non_fasting_glucose"), "incident_type_2_diabetes", pheno[fasting_time >= 3]),
  "Glucose (fasting time ≥ 8 hours)"=cox.test(sprintf(mf, "non_fasting_glucose"), "incident_type_2_diabetes", pheno[fasting_time >= 8]),
  "Adjusted glucose"=cox.test(sprintf(mf, "fasting_glucose"), "incident_type_2_diabetes", pheno)
)

# Write out model fits
fwrite(prev, sep="\t", quote=FALSE, file="output/UKB_tests/glucose_vs_prevalent_T2D.txt")
fwrite(inci, sep="\t", quote=FALSE, file="output/UKB_tests/glucose_vs_incident_T2D.txt")

# Plot
prev_null <- prev[model == "Reference"][1]
inci_null <- inci[model == "Reference"][1]
prev <- prev[coefficient %like% "glucose"]
inci <- inci[coefficient %like% "glucose"]

g_cind <- ggplot(inci) +
  aes(x=C.index, xmin=C.L95, xmax=C.U95, y=model) +
  geom_rect(inherit.aes=FALSE, data=inci_null,
          aes(xmin=C.L95, xmax=C.U95),
          ymin=-Inf, ymax=Inf, fill="#bdbdbd", alpha=0.3) +
  geom_vline(data=inci_null, aes(xintercept=C.index), linetype=2) +
  geom_errorbarh(height=0, alpha=0.8, color="#08519c") +
  geom_point(shape=19, color="#08519c") +
  xlab("C-index (95% CI)") +
  ylab("") +
  theme_bw() +
  theme(
    legend.position="none",
    panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    axis.text=element_text(size=8),
    axis.title=element_text(size=10)
  )

g_hr <- ggplot(inci) +
  aes(x=HR, xmin=L95, xmax=U95, y=model) +
  geom_vline(xintercept=1, linetype=2) +
  geom_errorbarh(height=0, alpha=0.8, color="#08519c") +
  geom_point(shape=19, color="#08519c") +
  xlab("Hazard Ratio (95% CI)") +
  ylab("") +
  theme_bw() +
  theme(
    legend.position="none",
    panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    axis.text=element_text(size=8),
    axis.title=element_text(size=10)
  )

g_auc <- ggplot(prev) +
  aes(x=AUC, xmin=AUC.L95, xmax=AUC.U95, y=model) +
  geom_rect(inherit.aes=FALSE, data=prev_null,
          aes(xmin=AUC.L95, xmax=AUC.U95),
          ymin=-Inf, ymax=Inf, fill="#bdbdbd", alpha=0.3) +
  geom_vline(data=prev_null, aes(xintercept=AUC), linetype=2) +
  geom_errorbarh(height=0, alpha=0.8, color="#08519c") +
  geom_point(shape=19, color="#08519c") +
  xlab("AUC (95% CI)") +
  ylab("") +
  theme_bw() +
  theme(
    legend.position="none",
    panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    axis.text=element_text(size=8),
    axis.title=element_text(size=10)
  )

g_or <- ggplot(prev) +
  aes(x=OR, xmin=OR.L95, xmax=OR.U95, y=model) +
  geom_vline(xintercept=1, linetype=2) +
  geom_errorbarh(height=0, alpha=0.8, color="#08519c") +
  geom_point(shape=19, color="#08519c") +
  xlab("Odds Ratio (95% CI)") +
  ylab("") +
  theme_bw() +
  theme(
    legend.position="none",
    panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    axis.text=element_text(size=8),
    axis.title=element_text(size=10)
  )

g <- plot_grid(g_cind, g_hr, g_auc, g_or, nrow=2, align="hv")
ggsave(g, width=13, height=6, units="in", file="output/slide_deck_plots/UKB_glucose_adjustment_comparison.png")

