library(data.table)
library(ggplot2)
library(cowplot)
library(ggforce)

# Make output directory
system("mkdir -p output/slide_deck_plots/2023-02-10-highlights/", wait=TRUE)

# Load association results
prev <- fread("output/UKB_tests/prevalent_T2D_associations.txt")
inci <- fread("output/UKB_tests/incident_T2D_associations.txt")

# Extract null model fit
prev_null <- prev[model == "age + sex + follow_diff"][1]
inci_null <- inci[model == "age + sex + follow_diff"][1]

# Extract PGS associations and modle fit
prev_pgs <- prev[coefficient %like% "PGS" | coefficient %like% "GRS"]
inci_pgs <- inci[coefficient %like% "PGS" | coefficient %like% "GRS"]

# Order
prev_pgs <- prev_pgs[order(-AUC)]
prev_pgs[, coefficient := factor(coefficient, levels=unique(coefficient))]

inci_pgs <- inci_pgs[order(-C.index)]
inci_pgs[, coefficient := factor(coefficient, levels=unique(coefficient))]

# Drop older versions of metaGRS
prev_pgs <- prev_pgs[!(coefficient %like% "Gad")]
inci_pgs <- inci_pgs[!(coefficient %like% "Gad")]

# Flag metaGRS (for plot colours)
prev_pgs[, metaGRS := fcase(
  coefficient == "T2D_metaGRS", "current",
  coefficient %like% "Gad", "previous",
  default = "no")]

inci_pgs[, metaGRS := fcase(
  coefficient == "T2D_metaGRS", "current",
  coefficient %like% "Gad", "previous",
  default = "no")]

# Build plots
g_auc <- ggplot(prev_pgs) +
  aes(y=AUC, ymin=AUC.L95, ymax=AUC.U95, x=coefficient, color=metaGRS) +
	geom_rect(inherit.aes=FALSE, data=prev_null,
					aes(ymin=AUC.L95, ymax=AUC.U95),
					xmin=-Inf, xmax=Inf, fill="#bdbdbd", alpha=0.3) +
  geom_hline(data=prev_null, aes(yintercept=AUC), linetype=2) +
  geom_errorbar(width=0, alpha=0.8) +
  geom_point(shape=19) +
  scale_color_manual(values=c("current"="#ae017e", "previous"="#000000", "no"="#045a8d")) +
  xlab("") +
  ylab("AUC (95% CI)") +
  theme_bw() +
  theme(
    legend.position="none",
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=0.5),
    axis.text.y=element_text(size=8),
    axis.title=element_text(size=10)
  )

g_or <- ggplot(prev_pgs) +
  aes(y=OR, ymin=OR.L95, ymax=OR.U95, x=coefficient, color=metaGRS) +
  geom_hline(yintercept=1, linetype=2) +
  geom_errorbar(width=0, alpha=0.8) +
  geom_point(shape=19) +
  scale_color_manual(values=c("current"="#ae017e", "previous"="#000000", "no"="#045a8d")) +
  xlab("") +
  ylab("Odds Ratio (95% CI)") +
  theme_bw() +
  theme(
    legend.position="none",
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=0.5),
    axis.text.y=element_text(size=8),
    axis.title=element_text(size=10)
  )

g_cind <- ggplot(inci_pgs) +
  aes(y=C.index, ymin=C.L95, ymax=C.U95, x=coefficient, color=metaGRS) +
	geom_rect(inherit.aes=FALSE, data=inci_null,
					aes(ymin=C.L95, ymax=C.U95),
					xmin=-Inf, xmax=Inf, fill="#bdbdbd", alpha=0.3) +
  geom_hline(data=inci_null, aes(yintercept=C.index), linetype=2) +
  geom_errorbar(width=0, alpha=0.8) +
  geom_point(shape=19) +
  scale_color_manual(values=c("current"="#ae017e", "previous"="#000000", "no"="#045a8d")) +
  xlab("") +
  ylab("C-index (95% CI)") +
  theme_bw() +
  theme(
    legend.position="none",
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=0.5),
    axis.text.y=element_text(size=8),
    axis.title=element_text(size=10)
  )

g_hr <- ggplot(inci_pgs) +
  aes(y=HR, ymin=L95, ymax=U95, x=coefficient, color=metaGRS) +
  geom_hline(yintercept=1, linetype=2) +
  geom_errorbar(width=0, alpha=0.8) +
  geom_point(shape=19) +
  scale_color_manual(values=c("current"="#ae017e", "previous"="#000000", "no"="#045a8d")) +
  xlab("") +
  ylab("Hazard Ratio (95% CI)") +
  theme_bw() +
  theme(
    legend.position="none",
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=0.5),
    axis.text.y=element_text(size=8),
    axis.title=element_text(size=10)
  )

g <- plot_grid(g_auc, g_cind, g_or, g_hr, nrow=2, align="hv")
ggsave(g, width=7.2, height=7.2, units="in", file="output/slide_deck_plots/2023-02-10-highlights/UKB_pgs_comparison.png")

# Compare QDiabetes scores
qdiab <- fread("output/UKB_tests/incident_T2D_associations.txt")
qdiab <- qdiab[coefficient %like% "QDiabetes"]
qdiab[model == "QDiabetes2018C", coefficient := "QDiabetes 2018 model C (HbA1c)"]
qdiab <- qdiab[c(5, 1:3, 6:7, 4)]
qdiab[, coefficient := factor(coefficient, levels=rev(unique(coefficient)))]

g_cind <- ggplot(qdiab) +
  aes(x=C.index, xmin=C.L95, xmax=C.U95, y=coefficient) +
	geom_rect(inherit.aes=FALSE, data=inci_null,
					aes(xmin=C.L95, xmax=C.U95),
					ymin=-Inf, ymax=Inf, fill="#bdbdbd", alpha=0.3) +
  geom_vline(data=inci_null, aes(xintercept=C.index), linetype=2) +
  geom_errorbarh(height=0, alpha=0.8, color="#e31a1c") +
  geom_point(shape=19, color="#e31a1c") +
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

g_hr <- ggplot(qdiab) +
  aes(x=HR, xmin=L95, xmax=U95, y=coefficient) +
  geom_vline(xintercept=1, linetype=2) +
  geom_errorbarh(height=0, alpha=0.8, color="#e31a1c") +
  geom_point(shape=19, color="#e31a1c") +
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

g <- plot_grid(g_cind, g_hr, nrow=1, align="hv")
ggsave(g, width=13, height=3, units="in", file="output/slide_deck_plots/2023-02-10-highlights/UKB_qdiab_comparison.png")

# Plot again in QDiabetes model C subset
qdiab <- fread("output/UKB_tests/incident_T2D_associations_QDiabetes2018C_subset.txt")
qdiab <- qdiab[coefficient %like% "QDiabetes"]
qdiab[model == "QDiabetes2018C", coefficient := "QDiabetes 2018 model C (HbA1c)"]
qdiab <- qdiab[c(5, 1:3, 6:7, 4)]
qdiab[, coefficient := factor(coefficient, levels=rev(unique(coefficient)))]

g_cind <- ggplot(qdiab) +
  aes(x=C.index, xmin=C.L95, xmax=C.U95, y=coefficient) +
	geom_rect(inherit.aes=FALSE, data=inci_null,
					aes(xmin=C.L95, xmax=C.U95),
					ymin=-Inf, ymax=Inf, fill="#bdbdbd", alpha=0.3) +
  geom_vline(data=inci_null, aes(xintercept=C.index), linetype=2) +
  geom_errorbarh(height=0, alpha=0.8, color="#e31a1c") +
  geom_point(shape=19, color="#e31a1c") +
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

g_hr <- ggplot(qdiab) +
  aes(x=HR, xmin=L95, xmax=U95, y=coefficient) +
  geom_vline(xintercept=1, linetype=2) +
  geom_errorbarh(height=0, alpha=0.8, color="#e31a1c") +
  geom_point(shape=19, color="#e31a1c") +
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

g <- plot_grid(g_cind, g_hr, nrow=1, align="hv")
ggsave(g, width=13, height=3, units="in", file="output/slide_deck_plots/2023-02-10-highlights/UKB_qdiab_comparison_modelC_subset.png")

# Plot C-index for GRS vs. QDiabetes and risk factors
inci <- fread("output/UKB_tests/incident_T2D_associations_QDiabetes2018C_subset.txt")
fit <- inci[, .SD[1], by=.(model_type, model)]
fit <- fit[!(model_type == "reference") | (model %like% " \\+ follow_diff")]
fit <- fit[!(model %like% "strata")]
fit <- fit[!(model_type == "QDiabetes") | (model == "QDiabetes2018C")]
fit <- fit[model != "fasting_glucose"]
fit[model == "age + sex + follow_diff", coefficient := "Age + Sex"]
fit[model == "smoking_status", coefficient := "Smoking Status"]
fit[, model := coefficient]

multi <- fread("output/UKB_tests/incident_T2D_multivariate_associations.txt")
multi <- multi[,.SD[1],by=.(model_type, model)]
multi <- multi[coefficient == "T2D_metaGRS"]
multi[, model := gsub("PGS", "T2D_metaGRS", model_type)]
multi[, model_type := "multivariate"]

fit <- rbind(fit, multi)

fit <- fit[!(coefficient %like% "Gad")]

fit <- fit[order(-C.index)]
fit[, model := factor(model, levels=unique(model))]
fit[, model_type := fcase(
  model_type == "reference", "Baseline characteristics",
  model_type == "PGS", "Polygenic risk scores",
  model_type == "risk factor", "QDiabetes 2018 risk factors",
  model_type == "QDiabetes", "QDiabetes 2018 model C",
  model_type == "multivariate", "Multivariate"
)]
fit[, model_type := factor(model_type, levels=c(
  "Baseline characteristics", "Polygenic risk scores", "QDiabetes 2018 risk factors",
  "QDiabetes 2018 model C", "Multivariate" 
))]

fit_ref <- fit[model == "Age + Sex", .(C.index, C.L95, C.U95)]
fit_metaGRS <- fit[model == "T2D_metaGRS", .(C.index, C.L95, C.U95)]
fit_qdiab <- fit[model == "QDiabetes 2018 model C", .(C.index, C.L95, C.U95)] 

g <- ggplot(fit) +
  aes(y=model, x=C.index, xmin=C.L95, xmax=C.U95) +
  geom_rect(inherit.aes=FALSE, data=fit_ref, aes(xmin=C.L95, xmax=C.U95),
          ymin=-Inf, ymax=Inf, fill="#bdbdbd", alpha=0.3) +
  geom_vline(data=fit_ref, aes(xintercept=C.index), linetype=2) +
  geom_rect(inherit.aes=FALSE, data=fit_metaGRS, aes(xmin=C.L95, xmax=C.U95),
          ymin=-Inf, ymax=Inf, fill="#6a51a3", alpha=0.3) +
  geom_vline(data=fit_metaGRS, aes(xintercept=C.index), linetype=2, color="#3f007d") +
  geom_rect(inherit.aes=FALSE, data=fit_qdiab, aes(xmin=C.L95, xmax=C.U95),
          ymin=-Inf, ymax=Inf, fill="#fc4e2a", alpha=0.3) +
  geom_vline(data=fit_qdiab, aes(xintercept=C.index), linetype=2, color="#e31a1c") +
  geom_errorbarh(height=0) +
  geom_point(shape=21, fill="#bdbdbd", color="black") +
  facet_col(facets = vars(model_type), scales="free_y", space="free") +
  xlab("C-index (95% CI)") +
  ylab("") + 
  theme_bw() +
  theme(
    panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    axis.text.y=element_text(size=8, color="black"),
    axis.text.x=element_text(size=8),
    axis.title=element_text(size=10),
    strip.text=element_text(size=8, face=2), 
    strip.background=element_blank()
  )
ggsave(g, width=7.2, height=7.2, units="in", file="output/slide_deck_plots/2023-02-10-highlights/UKB_prs_risk_factor_cindex_comparison.png")

# Show hazard ratios
inci <- fread("output/UKB_tests/incident_T2D_associations_QDiabetes2018C_subset.txt")

age_sex <- inci[model == "age + sex + follow_diff"][1:2]
age_sex[, model := coefficient]
uni <- inci[model_type != "reference", .SD[1], by=.(model)]
uni <- uni[!(model_type == "QDiabetes") | (model == "QDiabetes2018C")]
uni <- uni[model != "fasting_glucose" & model != "smoking_status"]
uni[, model := coefficient]
smoking <- inci[model == "smoking_status" & coefficient %like% "Smoking"]
smoking[, model := "Smoking status"]
uni <- rbind(age_sex, uni, smoking)
uni[, model_type := "Univariate"]

multi <- fread("output/UKB_tests/incident_T2D_multivariate_associations.txt")
multi <- multi[model == "T2D_metaGRS"]
multi <- multi[!(coefficient %like% "Assessment" | coefficient %like% "hospital" | coefficient %like% "Sex" | coefficient == "Age")]
multi[, model := coefficient]
multi[, model_type := gsub("PGS", "T2D_metaGRS", model_type)]

hrs <- rbind(uni, multi)

hrs <- hrs[!(model %like% "Gad")]

hrs <- hrs[order(-abs(logHR))]

hrs[, model_type := factor(model_type, levels=c("Univariate", "T2D_metaGRS + BMI + HbA1c", "T2D_metaGRS + QDiabetes2018C"))]
hrs[, model := factor(model, levels=unique(model))]

hrs_metaGRS <- hrs[model_type == "Univariate" & coefficient == "T2D_metaGRS", .(HR, L95, U95)]
hrs_qdiab <- hrs[model_type == "Univariate" & coefficient == "QDiabetes 2018 model C", .(HR, L95, U95)]

g <- ggplot(hrs) +
  aes(y=model, x=HR, xmin=L95, xmax=U95, by=coefficient) +
  geom_vline(xintercept=1, linetype=2) +
  geom_rect(inherit.aes=FALSE, data=hrs_metaGRS, aes(xmin=L95, xmax=U95),
          ymin=-Inf, ymax=Inf, fill="#6a51a3", alpha=0.3) +
  geom_vline(data=hrs_metaGRS, aes(xintercept=HR), linetype=2, color="#3f007d") +
  geom_rect(inherit.aes=FALSE, data=hrs_qdiab, aes(xmin=L95, xmax=U95),
          ymin=-Inf, ymax=Inf, fill="#fc4e2a", alpha=0.3) +
  geom_vline(data=hrs_qdiab, aes(xintercept=HR), linetype=2, color="#e31a1c") +
  geom_errorbarh(height=0, position=position_dodge(width=0.8)) +
  geom_point(shape=21, fill="#bdbdbd", color="black", position=position_dodge(width=0.8)) +
  facet_col(facets = vars(model_type), scales="free_y", space="free") +
  scale_x_sqrt("Hazard Ratio (95% CI)", breaks=1:10) +
  ylab("") + 
  coord_cartesian(xlim = c(0.5, 10)) +
  theme_bw() +
  theme(
    panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    axis.text.y=element_text(size=8, color="black"),
    axis.text.x=element_text(size=8),
    axis.title=element_text(size=10),
    strip.text=element_text(size=8, face=2), 
    strip.background=element_blank()
  )
ggsave(g, width=7.2, height=7.2, units="in", file="output/slide_deck_plots/2023-02-10-highlights/UKB_prs_risk_factor_hazard_ratio_comparison.png")

