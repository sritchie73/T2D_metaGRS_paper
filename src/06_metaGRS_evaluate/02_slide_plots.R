library(data.table)
library(ggplot2)
library(cowplot)

# Make output directory
system("mkdir -p output/slide_deck_plots/", wait=TRUE)

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
ggsave(g, width=7.2, height=7.2, units="in", file="output/slide_deck_plots/UKB_pgs_comparison.png")

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
ggsave(g, width=13, height=3, units="in", file="output/slide_deck_plots/UKB_qdiab_comparison.png")

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
ggsave(g, width=13, height=3, units="in", file="output/slide_deck_plots/UKB_qdiab_comparison_modelC_subset.png")
