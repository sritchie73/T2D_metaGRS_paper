library(data.table)
library(ggplot2)
library(cowplot)

# Make output directory
system("mkdir -p output/slide_deck_plots/", wait=TRUE)

# Load associations across different ancestry groups
prev <- fread("output/UKB_tests/non_european_prevalent_T2D_associations.txt")
inci <- fread("output/UKB_tests/non_european_incident_T2D_associations.txt")

# Get null model associations
prev_null <- prev[model == "age + sex", .SD[1], by=.(grouping, group_column)]
inci_null <- inci[model == "age + sex", .SD[1], by=.(grouping, group_column)]

# Extract PGS associations and modle fit
prev_pgs <- prev[coefficient %like% "PGS" | coefficient %like% "GRS"]
inci_pgs <- inci[coefficient %like% "PGS" | coefficient %like% "GRS"]

# Order:
prev_pgs <- prev_pgs[order(-AUC)]
prev_pgs[, AUC_rank := 1:.N, by=.(grouping, group_column)]
pgs_rank <- prev_pgs[,.(ranksum=sum(AUC_rank)),by=.(coefficient)]
prev_pgs <- prev_pgs[pgs_rank[order(ranksum), .(coefficient)], on = .(coefficient)]
prev_pgs[, coefficient := factor(coefficient, levels=unique(coefficient))]

inci_pgs <- inci_pgs[order(-C.index)]
inci_pgs[, Cindex_rank := 1:.N, by=.(grouping, group_column)]
pgs_rank <- inci_pgs[,.(ranksum=sum(Cindex_rank)),by=.(coefficient)]
inci_pgs <- inci_pgs[pgs_rank[order(ranksum), .(coefficient)], on = .(coefficient)]
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

# Compare PGS across ancestry groups
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
  facet_wrap(~ group_column + grouping, scales="free_y", nrow=2) +
  theme_bw() +
  theme(
    legend.position="none",
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=0.5),
    axis.text.y=element_text(size=8),
    axis.title=element_text(size=10),
    strip.text=element_text(size=8)
  )
ggsave(g_auc, width=13.5, height=6, units="in", file="output/slide_deck_plots/UKB_non_european_PGS_AUC_compare.png")

g_or <- ggplot(prev_pgs) +
  aes(y=OR, ymin=OR.L95, ymax=OR.U95, x=coefficient, color=metaGRS) +
  geom_hline(yintercept=1, linetype=2) +
  geom_errorbar(width=0, alpha=0.8) +
  geom_point(shape=19) +
  scale_color_manual(values=c("current"="#ae017e", "previous"="#000000", "no"="#045a8d")) +
  xlab("") +
  ylab("Odds Ratio (95% CI)") +
  facet_wrap(~ group_column + grouping, scales="free_y", nrow=2) +
  theme_bw() +
  theme(
    legend.position="none",
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=0.5),
    axis.text.y=element_text(size=8),
    axis.title=element_text(size=10),
    strip.text=element_text(size=8)
  )
ggsave(g_auc, width=13.5, height=6, units="in", file="output/slide_deck_plots/UKB_non_european_PGS_OR_compare.png")

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
  facet_wrap(~ group_column + grouping, scales="free_y", nrow=2) +
  theme_bw() +
  theme(
    legend.position="none",
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=0.5),
    axis.text.y=element_text(size=8),
    axis.title=element_text(size=10),
    strip.text=element_text(size=8)
  )
ggsave(g_auc, width=13.5, height=6, units="in", file="output/slide_deck_plots/UKB_non_european_PGS_Cindex_compare.png")

g_or <- ggplot(inci_pgs) +
  aes(y=HR, ymin=L95, ymax=U95, x=coefficient, color=metaGRS) +
  geom_hline(yintercept=1, linetype=2) +
  geom_errorbar(width=0, alpha=0.8) +
  geom_point(shape=19) +
  scale_color_manual(values=c("current"="#ae017e", "previous"="#000000", "no"="#045a8d")) +
  xlab("") +
  ylab("Hazard Ratio (95% CI)") +
  facet_wrap(~ group_column + grouping, scales="free_y", nrow=2) +
  theme_bw() +
  theme(
    legend.position="none",
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=0.5),
    axis.text.y=element_text(size=8),
    axis.title=element_text(size=10),
    strip.text=element_text(size=8)
  )
ggsave(g_auc, width=13.5, height=6, units="in", file="output/slide_deck_plots/UKB_non_european_PGS_HR_compare.png")
