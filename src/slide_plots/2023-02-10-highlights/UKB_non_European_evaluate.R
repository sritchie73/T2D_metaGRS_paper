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

# Drop older versions of metaGRS
prev_pgs <- prev_pgs[!(coefficient %like% "Gad")]
inci_pgs <- inci_pgs[!(coefficient %like% "Gad")]

# Filter to groupings we want to keep for simplified plots
prev_pgs <- prev_pgs[group_column == "ethnicity_grouping" | grouping == "Other ethnic group"]
inci_pgs <- inci_pgs[group_column == "ethnicity_grouping" | grouping == "Other ethnic group"]
prev_pgs[grouping == "Other ethnic group", grouping := "Other non-European"]
inci_pgs[grouping == "Other ethnic group", grouping := "Other non-European"]

prev_null <- prev_null[group_column == "ethnicity_grouping" | grouping == "Other ethnic group"]
inci_null <- inci_null[group_column == "ethnicity_grouping" | grouping == "Other ethnic group"]
prev_null[grouping == "Other ethnic group", grouping := "Other non-European"]
inci_null[grouping == "Other ethnic group", grouping := "Other non-European"]

# Flag metaGRS (for plot colours)
prev_pgs[, metaGRS := fcase(
  coefficient == "T2D_metaGRS", "current",
  coefficient %like% "Gad", "previous",
  default = "no")]

inci_pgs[, metaGRS := fcase(
  coefficient == "T2D_metaGRS", "current",
  coefficient %like% "Gad", "previous",
  default = "no")]

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

prev_pgs[, grouping := factor(grouping, levels=c("South Asian", "East Asian", "Other non-European", "African"))]
inci_pgs[, grouping := factor(grouping, levels=c("South Asian", "East Asian", "Other non-European", "African"))]
prev_null[, grouping := factor(grouping, levels=c("South Asian", "East Asian", "Other non-European", "African"))]
inci_null[, grouping := factor(grouping, levels=c("South Asian", "East Asian", "Other non-European", "African"))]

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
  facet_wrap(~ grouping, scales="free_y", nrow=1) +
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
ggsave(g_auc, width=9, height=4, units="in", file="output/slide_deck_plots/UKB_non_european_PGS_AUC_compare.png")

g_or <- ggplot(prev_pgs) +
  aes(y=OR, ymin=OR.L95, ymax=OR.U95, x=coefficient, color=metaGRS) +
  geom_hline(yintercept=1, linetype=2) +
  geom_errorbar(width=0, alpha=0.8) +
  geom_point(shape=19) +
  scale_color_manual(values=c("current"="#ae017e", "previous"="#000000", "no"="#045a8d")) +
  xlab("") +
  ylab("Odds Ratio (95% CI)") +
  facet_wrap(~ grouping, scales="free_y", nrow=1) +
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
ggsave(g_or, width=9, height=4, units="in", file="output/slide_deck_plots/UKB_non_european_PGS_OR_compare.png")

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
  facet_wrap(~ grouping, scales="free_y", nrow=1) +
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
ggsave(g_cind, width=9, height=4, units="in", file="output/slide_deck_plots/UKB_non_european_PGS_Cindex_compare.png")

g_hr <- ggplot(inci_pgs) +
  aes(y=HR, ymin=L95, ymax=U95, x=coefficient, color=metaGRS) +
  geom_hline(yintercept=1, linetype=2) +
  geom_errorbar(width=0, alpha=0.8) +
  geom_point(shape=19) +
  scale_color_manual(values=c("current"="#ae017e", "previous"="#000000", "no"="#045a8d")) +
  xlab("") +
  ylab("Hazard Ratio (95% CI)") +
  facet_wrap(~ grouping, scales="free_y", nrow=1) +
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
ggsave(g_hr, width=9, height=4, units="in", file="output/slide_deck_plots/UKB_non_european_PGS_HR_compare.png")

# Examine transferrability.

# Filter to T2D metaGRS
prev_pgs <- prev_pgs[coefficient == "T2D_metaGRS"]
inci_pgs <- inci_pgs[coefficient == "T2D_metaGRS"]

# Drop extraneous columns
prev_pgs[, AUC_rank := NULL]
prev_pgs[, metaGRS := NULL]

inci_pgs[, Cindex_rank := NULL]
inci_pgs[, metaGRS := NULL]

# Add in White British cohort results
prev_wb <- fread("output/UKB_tests/prevalent_T2D_associations.txt")
prev_wb <- prev_wb[coefficient == "T2D_metaGRS"]
prev_wb[, model_type := NULL]
prev_wb[, grouping := "White British"]
prev_wb[, group_column := "metaGRS_test_set"]
prev_pgs <- rbind(prev_pgs, prev_wb)

inci_wb <- fread("output/UKB_tests/incident_T2D_associations.txt")
inci_wb <- inci_wb[coefficient == "T2D_metaGRS"]
inci_wb[, model_type := NULL]
inci_wb[, grouping := "White British"]
inci_wb[, group_column := "metaGRS_test_set"]
inci_pgs <- rbind(inci_pgs, inci_wb)

# Relabel so axes are unique
##' prev_pgs[, xlabel := fcase(
##' 	grouping == "White British" & group_column == "metaGRS_test_set", "White British (test set)",
##' 	grouping == "South Asian" & group_column == "ethnicity_grouping", "South Asian (supergroup)",
##' 	grouping == "Indian" & group_column == "ethnicity_subgroup", "Indian",
##' 	grouping == "Bangladeshi" & group_column == "ethnicity_subgroup", "Bangladeshi",
##' 	grouping == "Pakistani" & group_column == "ethnicity_subgroup", "Pakistani",
##' 	grouping == "East Asian" & group_column == "ethnicity_grouping", "East Asian (supergroup)",
##' 	grouping == "Chinese" & group_column == "ethnicity_subgroup",  "Chinese",
##' 	grouping == "Any other Asian background" & group_column == "ethnicity_subgroup", "Any other Asian background",
##' 	grouping == "African" & group_column == "ethnicity_grouping", "African (supergroup)",
##' 	grouping == "African" & group_column == "ethnicity_subgroup", "African",
##' 	grouping == "Caribbean" & group_column == "ethnicity_subgroup", "Caribbean",
##' 	grouping == "Other ethnic group" & group_column == "ethnicity_subgroup", "Other ethnic group",
##' 	grouping == "Other" & group_column == "ethnicity_qdiabetes_groups", "Other (QDiabetes supergroup)",
##' 	grouping == "OtherAsian" & group_column == "ethnicity_qdiabetes_groups", "Other Asian (QDiabetes supergroup)"
##' )]
##' prev_pgs[, xlabel := factor(xlabel, levels=c(
##' 	"White British (test set)",
##'   "South Asian (supergroup)", "Indian", "Bangladeshi", "Pakistani",
##'   "East Asian (supergroup)", "Chinese", "Any other Asian background",
##'   "African (supergroup)", "African", "Caribbean",
##'   "Other ethnic group",
##'   "Other (QDiabetes supergroup)",
##'   "Other Asian (QDiabetes supergroup)"
##' ))]
##' 
##' inci_pgs[, xlabel := fcase(
##' 	grouping == "White British" & group_column == "metaGRS_test_set", "White British (test set)",
##' 	grouping == "South Asian" & group_column == "ethnicity_grouping", "South Asian (supergroup)",
##' 	grouping == "Indian" & group_column == "ethnicity_subgroup", "Indian",
##' 	grouping == "Bangladeshi" & group_column == "ethnicity_subgroup", "Bangladeshi",
##' 	grouping == "Pakistani" & group_column == "ethnicity_subgroup", "Pakistani",
##' 	grouping == "East Asian" & group_column == "ethnicity_grouping", "East Asian (supergroup)",
##' 	grouping == "Chinese" & group_column == "ethnicity_subgroup",  "Chinese",
##' 	grouping == "Any other Asian background" & group_column == "ethnicity_subgroup", "Any other Asian background",
##' 	grouping == "African" & group_column == "ethnicity_grouping", "African (supergroup)",
##' 	grouping == "African" & group_column == "ethnicity_subgroup", "African",
##' 	grouping == "Caribbean" & group_column == "ethnicity_subgroup", "Caribbean",
##' 	grouping == "Other ethnic group" & group_column == "ethnicity_subgroup", "Other ethnic group",
##' 	grouping == "Other" & group_column == "ethnicity_qdiabetes_groups", "Other (QDiabetes supergroup)",
##' 	grouping == "OtherAsian" & group_column == "ethnicity_qdiabetes_groups", "Other Asian (QDiabetes supergroup)"
##' )]
##' inci_pgs[, xlabel := factor(xlabel, levels=c(
##' 	"White British (test set)",
##'   "South Asian (supergroup)", "Indian", "Bangladeshi", "Pakistani",
##'   "East Asian (supergroup)", "Chinese", "Any other Asian background",
##'   "African (supergroup)", "African", "Caribbean",
##'   "Other ethnic group",
##'   "Other (QDiabetes supergroup)",
##'   "Other Asian (QDiabetes supergroup)"
##' ))]
prev_pgs[, xlabel := factor(grouping, levels=c("White British", "South Asian", "East Asian", "Other non-European", "African"))]
inci_pgs[, xlabel := factor(grouping, levels=c("White British", "South Asian", "East Asian", "Other non-European", "African"))]

# Build plots
g_or <- ggplot(prev_pgs) +
  aes(y=OR, ymin=OR.L95, ymax=OR.U95, x=xlabel) +
  geom_hline(yintercept=1, linetype=2) +
  geom_errorbar(width=0, alpha=0.8, color="#ae017e") +
  geom_point(shape=19, color="#ae017e") +
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

g_hr <- ggplot(inci_pgs) +
  aes(y=HR, ymin=L95, ymax=U95, x=xlabel) +
  geom_hline(yintercept=1, linetype=2) +
  geom_errorbar(width=0, alpha=0.8, color="#ae017e") +
  geom_point(shape=19, color="#ae017e") +
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

g <- plot_grid(g_or, g_hr, align="hv")
ggsave(g, width=4, height=7.2, units="in", file="output/slide_deck_plots/UKB_T2D_metaGRS_transferrability.png")


