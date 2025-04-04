library(data.table)
library(ggplot2)
library(ggstance)
library(cowplot)

system("mkdir -p output/paper_outputs/{figures,tables}", wait=TRUE)

# Load C-indice information and prepare for plotting
cinds <- fread("output/UKB_tests/QDiabetes_comparison_by_age_tertile.txt")
cinds <- cinds[,.SD[1],by=.(age_tertile, model_type)] # multiple rows per model - HRs also reported
cinds[, QDiabetes_model := fcase(
  model_type %like% "^T2D_metaGRS", "T2D metaPRS",
  model_type %like% "2018A", "QDiabetes A",
  model_type %like% "2018B", "QDiabetes B",
  model_type %like% "2018C", "QDiabetes C"
)]
cinds[, model := fcase(
  coefficient %like% "QDiabetes", "QDiabetes",
  coefficient == "T2D_metaGRS" & model_type %like% "^T2D_metaGRS", "T2D metaPRS", 
  default = "QDiabetes + T2D metaPRS"
)]
cinds[, model := factor(model, levels=c("QDiabetes + T2D metaPRS", "QDiabetes", "T2D metaPRS"))]
cinds[, age_tertile := factor(age_tertile, levels=c("(61,72]", "(53,61]", "[39,53]"))]
cinds <- cinds[!(model_type %in% c("T2D_metaGRS (QDiabetes 2018A subcohort)", "T2D_metaGRS (QDiabetes 2018B subcohort)"))]

# Load delta C-index information and preprare for plotting
deltaC <- fread("output/UKB_tests/age_tertiles_delta_cind_bootstraps.txt")
deltaC <- deltaC[prs_added_to != "QDiabetes2013"]
deltaC[, QDiabetes_model := fcase(
  prs_added_to %like% "2018A", "QDiabetes A",
  prs_added_to %like% "2018B", "QDiabetes B",
  prs_added_to %like% "2018C", "QDiabetes C"
)]
deltaC[, age_tertile := factor(age_tertile, levels=c("(61,72]", "(53,61]", "[39,53]"))]

# Make plots
g1 <- ggplot(cinds) +
  aes(x=C.index, xmin=C.L95, xmax=C.U95, y=age_tertile, color=model) +
  facet_wrap(~ QDiabetes_model, nrow=1) +
  geom_errorbarh(height=0, position=position_dodgev(heigh=0.6)) +
  geom_point(shape=23, fill="white", position=position_dodgev(heigh=0.6)) +
  scale_color_manual(values=c("QDiabetes + T2D metaPRS"="#bd0026", "QDiabetes"="#08306b", "T2D metaPRS"="#810f7c")) +
  xlab("C-index (95% CI)") +
  ylab("Age tertile") +
  theme_bw() +
  theme(
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
    axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
    strip.background=element_blank(), strip.text=element_text(size=7, face="bold"),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="bottom", legend.text=element_text(size=7), legend.title=element_blank()
  )

g2 <- ggplot(deltaC) + 
  aes(x=delta_C, xmin=boot_L95, xmax=boot_U95, y=age_tertile) +
  facet_wrap(~ QDiabetes_model) +
  geom_vline(xintercept=0, linetype=2) +
  geom_errorbarh(height=0, color="#bd0026") +
  geom_point(shape=23, fill="white", color="#bd0026") +
  xlab("Change in C-index (95% CI)") +
  ylab("Age tertile") +
  theme_bw() +
  theme(
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
    axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
    strip.background=element_blank(), strip.text=element_text(size=7, face="bold"),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank()
  )

g <- plot_grid(g1, g2, rel_heights=c(0.6, 0.4), align="hv", axis="r", nrow=2)

ggsave(g, width=7.2, height=4, file="output/paper_outputs/figures/age_tertile_analysis.pdf")

# Build output table for supp
dt <- fread("output/UKB_tests/QDiabetes_comparison_by_age_tertile.txt")
dt[, model_type := gsub("2018", " 2018 model ", model_type)]
dt[, model_type := gsub("T2D_metaGRS", "T2D metaPRS", model_type)]
dt[, coefficient := gsub("T2D_metaGRS", "T2D metaPRS", coefficient)]
dt[, age_tertile := fcase(
  age_tertile == "[39,53]", "39-53",
  age_tertile == "(53,61]", "54-61",
  age_tertile == "(61,72]", "62-72"
)]
dt <- dt[,.(model=model_type, age_tertile, Samples, Cases, C.index, C.L95, C.U95, coefficient, HR, L95, U95, Pvalue)]
dt <- dt[order(model)]

deltaC[, model_type := fcase(
  prs_added_to == "QDiabetes2018A", "QDiabetes 2018 model A + T2D metaPRS",
  prs_added_to == "QDiabetes2018B_non_fasting", "QDiabetes 2018 model B + T2D metaPRS",
  prs_added_to == "QDiabetes2018C", "QDiabetes 2018 model C + T2D metaPRS"
)]
deltaC[,age_tertile := fcase(
  age_tertile == "[39,53]", "39-53",
  age_tertile == "(53,61]", "54-61",
  age_tertile == "(61,72]", "62-72"
)]
dt[deltaC, on = .(model=model_type, age_tertile), c("deltaC", "deltaC.L95", "deltaC.U95", "deltaC.Pval") := .(delta_C, boot_L95, boot_U95, boot_pval)]

fwrite(dt, sep="\t", quote=FALSE, file="output/paper_outputs/tables/age_tertile_qdiabetes.txt")



  
