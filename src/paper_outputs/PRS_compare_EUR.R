library(data.table)
library(ggplot2)
library(patchwork)

system("mkdir -p output/paper_outputs/{figures,tables}", wait=TRUE)

# Create unified theme for plots
theme_prs <- function(...) {
  theme_bw() +
  theme(...) %+replace% theme(
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
    axis.text.y=element_text(size=6, color="black"), axis.title.y=element_blank(),
    plot.title=element_text(size=7, face="bold"), plot.title.position="plot",
    plot.subtitle=element_text(size=6, color="#525252"),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="none"
  )
}

# Compare PRS performance for prevalent T2D stratification in UKB
prev_ukb <- fread("output/UKB_tests/prevalent_T2D_associations.txt")
prev_ukb <- prev_ukb[model == coefficient & !(model %like% "^Gad")]
prev_ukb <- prev_ukb[model != "Mars2022_AJHG_PGS002771"] # Developed based on Mahajan 2018 GWAS in UKB - hence odds ratio of 4.3!!
prev_ukb <- prev_ukb[order(OR)]
prev_ukb[, model := factor(model, levels=model)]
prev_ukb[, color_anno := model == "T2D_metaGRS"]

g1 <- ggplot(prev_ukb) +
  aes(x=OR, xmin=OR.L95, xmax=OR.U95, y=model, color=color_anno) +
  geom_vline(xintercept=1, linetype=2) +
  geom_point(shape=23, fill="white") + 
  geom_errorbarh(height=0) +
  scale_color_manual(values=c("TRUE"="#bd0026", "FALSE"="#08306b")) +
  scale_x_continuous("Odds Ratio (95% CI)", limits=c(1, 2.5), breaks=c(1,1.5,2,2.5)) +
  ggtitle("Prediction of T2D prevalence in UK Biobank", subtitle="N=245,177; 11,080 prevalent T2D cases") +
  theme_prs()

# Compare PRS performance for incident T2D risk prediction in UKB
inci_ukb <- fread("output/UKB_tests/incident_T2D_associations.txt")
inci_ukb <- inci_ukb[model == coefficient & !(model %like% "^Gad")]
inci_ukb <- inci_ukb[model != "Mars2022_AJHG_PGS002771"] # Developed based on Mahajan 2018 GWAS in UKB
inci_ukb <- inci_ukb[order(HR)]
inci_ukb[, model := factor(model, levels=model)]
inci_ukb[, color_anno := model == "T2D_metaGRS"]

g2 <- ggplot(inci_ukb) +
  aes(x=HR, xmin=L95, xmax=U95, y=model, color=color_anno) +
  geom_vline(xintercept=1, linetype=2) +
  geom_point(shape=23, fill="white") + 
  geom_errorbarh(height=0) +
  scale_color_manual(values=c("TRUE"="#bd0026", "FALSE"="#08306b")) +
  scale_x_continuous("Hazard Ratio (95% CI)", limits=c(1, 2.5), breaks=c(1,1.5,2,2.5)) +
  ggtitle("Prediction of incident T2D risk in UK Biobank", subtitle="N=232,808; 6,016 T2D cases within 10 years of follow-up") +
  theme_prs()

# Compare PRS performance for incident T2D risk prediction in INTERVAL
inci_interval <- fread("output/INTERVAL_tests/incident_diabetes_associations.txt")
inci_interval <- inci_interval[model == coefficient & !(model %like% "^Gad")]
inci_interval <- inci_interval[order(HR)]
inci_interval[, model := factor(model, levels=model)]
inci_interval[, color_anno := model == "T2D_metaGRS"]

g3 <- ggplot(inci_interval) +
  aes(x=HR, xmin=L95, xmax=U95, y=model, color=color_anno) +
  geom_vline(xintercept=1, linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, fill="white") + 
  scale_color_manual(values=c("TRUE"="#bd0026", "FALSE"="#08306b")) +
  scale_x_continuous("Hazard Ratio (95% CI)", limits=c(1, 2.5), breaks=c(1,1.5,2,2.5)) +
  ggtitle("Prediction of incident T2D risk in INTERVAL", subtitle="N=38,941; 726 T2D cases within 11 years of follow-up") +
  theme_prs()

# Compare PRS performance for T2D prediction in All of Us
cc_aou <- fread("data/Henry_Taylor_AoU_cohort/glm_results__pgs_covariates.tsv")
cc_aou <- cc_aou[ancestry == "eur"]

prs_map <- cc_aou[,.(AoU_name=grs)]
prs_map[, PGS_ID := gsub(".*_", "", AoU_name)]
prs_map[PGS_ID == "chagoya", PGS_ID := "PGS00344X"]
prs_map2 <- inci_interval[,.(name=model)]
prs_map2[, PGS_ID := gsub(".*_", "", name)]
prs_map <- merge(prs_map, prs_map2, by="PGS_ID", all=TRUE)
cc_aou[prs_map, on = .(grs=AoU_name), model := i.name]

cc_aou <- cc_aou[order(OR)]
cc_aou[, model := factor(model, levels=model)]
cc_aou[, color_anno := model == "T2D_metaGRS"]

g4 <- ggplot(cc_aou) +
  aes(x=OR, xmin=OR.L95, xmax=OR.U95, y=model, color=color_anno) +
  geom_vline(xintercept=1, linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, fill="white") +
  scale_color_manual(values=c("TRUE"="#bd0026", "FALSE"="#08306b")) +
  scale_x_continuous("Odds Ratio (95% CI)", limits=c(1, 2.5), breaks=c(1,1.5,2,2.5)) +
  ggtitle("Prediction of T2D prevalence in All of Us", subtitle="N=109,021; 10,069 prevalent T2D cases") +
  theme_prs()

# Arrange into single plot
g <- (g1 / g2) | g3 | g4
ggsave(g, width=7.2, height=3.5, file="output/paper_outputs/figures/PRS_compare_EUR.pdf")

# Write out tables
prev_dt <- rbind(idcol="cohort",
  "UK Biobank"=prev_ukb[,.(ancestry="EUR", samples, cases, PRS=model, OR, OR.L95, OR.U95, P.value, AUC, AUC.L95, AUC.U95)],
  "All of Us"=cc_aou[ancestry == "eur", .(ancestry="EUR", samples, cases, PRS=model, OR, OR.L95, OR.U95, P.value, AUC, AUC.L95, AUC.U95)]
)
prev_dt <- prev_dt[order(-OR)][order(cohort)]

inci_dt <- rbind(idcol="cohort",
  "UK Biobank"=inci_ukb[,.(ancestry="EUR", samples=Samples, cases=Cases, PRS=model, HR, HR.L95=L95, HR.U95=U95, P.value=Pvalue, C.index, C.L95, C.U95)],
  "INTERVAL"=inci_interval[,.(ancestry="EUR", samples=Samples, cases=Cases, PRS=model, HR, HR.L95=L95, HR.U95=U95, P.value=Pvalue, C.index, C.L95, C.U95)]
)
inci_dt <- inci_dt[order(-HR)][order(cohort)]

fwrite(prev_dt, sep="\t", quote=FALSE, file="output/paper_outputs/tables/PRS_compare_EUR_prevalent.txt")
fwrite(inci_dt, sep="\t", quote=FALSE, file="output/paper_outputs/tables/PRS_compare_EUR_incident.txt")





