library(data.table)
library(openxlsx)
library(ggplot2)
library(gghalves)
library(ggbeeswarm)
library(patchwork)

system("mkdir -p output/paper_outputs/{figures,tables}", wait=TRUE)

# Unified theme for plotting
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

# Load prevalent UK Biobank analyses in non-Europeans
prev_ukb <- fread("output/UKB_tests/non_european_prevalent_T2D_associations.txt")
prev_ukb <- prev_ukb[model == coefficient & !(model %like% "^Gad")]
prev_ukb[model == "Aly2021_PGS000864", model := "MonsourAly2021_PGS000864"]

# Extract three major ancestry groupings
prev_ukb_SAS <- prev_ukb[ancestry == "SAS"]
prev_ukb_AFR <- prev_ukb[ancestry == "AFR"]
prev_ukb_EAS <- prev_ukb[ancestry == "EAS"]

# Order by PRS performance
prev_ukb_SAS <- prev_ukb_SAS[order(-OR)]
prev_ukb_AFR <- prev_ukb_AFR[order(-OR)]
prev_ukb_EAS <- prev_ukb_EAS[order(-OR)]

# Create factor level ordering for PRS for plotting
prev_ukb_SAS[, PRS := factor(model, levels=rev(unique(model)))]
prev_ukb_AFR[, PRS := factor(model, levels=rev(unique(model)))]
prev_ukb_EAS[, PRS := factor(model, levels=rev(unique(model)))]

# Add annotation column for coloring T2D metaPRS
prev_ukb_SAS[, color_anno := ifelse(PRS == "T2D_metaGRS", TRUE, FALSE)]
prev_ukb_AFR[, color_anno := ifelse(PRS == "T2D_metaGRS", TRUE, FALSE)]
prev_ukb_EAS[, color_anno := ifelse(PRS == "T2D_metaGRS", TRUE, FALSE)]

# Create plots
gg_prev_ukb_SAS <- ggplot(prev_ukb_SAS[1:6]) +
  aes(x=OR, xmin=OR.L95, xmax=OR.U95, y=PRS, color=color_anno) +
  geom_vline(xintercept=1, linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, fill="white") +
  scale_color_manual(values=c("TRUE"="#bd0026", "FALSE"="#08306b")) +
  scale_x_continuous("Odds Ratio (95% CI)", limits=c(1, 2.5), breaks=c(1,1.5,2,2.5)) +
  ggtitle("Prediction of T2D prevalence in UK Biobank", subtitle="N=6,992; 1,253 prevalent T2D cases") +
  theme_prs()

gg_prev_ukb_AFR <- ggplot(prev_ukb_AFR[1:6]) +
  aes(x=OR, xmin=OR.L95, xmax=OR.U95, y=PRS, color=color_anno) +
  geom_vline(xintercept=1, linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, fill="white") +
  scale_color_manual(values=c("TRUE"="#bd0026", "FALSE"="#08306b")) +
  scale_x_continuous("Odds Ratio (95% CI)", limits=c(1, 2.5), breaks=c(1,1.5,2,2.5)) +
  ggtitle("Prediction of T2D prevalence in UK Biobank", subtitle="N=6,871; 766 prevalent T2D cases") +
  theme_prs()

gg_prev_ukb_EAS <- ggplot(prev_ukb_EAS[1:6]) +
  aes(x=OR, xmin=OR.L95, xmax=OR.U95, y=PRS, color=color_anno) +
  geom_vline(xintercept=1, linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, fill="white") +
  scale_color_manual(values=c("TRUE"="#bd0026", "FALSE"="#08306b")) +
  scale_x_continuous("Odds Ratio (95% CI)", limits=c(1, 2.5), breaks=c(1,1.5,2,2.5)) +
  ggtitle("Prediction of T2D prevalence in UK Biobank", subtitle="N=1,432; 76 prevalent T2D cases") +
  theme_prs()

# Load incident UK Biobank analyses in non-Europeans
inci_ukb <- fread("output/UKB_tests/non_european_incident_T2D_associations.txt")
inci_ukb <- inci_ukb[model == coefficient & !(model %like% "^Gad")]
inci_ukb[model == "Aly2021_PGS000864", model := "MonsourAly2021_PGS000864"]

# Extract three major ancestry groupings
inci_ukb_SAS <- inci_ukb[ancestry == "SAS"]
inci_ukb_AFR <- inci_ukb[ancestry == "AFR"]
inci_ukb_EAS <- inci_ukb[ancestry == "EAS"]

# Order by PRS performance
inci_ukb_SAS <- inci_ukb_SAS[order(-HR)]
inci_ukb_AFR <- inci_ukb_AFR[order(-HR)]
inci_ukb_EAS <- inci_ukb_EAS[order(-HR)]

# Create factor level ordering for PRS for plotting
inci_ukb_SAS[, PRS := factor(model, levels=rev(unique(model)))]
inci_ukb_AFR[, PRS := factor(model, levels=rev(unique(model)))]
inci_ukb_EAS[, PRS := factor(model, levels=rev(unique(model)))]

# Add annotation column for coloring T2D metaPRS
inci_ukb_SAS[, color_anno := ifelse(PRS == "T2D_metaGRS", TRUE, FALSE)]
inci_ukb_AFR[, color_anno := ifelse(PRS == "T2D_metaGRS", TRUE, FALSE)]
inci_ukb_EAS[, color_anno := ifelse(PRS == "T2D_metaGRS", TRUE, FALSE)]

# Create plots
gg_inci_ukb_SAS <- ggplot(inci_ukb_SAS[1:6]) +
  aes(x=HR, xmin=L95, xmax=U95, y=PRS, color=color_anno) +
  geom_vline(xintercept=1, linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, fill="white") +
  scale_color_manual(values=c("TRUE"="#bd0026", "FALSE"="#08306b")) +
  scale_x_continuous("Hazard Ratio (95% CI)", limits=c(1, 2.5), breaks=c(1,1.5,2,2.5)) +
  ggtitle("Prediction of incident T2D risk in UK Biobank", subtitle="N=5,685; 513 T2D cases within 10 years of follow-up") +
  theme_prs()

gg_inci_ukb_AFR <- ggplot(inci_ukb_AFR[1:6]) +
  aes(x=HR, xmin=L95, xmax=U95, y=PRS, color=color_anno) +
  geom_vline(xintercept=1, linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, fill="white") +
  scale_color_manual(values=c("TRUE"="#bd0026", "FALSE"="#08306b")) +
  scale_x_continuous("Hazard Ratio (95% CI)", limits=c(1, 2.5), breaks=c(1,1.5,2,2.5)) +
  ggtitle("Prediction of incident T2D risk in UK Biobank", subtitle="N=6,019; 395 T2D cases within 10 years of follow-up") +
  theme_prs()

gg_inci_ukb_EAS <- ggplot(inci_ukb_EAS[1:6]) +
  aes(x=HR, xmin=L95, xmax=U95, y=PRS, color=color_anno) +
  geom_vline(xintercept=1, linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, fill="white") +
  scale_color_manual(values=c("TRUE"="#bd0026", "FALSE"="#08306b")) +
  scale_x_continuous("Hazard Ratio (95% CI)", limits=c(1, 2.5), breaks=c(1,1.5,2,2.5)) +
  ggtitle("Prediction of incident T2D risk in UK Biobank", subtitle="N=1,350; 37 T2D cases within 10 years of follow-up") +
  theme_prs()

# Load All of Us analyses
prev_aou <- fread("data/Henry_Taylor_AoU_cohort/glm_results__pgs_covariates.tsv")

# Rename PRSs for consistency
prs_map <- prev_aou[,.(AoU_name=grs)]
prs_map[, PGS_ID := gsub(".*_", "", AoU_name)]
prs_map[PGS_ID == "chagoya", PGS_ID := "PGS00344X"]
prs_map2 <- inci_ukb_EAS[,.(name=PRS)]
prs_map2[, PGS_ID := gsub(".*_", "", name)]
prs_map <- merge(prs_map, prs_map2, by="PGS_ID", all=TRUE)
prev_aou[prs_map, on = .(grs=AoU_name), PRS := i.name]
prev_aou[grs == "pgs_PGS003867", PRS := "Shim2023_PGS003867"]
prev_aou[PRS == "Aly2021_PGS000864", PRS := "MonsourAly2021_PGS000864"]

# Extract major non-European ancestry groupings
prev_aou_AMR <- prev_aou[ancestry == "amr"]
prev_aou_AFR <- prev_aou[ancestry == "afr"]

# Order by PRS performance
prev_aou_AMR <- prev_aou_AMR[order(-OR)]
prev_aou_AFR <- prev_aou_AFR[order(-OR)]

# Create factor level ordering for PRS for plotting
prev_aou_AMR[, PRS := factor(PRS, levels=rev(unique(PRS)))]
prev_aou_AFR[, PRS := factor(PRS, levels=rev(unique(PRS)))]

# Add annotation column for coloring T2D metaPRS
prev_aou_AMR[, color_anno := ifelse(PRS == "T2D_metaGRS", TRUE, FALSE)]
prev_aou_AFR[, color_anno := ifelse(PRS == "T2D_metaGRS", TRUE, FALSE)]

# Create plots
gg_prev_aou_AMR <- ggplot(prev_aou_AMR[1:6]) +
  aes(x=OR, xmin=OR.L95, xmax=OR.U95, y=PRS, color=color_anno) +
  geom_vline(xintercept=1, linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, fill="white") +
  scale_color_manual(values=c("TRUE"="#bd0026", "FALSE"="#08306b")) +
  scale_x_continuous("Odds Ratio (95% CI)", limits=c(1, 2.5), breaks=c(1,1.5,2,2.5)) +
  ggtitle("Prediction of T2D prevalence in All of Us", subtitle="N=33,652; 4,033 prevalent T2D cases") +
  theme_prs()

gg_prev_aou_AFR <- ggplot(prev_aou_AFR[1:6]) +
  aes(x=OR, xmin=OR.L95, xmax=OR.U95, y=PRS, color=color_anno) +
  geom_vline(xintercept=1, linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, fill="white") +
  scale_color_manual(values=c("TRUE"="#bd0026", "FALSE"="#08306b")) +
  scale_x_continuous("Odds Ratio (95% CI)", limits=c(1, 2.5), breaks=c(1,1.5,2,2.5)) +
  ggtitle("Prediction of T2D prevalence in All of Us", subtitle="N=44,364; 5,663 prevalent T2D cases") +
  theme_prs()

# Load in MEC analyses
inci_mec_EAS <- read.xlsx("data/Yujian_Liang_MEC_cohort/metaPRS.MEC_validation.logistic_res.CHS.xlsx")
inci_mec_SAS <- read.xlsx("data/Yujian_Liang_MEC_cohort/metaPRS.MEC_validation.logistic_res.INS.xlsx")
inci_mec_ASN <- read.xlsx("data/Yujian_Liang_MEC_cohort/metaPRS.MEC_validation.logistic_res.MAS.xlsx")
setDT(inci_mec_EAS)
setDT(inci_mec_SAS)
setDT(inci_mec_ASN)
inci_mec <- rbind(inci_mec_EAS, inci_mec_SAS, inci_mec_ASN)

# Rename PRSs for consistency
inci_mec[, PRS := Study_ID]
inci_mec[Study_ID == "T2D_metaGRS_v3_9e639289_hmPOSGRCh38", PRS := "T2D_metaGRS"]
inci_mec[Study_ID == "HuertaChagoya2023", PRS := "HuertaChagoya2023_PGS00344X"]
inci_mec[PRS == "Aly2021_PGS000864", PRS := "MonsourAly2021_PGS000864"]

# Correct sign of Ye2021 PRS
inci_mec[PRS == "Ye2021_PGS001357", 
  c("logOR", "logOR.L95", "logOR.U95", "OR", "OR.L95", "OR.U95", "Z.score") := 
  .(-logOR, -logOR.L95, -logOR.U95, exp(-logOR), exp(-logOR.U95), exp(-logOR.L95), -Z.score)
]

# Define and extract ancestry groups
inci_mec_SAS <- inci_mec[Genetic.ancestry == "Indian"]
inci_mec_EAS <- inci_mec[Genetic.ancestry == "Chinese"]
inci_mec_ASN <- inci_mec[Genetic.ancestry == "Malay"]

# Order by PRS performance
inci_mec_SAS <- inci_mec_SAS[order(-OR)]
inci_mec_EAS <- inci_mec_EAS[order(-OR)]
inci_mec_ASN <- inci_mec_ASN[order(-OR)]

# Create factor level ordering for PRS for plotting
inci_mec_SAS[, PRS := factor(PRS, levels=rev(unique(PRS)))]
inci_mec_EAS[, PRS := factor(PRS, levels=rev(unique(PRS)))]
inci_mec_ASN[, PRS := factor(PRS, levels=rev(unique(PRS)))]

# Add annotation column for coloring T2D metaPRS
inci_mec_SAS[, color_anno := ifelse(PRS == "T2D_metaGRS", TRUE, FALSE)]
inci_mec_EAS[, color_anno := ifelse(PRS == "T2D_metaGRS", TRUE, FALSE)]
inci_mec_ASN[, color_anno := ifelse(PRS == "T2D_metaGRS", TRUE, FALSE)]

# Create plots
gg_inci_mec_SAS <- ggplot(inci_mec_SAS[1:6]) +
  aes(x=OR, xmin=OR.L95, xmax=OR.U95, y=PRS, color=color_anno) +
  geom_vline(xintercept=1, linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, fill="white") +
  scale_color_manual(values=c("TRUE"="#bd0026", "FALSE"="#08306b")) +
  scale_x_continuous("Odds Ratio (95% CI)", limits=c(1, 2.5), breaks=c(1,1.5,2,2.5)) +
  ggtitle("Prediction of incident T2D risk in MEC", subtitle="N=852; 194 T2D cases within 7 years of follow-up") +
  theme_prs()

gg_inci_mec_EAS <- ggplot(inci_mec_EAS[1:6]) +
  aes(x=OR, xmin=OR.L95, xmax=OR.U95, y=PRS, color=color_anno) +
  geom_vline(xintercept=1, linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, fill="white") +
  scale_color_manual(values=c("TRUE"="#bd0026", "FALSE"="#08306b")) +
  scale_x_continuous("Odds Ratio (95% CI)", limits=c(1, 2.5), breaks=c(1,1.5,2,2.5)) +
  ggtitle("Prediction of incident T2D risk in MEC", subtitle="N=1,149; 205 T2D cases within 7 years of follow-up") +
  theme_prs()

gg_inci_mec_ASN <- ggplot(inci_mec_ASN[1:6]) +
  aes(x=OR, xmin=OR.L95, xmax=OR.U95, y=PRS, color=color_anno) +
  geom_vline(xintercept=1, linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, fill="white") +
  scale_color_manual(values=c("TRUE"="#bd0026", "FALSE"="#08306b")) +
  scale_x_continuous("Odds Ratio (95% CI)", limits=c(1, 2.5), breaks=c(1,1.5,2,2.5)) +
  ggtitle("Prediction of incident T2D risk in MEC", subtitle="N=870; 187 T2D cases within 7 years of follow-up") +
  theme_prs()


# Create combined plot
g <- (gg_prev_ukb_SAS / gg_inci_ukb_SAS / gg_inci_mec_SAS / gg_inci_mec_ASN) | 
     (gg_prev_ukb_EAS / gg_inci_ukb_EAS / gg_inci_mec_EAS / gg_prev_aou_AMR) |
     (gg_prev_ukb_AFR / gg_inci_ukb_AFR / gg_prev_aou_AFR / gg_prev_aou_AFR) 

ggsave(g, width=7.2, height=5.282, file="output/paper_outputs/figures/PRS_compare_non_EUR.pdf")

# Load in EUR results
prev_eur <- fread("output/paper_outputs/tables/PRS_compare_EUR_prevalent.txt")
inci_eur <- fread("output/paper_outputs/tables/PRS_compare_EUR_incident.txt")

# Create plot comparing the overall ranking of PRSs across multiple ancestries
ranks <- rbind(
  prev_ukb_SAS[,.(cohort="UK Biobank (prevalent T2D)", ancestry="SAS", type="prevalent", PRS, beta=log(OR))],
  prev_ukb_EAS[,.(cohort="UK Biobank (prevalent T2D)", ancestry="EAS", type="prevalent", PRS, beta=log(OR))],
  prev_ukb_AFR[,.(cohort="UK Biobank (prevalent T2D)", ancestry="AFR", type="prevalent", PRS, beta=log(OR))],
  inci_ukb_SAS[,.(cohort="UK Biobank (incident T2D)", ancestry="SAS", type="incident", PRS, beta=log(HR))],
  inci_ukb_EAS[,.(cohort="UK Biobank (incident T2D)", ancestry="EAS", type="incident", PRS, beta=log(HR))],
  inci_ukb_AFR[,.(cohort="UK Biobank (incident T2D)", ancestry="AFR", type="incident", PRS, beta=log(HR))],
  prev_aou_AMR[,.(cohort="All of Us", ancestry="AMR", type="prevalent", PRS, beta=log(OR))],
  prev_aou_AFR[,.(cohort="All of Us", ancestry="AFR", type="prevalent", PRS, beta=log(OR))],
  inci_mec_SAS[,.(cohort="Singapore Multi-Ethnic Cohort", ancestry="SAS", type="incident", PRS, beta=log(OR))],
  inci_mec_EAS[,.(cohort="Singapore Multi-Ethnic Cohort", ancestry="EAS", type="incident", PRS, beta=log(OR))],
  inci_mec_ASN[,.(cohort="Singapore Multi-Ethnic Cohort", ancestry="ASN", type="incident", PRS, beta=log(OR))],
  prev_eur[cohort == "UK Biobank", .(cohort="UK Biobank (prevalent T2D)", ancestry="EUR", type="prevalent", PRS, beta=log(OR))],
  prev_eur[cohort == "All of Us", .(cohort, ancestry="EUR", type="prevalent", PRS, beta=log(OR))],
  inci_eur[cohort == "UK Biobank", .(cohort="UK Biobank (incident T2D)", ancestry="EUR", type="incident", PRS, beta=log(HR))],
  inci_eur[cohort == "INTERVAL", .(cohort, ancestry="EUR", type="incident", PRS, beta=log(HR))]
)
ranks[PRS == "Aly2021_PGS000864", PRS := "MonsourAly2021_PGS000864"]

ranks[ranks[,.SD[which.max(beta)], by=PRS], on = .(PRS), portability := beta / i.beta] # How portable is each PRS with respect to its own max OR/HR?
ranks[, rel_EUR_metaPRS := beta / ranks[which.max(beta), beta]] # How strong / portable is each PRS with respect to the metaPRS in EUR UKB?
ranks[ranks[,.SD[which.max(beta)],by=.(cohort, ancestry, type)], on = .(cohort, ancestry, type), rel_cohort := beta / i.beta] # standardise to maximum log OR or log HR per panel

# Plot relative ranking of PRS within each cohort and ancestry
ggdt <- copy(ranks)
prs_order <- ggdt[,.(median=median(portability)), by=PRS][order(-median)] # Determine plot order for PRS
ggdt[, PRS := factor(PRS, levels=prs_order$PRS)]

g <- ggplot(ggdt) + 
  aes(x=PRS, y=portability) + 
  geom_beeswarm(aes(fill=ancestry, shape=cohort), side=1, size=1, color="black", stroke=0.5, cex=0.8) +
  geom_half_violin(side="l", color="black", scale="count", width=1) +
  geom_hline(yintercept=0, linetype=2) +
  scale_fill_manual("Ancestry", values=c("AFR"="#4daf4a", "AMR"="#377eb8", "EAS"="#ff7f00", "SAS"="#984ea3", "ASN"="#ffff33", "EUR"="#e41a1c")) +
  scale_shape_manual("Cohort", values=c("UK Biobank (prevalent T2D)"=21, "UK Biobank (incident T2D)"=24, "All of Us"=22, "Singapore Multi-Ethnic Cohort"=23, "INTERVAL"=25)) +
  guides(fill = guide_legend(override.aes = list(shape = 21), nrow=1)) +
  ylab("Portability") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
    axis.text.x=element_text(size=6, color="black", angle=45, hjust=1, vjust=1), 
    axis.title.x=element_blank(),
    panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
    legend.position="bottom", legend.text=element_text(size=7), legend.box="vertical",
    legend.margin=margin(0,0,0,0), legend.title=element_text(size=8)
  )
ggsave(g, width=7.2, height=4, file="output/paper_outputs/figures/PRS_portability_within_PRS.pdf")

# Plot relative ranking of PRS within each cohort and ancestry
ggdt <- copy(ranks)
prs_order <- ggdt[,.(median=median(rel_EUR_metaPRS)), by=PRS][order(-median)] # Determine plot order for PRS
ggdt[, PRS := factor(PRS, levels=prs_order$PRS)]

g <- ggplot(ggdt) + 
  aes(x=PRS, y=rel_EUR_metaPRS) + 
  geom_beeswarm(aes(fill=ancestry, shape=cohort), side=1, size=1, color="black", stroke=0.5, cex=0.8) +
  geom_half_violin(side="l", color="black", scale="count", width=1) +
  geom_hline(yintercept=0, linetype=2) +
  scale_fill_manual("Ancestry", values=c("AFR"="#4daf4a", "AMR"="#377eb8", "EAS"="#ff7f00", "SAS"="#984ea3", "ASN"="#ffff33", "EUR"="#e41a1c")) +
  scale_shape_manual("Cohort", values=c("UK Biobank (prevalent T2D)"=21, "UK Biobank (incident T2D)"=24, "All of Us"=22, "Singapore Multi-Ethnic Cohort"=23, "INTERVAL"=25)) +
  guides(fill = guide_legend(override.aes = list(shape = 21), nrow=1)) +
  ylab("Relative to metaPRS in EUR UKB") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
    axis.text.x=element_text(size=6, color="black", angle=45, hjust=1, vjust=1), 
    axis.title.x=element_blank(),
    panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
    legend.position="bottom", legend.text=element_text(size=7), legend.box="vertical",
    legend.margin=margin(0,0,0,0), legend.title=element_text(size=8)
  )
ggsave(g, width=7.2, height=4, file="output/paper_outputs/figures/PRS_portability_relative_to_EUR_UKB_metaPRS.pdf")

# Plot relative ranking of PRS within each cohort and ancestry
ggdt <- copy(ranks)
prs_order <- ggdt[,.(median=median(rel_cohort)), by=PRS][order(-median)] # Determine plot order for PRS
ggdt[, PRS := factor(PRS, levels=prs_order$PRS)]

g <- ggplot(ggdt) + 
  aes(x=PRS, y=rel_cohort) + 
  geom_beeswarm(aes(fill=ancestry, shape=cohort), side=1, size=1, color="black", stroke=0.5, cex=0.8) +
  geom_half_violin(side="l", color="black", scale="count", width=1) +
  geom_hline(yintercept=0, linetype=2) +
  scale_fill_manual("Ancestry", values=c("AFR"="#4daf4a", "AMR"="#377eb8", "EAS"="#ff7f00", "SAS"="#984ea3", "ASN"="#ffff33", "EUR"="#e41a1c")) +
  scale_shape_manual("Cohort", values=c("UK Biobank (prevalent T2D)"=21, "UK Biobank (incident T2D)"=24, "All of Us"=22, "Singapore Multi-Ethnic Cohort"=23, "INTERVAL"=25)) +
  guides(fill = guide_legend(override.aes = list(shape = 21), nrow=1)) +
  ylab("Relative prediction") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
    axis.text.x=element_text(size=6, color="black", angle=45, hjust=1, vjust=1), 
    axis.title.x=element_blank(),
    panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
    legend.position="bottom", legend.text=element_text(size=7), legend.box="vertical",
    legend.margin=margin(0,0,0,0), legend.title=element_text(size=8)
  )
ggsave(g, width=7.2, height=4, file="output/paper_outputs/figures/PRS_transferrability_within_cohort_ranking.pdf")

# create tables
inci_mec[, ancestry := fcase(
  Genetic.ancestry == "Indian", "SAS",
  Genetic.ancestry == "Chinese", "EAS",
  Genetic.ancestry == "Malay", "ASN"
)]
inci_mec[, samples := fcase(
  Genetic.ancestry == "Indian", 852,
  Genetic.ancestry == "Chinese", 1149,
  Genetic.ancestry == "Malay", 870
)]

logistic_dt <- rbind(idcol="cohort",
  "UK Biobank"=prev_ukb[,.(ancestry, samples, cases, PRS=model, OR, OR.L95, OR.U95, P.value, AUC, AUC.L95, AUC.U95)],
  "All of Us"=prev_aou[ancestry != "eur", .(ancestry=toupper(ancestry), samples, cases, PRS, OR, OR.L95, OR.U95, P.value, AUC, AUC.L95, AUC.U95)],
  "Singapore MEC"=inci_mec[,.(ancestry, samples, cases=Cases.Number, PRS, OR, OR.L95, OR.U95, P.value, AUC, AUC.L95, AUC.U95)]
)
logistic_dt[, ancestry := paste0("1KG-", ancestry, "-like")]
logistic_dt[ancestry == "1KG-ASN-like", ancestry := "ASN-like"]
logistic_dt <- logistic_dt[order(-abs(OR))][order(P.value)][order(ancestry)][order(cohort)]

cox_dt <- rbind(idcol="cohort",
  "UK Biobank"=inci_ukb[,.(ancestry, samples=Samples, cases=Cases, PRS=model, HR, HR.L95=L95, HR.U95=U95, P.value=Pvalue, C.index, C.L95, C.U95)]
)
cox_dt[, ancestry := paste0("1KG-", ancestry, "-like")]
cox_dt <- cox_dt[order(-abs(HR))][order(P.value)][order(ancestry)][order(cohort)]

fwrite(logistic_dt, sep="\t", quote=FALSE, file="output/paper_outputs/tables/PRS_compare_non_EUR_logistic.txt")
fwrite(cox_dt, sep="\t", quote=FALSE, file="output/paper_outputs/tables/PRS_compare_non_EUR_coxph.txt")

