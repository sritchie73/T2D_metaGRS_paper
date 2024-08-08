library(data.table)
library(ggplot2)
library(ggstance)
library(scales)

system("mkdir -p output/paper_outputs/{figures,tables}", wait=TRUE)

# Create categorical NRI plot
nri_estimates <- fread("output/UKB_tests/categorical_nri_estimates.txt")
nri_estimates[, risk_threshold := factor(risk_threshold, levels=c("15%", "10%", "5%"))]
nri_estimates <- nri_estimates[metric %in% c("NRI+", "NRI-")]
nri_estimates[, metric := ifelse(metric == "NRI+", "T2D cases", "Non-cases")]
nri_estimates[, metric := factor(metric, levels=c("Non-cases", "T2D cases"))]
nri_estimates[, QDiabetes_model := gsub(" model", "", QDiabetes_model)]

g <- ggplot(nri_estimates) +
  aes(x=Estimate, xmin=L95, xmax=U95, y=risk_threshold, color=metric) +
  facet_wrap(~ QDiabetes_model) +
  geom_errorbarh(height=0, position=position_dodgev(heigh=0.6)) +
  geom_point(shape=23, fill="white", position=position_dodgev(height=0.6)) +
  scale_color_manual(values=c("T2D cases"="#c51b7d", "Non-cases"="#4d9221")) +
  geom_vline(xintercept=0, linetype=2) +
  scale_x_continuous("% reclassified (95% CI)", labels=scales::percent, limits=c(-0.05, 0.15)) +
  scale_y_discrete("Risk threshold") +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=7), strip.background=element_blank(),
    strip.text=element_text(size=7, face="bold"), legend.position="bottom",
    legend.text=element_text(size=6), legend.title=element_blank()
  )
ggsave(g, width=3.3, height=1.6, file="output/paper_outputs/figures/categorical_nri.pdf")

# Curate supp tables
supp <- copy(nri_estimates)
supp[, QDiabetes_model := paste(gsub(" ", " 2018 model ", QDiabetes_model), "+ T2D metaPRS")]
supp <- supp[,.(risk_threshold, model=QDiabetes_model, group=metric, pct_reclassified=Estimate, L95, U95, Pval)]

# More complicated than collating the NRI %, we now curate the total numbers reclassified
strata_num <- fread("output/UKB_tests/categorical_nri_reclassified.txt")
strata_num[, QDiabetes_model := gsub("model", "2018 model", QDiabetes_model)]

supp_strata <- strata_num[, .SD[1], by=.(risk_threshold, QDiabetes_model)]
supp_strata <- supp_strata[,.(risk_threshold, model=QDiabetes_model, Samples=Total_Samples, Cases=Total_Cases)] 

# N.b. "Samples", and "Cases" are placeholder column names, denoting different things depending on the group/strategy:
# - For Non-cases, the "Cases" column holds the total number of Non-cases
supp_strata <- rbind(
  supp_strata[,.(risk_threshold, model, group="T2D cases", Samples, Cases)],
  supp_strata[,.(risk_threshold, model=paste(model, "+ T2D metaPRS"), group="T2D cases", Samples, Cases)],
  supp_strata[,.(risk_threshold, model, group="Non-cases", Samples, Cases=Samples-Cases)],
  supp_strata[,.(risk_threshold, model=paste(model, "+ T2D metaPRS"), group="Non-cases", Samples, Cases=Samples-Cases)]
)

# Curate the number of cases allocated to low and high risk groups, as well as the number reclassified when adding the metaPRS
# - N.b. for controls, the risk groups are swapped in meaning, and cases are controls
low_risk <- rbind(
  strata_num[QDiabetes %like% "<" , .(Cases=sum(Cases), group="T2D cases"), by=.(risk_threshold, model=QDiabetes_model)],
  strata_num[QDiabetes_plus_PRS %like% "<", .(Cases=sum(Cases), group="T2D cases"), by=.(risk_threshold, model=paste(QDiabetes_model, "+ T2D metaPRS"))],
  strata_num[QDiabetes %like% "<" , .(Cases=sum(All - Cases), group="Non-cases"), by=.(risk_threshold, model=QDiabetes_model)],
  strata_num[QDiabetes_plus_PRS %like% "<", .(Cases=sum(All - Cases), group="Non-cases"), by=.(risk_threshold, model=paste(QDiabetes_model, "+ T2D metaPRS"))]
)
supp_strata[low_risk, on = .(risk_threshold, model, group), Low.Cases := i.Cases]

reclassified_low_risk <- rbind(
  strata_num[QDiabetes %like% ">=" & QDiabetes_plus_PRS %like% "<", .(risk_threshold, model=paste(QDiabetes_model, "+ T2D metaPRS"), group="T2D cases", Cases)],
  strata_num[QDiabetes %like% ">=" & QDiabetes_plus_PRS %like% "<", .(risk_threshold, model=paste(QDiabetes_model, "+ T2D metaPRS"), group="Non-cases", Cases=All - Cases)]
)
supp_strata[reclassified_low_risk, on = .(risk_threshold, model, group), Reclassified.Low.Cases := i.Cases]

high_risk <- rbind(
  strata_num[QDiabetes %like% ">=" , .(Cases=sum(Cases), group="T2D cases"), by=.(risk_threshold, model=QDiabetes_model)],
  strata_num[QDiabetes_plus_PRS %like% ">=", .(Cases=sum(Cases), group="T2D cases"), by=.(risk_threshold, model=paste(QDiabetes_model, "+ T2D metaPRS"))],
  strata_num[QDiabetes %like% ">=" , .(Cases=sum(All - Cases), group="Non-cases"), by=.(risk_threshold, model=QDiabetes_model)],
  strata_num[QDiabetes_plus_PRS %like% ">=", .(Cases=sum(All - Cases), group="Non-cases"), by=.(risk_threshold, model=paste(QDiabetes_model, "+ T2D metaPRS"))]
)
supp_strata[high_risk, on = .(risk_threshold, model, group), High.Cases := i.Cases]

reclassified_high_risk <- rbind(
  strata_num[QDiabetes %like% "<" & QDiabetes_plus_PRS %like% ">=", .(risk_threshold, model=paste(QDiabetes_model, "+ T2D metaPRS"), group="T2D cases", Cases)],
  strata_num[QDiabetes %like% "<" & QDiabetes_plus_PRS %like% ">=", .(risk_threshold, model=paste(QDiabetes_model, "+ T2D metaPRS"), group="Non-cases", Cases=All - Cases)]
)
supp_strata[reclassified_high_risk, on = .(risk_threshold, model, group), Reclassified.High.Cases := i.Cases]

# Combine with NRI estimates
supp <- supp[supp_strata, on = .( risk_threshold, model, group)]

# Reorder rows
supp <- supp[order(model)][order(as.numeric(gsub("%.*", "", risk_threshold)))][order(group)]

# Split out cases and non-cases
supp_cases <- supp[group == "T2D cases"]
supp_cases[, group := NULL]

supp_controls <- supp[group == "Non-cases"]
supp_controls[, group := NULL]
setnames(supp_controls, gsub("Cases", "Controls", names(supp_controls)))

# Write out
fwrite(supp_cases, sep="\t", quote=FALSE, file="output/paper_outputs/tables/categorical_nri_cases.txt")
fwrite(supp_controls, sep="\t", quote=FALSE, file="output/paper_outputs/tables/categorical_nri_controls.txt")

