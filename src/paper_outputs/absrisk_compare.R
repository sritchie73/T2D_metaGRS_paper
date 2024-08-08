library(data.table)
library(ggplot2)
library(scales)

system("mkdir -p output/paper_outputs/{figures,tables}", wait=TRUE)

# Load QDiabetes + metaPRS absolute risks
absrisk <- fread("output/UKB_tests/QDiabetes_plus_metaPRS_absrisks.txt")

# Melt to long format for plotting eCDFs
absrisk <- melt(absrisk, id.vars=c("eid", "sex", "age", "age_tertile", "incident_type_2_diabetes", "incident_censor_years"), value.name="predicted", variable.name="model")
absrisk <- absrisk[model != "T2D_metaGRS"] # drop column - too lazy to write out full measure.vars
absrisk <- absrisk[!is.na(predicted)]
absrisk[, QDiabetes_model := fcase(
  model %in% c("QDiabetes2018A", "QDiabA_plus_PRS"), "QDiabetes model A",
  model %in% c("QDiabetes2018B_non_fasting", "QDiabB_plus_PRS"), "QDiabetes model B",
  model %in% c("QDiabetes2018C", "QDiabC_plus_PRS"), "QDiabetes model C"
)]
absrisk[, model := ifelse(model %like% "plus_PRS", "QDiabetes + metaPRS", "QDiabetes")]
absrisk[, group := ifelse(incident_type_2_diabetes, "Case", "Control")]

# Plot probability of exceeding X% risk, stratified by case/control status
absrisk[, eCDF := ecdf(predicted)(predicted), by=.(QDiabetes_model, model, incident_type_2_diabetes)]
g <- ggplot(absrisk) +
  aes(x=predicted, y=1-eCDF, color=model) +
  facet_grid(group ~ QDiabetes_model) +
  geom_line(linewidth=0.6) +
  scale_color_manual(values=c("QDiabetes"="#08306b", "QDiabetes + metaPRS"="#bd0026")) +
  scale_x_continuous("Predicted 10-year T2D risk", labels=scales::percent, limits=c(0,0.3), oob=oob_keep, expand=expansion(mult=0, add=c(0.01, 0.03))) +
  scale_y_continuous("Probability of 10-year T2D risk exceeding X%", labels=scales::percent) +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=7), strip.background=element_blank(),
    strip.text=element_text(size=7, face="bold"), legend.position="bottom",
    legend.text=element_text(size=7), legend.title=element_blank()
  )
ggsave(g, width=3.4, height=2.5, file="output/paper_outputs/figures/QDiabetes_plus_PRS_absrisks.pdf")

