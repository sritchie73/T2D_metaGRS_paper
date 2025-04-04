library(data.table)
library(ggplot2)

system("mkdir -p output/paper_outputs/{figures,tables}", wait=TRUE)

# Plot metaGRS training results
train <- fread("output/metaGRS/train/combined_basic_adjudicated/prefilter_none/enfit_model_fit.txt")

g <- ggplot(train) + 
  aes(x=lambda, y=fit_mean, ymin=fit_mean_minus_sd, ymax=fit_mean_plus_sd) +
  facet_grid(~ alpha, scales="free_x") +
  geom_ribbon(fill="#636363", alpha=0.6) +
  geom_line(linewidth=0.2) +
  geom_point(data=train[,.SD[which.max(fit_mean)]], shape=23, fill="white", color="#aa0000") +
  geom_errorbar(data=train[,.SD[which.max(fit_mean)]], color="#aa0000", width=0) +
  scale_x_continuous("Lambda penalty", trans="log10") +
  ylab("AUC (+/- SD)") +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=7),
    strip.background=element_blank(), strip.text=element_text(size=6, face="bold") 
  )

ggsave(g, width=7.2, height=1.8, file="output/paper_outputs/figures/metaGRS_training.pdf")

# Plot model coefficients
contrib <- fread("output/metaGRS/hyperparam_selection/T2D_metaGRS_component_GRS_contributions.txt")
contrib[, GRS := factor(GRS, levels=GRS)]

g <- ggplot(contrib) +
  aes(x=GRS, y=metaGRS_weight) +
  geom_hline(yintercept=0, linetype=2) +
  geom_point(shape=23, fill="white") +
  ylab("Elasticnet beta") + 
  xlab("") +
  theme_bw() +
  theme(
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    axis.text.x=element_text(size=7, angle=90, hjust=1, vjust=0.5, color="black"),
    axis.text.y=element_text(size=6),
    axis.title=element_text(size=8),
  )
ggsave(g, width=7.2, height=2.7, units="in", file="output/paper_outputs/figures/T2D_metaGRS_component_GRS_contributions.pdf")

# Write out per-PRS contributions to the metaPRS
fwrite(contrib, sep="\t", quote=FALSE, file="output/paper_outputs/tables/T2D_metaGRS_component_GRS_contributions.txt")




