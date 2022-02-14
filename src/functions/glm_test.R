library(data.table)                                                                                                                                                                                                                          
library(MASS)
library(pROC)

# Function to extract relevant details from logistic regression
# Calculating 95% ci can be quite slow with factors with many levels
# (e.g. assessmenet centre) so option to skip included for quick checks.
glm.test <- function(formula, event_col, data, skip.ci=FALSE) {
  # Detect if formula provided as string and convert
  if (length(formula) == 1 && !inherits(formula, "formula") && is.character(formula)) {
    formula <- as.formula(formula)
  }

  # Detect if response is a ordered factor with > 2 levels
  if (length(levels(data[[event_col]]) > 2) && is.ordered(data[[event_col]])) {
    dt <- ordinal.test(formula, event_col, data, skip.ci)
    setnames(dt, "T.value", "Z.score") # for consistency

    # simplify class counts into cases / controls
    dt[,.SD,.SDcols=names(dt)[!(names(dt) %like% event_col)]]
    dt[, controls := sum(data[[event_col]] == levels(data[[event_col]])[1])]
    dt[, cases := samples - controls]
   
    return(dt)
  }

  # Fit model and compute AUC
	g1 <- glm(formula, data=data, family="binomial")
	auc <- auc(data[[event_col]], predict(g1, data, type="response"))

  # Compute 95% confidence intervals
  if (!skip.ci) {
    ci <- confint(g1)
    ci.auc <- ci.auc(auc)
  }
 
  # harmonize coefficients
  cf <- coef(summary(g1))

  # Return coefficients
  if (skip.ci) {
    data.table(samples=data[,.N], cases=sum(data[[event_col]]), controls=data[,.N]-sum(data[[event_col]]),
               coefficient=rownames(cf), logOR=cf[,1], logOR.SE=cf[,2], OR=exp(cf[,1]), Z.score=cf[,3], P.value=cf[,4],
               AUC=c(auc))
  } else {
    data.table(samples=data[,.N], cases=sum(data[[event_col]]), controls=data[,.N]-sum(data[[event_col]]),
               coefficient=rownames(cf), logOR=cf[,1], logOR.SE=cf[,2], logOR.L95=ci[,1], logOR.U95=ci[,2],
               OR=exp(cf[,1]), OR.L95=exp(ci[,1]), OR.U95=exp(ci[,2]), Z.score=cf[,3], P.value=cf[,4], 
               AUC=c(auc), AUC.L95=c(ci.auc)[1], AUC.U95=c(ci.auc)[3])
  }
}

ordinal.test <- function(formula, event_col, data, skip.ci=FALSE) {
  # Detect if formula provided as string and convert
  if (length(formula) == 1 && !inherits(formula, "formula") && is.character(formula)) {
    formula <- as.formula(formula)
  }

  # Fit model and obtain AUC
	g1 <- polr(formula, data=data, Hess=TRUE)
  auc <- auc(multiclass.roc(data[[event_col]] ~ g1$lp))

  # Get Pvalue
  cf <- coef(summary(g1))
  pval <- pmin(1, 2*pt(abs(cf[,3]), g1$df.residual, lower.tail=FALSE))
  cf <- cbind(cf, pval)

  # Compute 95% confidence intervals
  if (!skip.ci) {
    ci <- confint(g1)
    # Add in missing values for 95% CI for intercept terms
    ci <- rbind(ci, matrix(NA_real_, ncol=ncol(ci), nrow=nrow(cf) - nrow(ci)))
  }
 
  # Get class numbers
  class_n <- dcast(data[,.N,by=hx_t2d], . ~ hx_t2d, value.var="N")
  class_n <- class_n[, .SD, .SDcols=levels(data[[event_col]])]
  setnames(class_n, sprintf("%s:%s", event_col, names(class_n)))
  info <- cbind(samples=data[,.N], class_n)

  # Return coefficients
  if (skip.ci) {
    cbind(info, coefficient=rownames(cf), logOR=cf[,1], logOR.SE=cf[,2], OR=exp(cf[,1]), T.value=cf[,3], 
          P.value=cf[,4], AUC=c(auc))
  } else {
    cbind(info, coefficient=rownames(cf), logOR=cf[,1], logOR.SE=cf[,2], logOR.L95=ci[,1], logOR.U95=ci[,2],
          OR=exp(cf[,1]), OR.L95=exp(ci[,1]), OR.U95=exp(ci[,2]), T.value=cf[,3], P.value=cf[,4], AUC=c(auc))
  }
}
