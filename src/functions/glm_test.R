library(data.table)                                                                                                                                                                                                                          
library(MASS)
library(pROC)

# Function to extract relevant details from logistic regression
# 
# Two methods for calculating the 95% confidence intervals of the
# odds ratios are provided: 'wald' and 'likelihood'. 
# 
# 'likelihood' is the standard approach used by 'confint()' when applied to
# a fitted logistic regression glm, but is quite slow.
#
# 'wald' is a much simpler approach that is very fast (10-20 times faster than
# the above), obtained by 'confint.default()', and is reasonably accurate given 
# certain assumptions are met. 
#
# For discussion, see:
# https://stats.stackexchange.com/questions/275416/computing-confidence-intervals-for-coefficients-in-logistic-regression
# https://thestatsgeek.com/2014/02/08/wald-vs-likelihood-ratio-test/
glm.test <- function(formula, event_col, data, ci.method="likelihood") {
  stopifnot(ci.method %in% c("wald", "likelihood"))

  # Detect if formula provided as string and convert
  if (length(formula) == 1 && !inherits(formula, "formula") && is.character(formula)) {
    formula <- as.formula(formula)
  }

  # Fit model and compute AUC
	g1 <- glm(formula, data=data, family="binomial")
  fitidx <- as.integer(names(g1$residuals))
	auc <- auc(data[[event_col]][fitidx], predict(g1, data[fitidx], type="response"))

  # Compute 95% confidence intervals
	if (ci.method == "wald") {
		ci <- confint.default(g1) 
	} else if (ci.method == "likelihood") {
		ci <- confint(g1)
	}
	ci.auc <- ci.auc(auc)
 
  # harmonize coefficients
  cf <- coef(summary(g1))

  # Return coefficients
	data.table(samples=data[fitidx,.N], cases=sum(data[[event_col]][fitidx]), 
             controls=data[fitidx,.N]-sum(data[[event_col]][fitidx]),
						 coefficient=rownames(cf), logOR=cf[,1], logOR.SE=cf[,2], logOR.L95=ci[,1], logOR.U95=ci[,2],
						 OR=exp(cf[,1]), OR.L95=exp(ci[,1]), OR.U95=exp(ci[,2]), Z.score=cf[,3], P.value=cf[,4], 
						 AUC=c(auc), AUC.L95=c(ci.auc)[1], AUC.U95=c(ci.auc)[3])
}

