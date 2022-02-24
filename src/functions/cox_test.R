library(survival)
library(data.table)

# Function to extract all relevant details from cox models
cox.test <- function(formula, event_col, data) {
  if (length(formula) == 1 && !inherits(formula, "formula") && is.character(formula)) {
    formula <- as.formula(formula)
  }

  cx <- coxph(formula, data=data, x=TRUE) # fit the cox model
  ci <- confint(cx) # get 95% confidence intervals
  cz <- cox.zph(cx, terms=FALSE) # test proportional hazards assumption
  cc <- as.integer(rownames(cx$x)) # get complete cases used to fit the model

  # now get C-index and its 95% CI
  cindex <- summary(cx)$concordance[1]
  cindex.se <- summary(cx)$concordance[2]
  cindex.l95 <- cindex - qnorm(0.975)*cindex.se
  cindex.u95 <- cindex + qnorm(0.975)*cindex.se

  # collate information about all model coefficients
  cx <- as.data.table(coef(summary(cx)), keep.rownames="coefficient")
  ci <- as.data.table(ci, keep.rownames="coefficient")
  cz <- as.data.table(cz$table, keep.rownames="coefficient")
  dt <- cx[ci, on=.(coefficient)][cz, on=.(coefficient), nomatch=0]

  dt <- dt[, .(coefficient, logHR=coef, SE=`se(coef)`, HR=`exp(coef)`,
              L95=exp(`2.5 %`), U95=exp(`97.5 %`), Pvalue=`Pr(>|z|)`,
              Proportionality.chisq=chisq, Proportionality.df=df, Proportionality.Pvalue=p)]

  # collate whole model stats
  wm <- data.table(C.index=cindex, C.SE=cindex.se, C.L95=cindex.l95, C.U95=cindex.u95,
                   Global.Proportionality.chisq=cz[coefficient == "GLOBAL", chisq],
                   Global.Proportionality.df=cz[coefficient == "GLOBAL", df],
                   Global.Proportionality.Pvalue=cz[coefficient == "GLOBAL", p],
                   Samples=data[cc, .N], Cases=sum(data[cc, ..event_col][[1]]))
  return(cbind(dt, wm))
}
