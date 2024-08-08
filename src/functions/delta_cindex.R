library(survival)

deltaC.coxph <- function(model1, model2, data, time, event) {
  # Estimates the change in C-index between two coxph model fits: `model1`
  # (full model) and `model2` (reduced model). Estimated using the `data` data
  # frame containing column names `time` and `event`.
  
  data <- data[!is.na(data[[time]]) & !is.na(data[[event]]),]
  
  # Calculate predicted responses
  pred1 <- -predict(model1)
  pred2 <- -predict(model2)
  
  # Contrast vector for model1-model2
  contrast1 <- c(1, -1)
  
  # Calculate concordance
  surv.object <- Surv(time = data[[time]], event = data[[event]])
  concord <- concordance(surv.object ~ pred1 + pred2, data)
  
  # Calculate Delta C-index and 95% CI
  diff <- as.numeric(contrast1 %*% coef(concord))
  sd <- as.numeric(sqrt(contrast1 %*% vcov(concord) %*% contrast1))
  l95 = diff - qnorm(0.975) * sd
  u95 = diff + qnorm(0.975) * sd
  p.val <- pmin(1, pnorm(abs(diff / sd), lower.tail = FALSE) * 2)
  
  data.table(deltaC=diff, deltaC.SE=sd, deltaC.L95=l95, deltaC.U95=u95, deltaC.pval=p.val)
}

deltaC.score <- function(survy, lp1, lp2, strata.term) {
  if (!missing(strata.term)) {
    c_compare <- concordance(survy ~ lp1 + lp2 + strata(strata.term))
  } else {
    c_compare <- concordance(survy ~ lp1 + lp2)
  }
  deltaC <- as.vector(c(1, -1) %*% coef(c_compare))
  SE <- as.vector(sqrt(c(1, -1) %*% vcov(c_compare) %*% c(1, -1)))
  L95 <- deltaC - qnorm(0.975) * SE
  U95 <- deltaC + qnorm(0.975) * SE
  pval <- pmin(1, pnorm(abs(deltaC / SE), lower.tail = FALSE) * 2)

  data.table(deltaC=deltaC, deltaC.SE=SE, deltaC.L95=L95, deltaC.U95=U95, deltaC.pval=pval)
}

