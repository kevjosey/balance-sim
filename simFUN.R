##################################################################
## TITLE: simFUN.R                                              ##     
## PURPOSE: Functions for simulating and fitting causal models  ##
##################################################################

require(cbal)
require(CBPS)
require(ATE)
require(survey)

# function for calculating inverse probabilty of treatment weights
ipw <- function(ps, treat, estimand = c("ATE", "ATT"), standardize = TRUE) {
  
  if (estimand == "ATE")
    w <- treat/ps + (1 - treat)/(1-ps)
  else
    w <- treat + ps*(1-treat)/(1-ps)
  
  if (standardize) { 
    w[treat == 0] <- w[treat == 0]/sum(w[treat == 0])
    w[treat == 1] <- w[treat == 1]/sum(w[treat == 1])
  }
  
  return(w)
  
}

# Fits the balancing weights using a variety of methods
simFit <- function(iter = 1, simDat) {
  
  dat <- simDat[,iter]
  formula <- as.formula(z ~ x1 + x2 + x3 + x4, env = environment(dat))
  cov_dat <- as.data.frame(dat[c("x1", "x2", "x3", "x4")])
  
  # cbps
  fit_cbps <- CBPS(formula, data = dat, ATT = 0, method = "exact", verbose = FALSE)
  wts_cbps <- fit_cbps$weights
  
  # sent
  fit_sent <- cbalance(formula, data = dat, distance = "shifted", estimand = "ATE")
  wts_sent <- fit_sent$weights
  
  # ate
  fit_ate <- ATE(Y = dat$y, Ti = dat$z, X = cov_dat, theta = 0, ATT = FALSE)
  wts_ate <- fit_ate$weights.p + fit_ate$weights.q

  # ent
  fit_ent <- cbalance(formula, data = dat, distance = "entropy", estimand = "ATE")
  wts_ent <- fit_ent$weights
  
  # bent
  fit_bent <- cbalance(formula, data = dat, distance = "binary", estimand = "ATE")
  wts_bent <- fit_bent$weights

  # results
  out <- as.data.frame(cbind(wts_cbps, wts_sent, wts_ate, wts_ent, wts_bent))
  
  return(out)
  
}

# wrapper for simEstimate()
simPerf <- function(iter = 1, weightsList, simDat, tau, type = c("point", "se", "coverage")) {
  
  data <- as.data.frame(simDat[,iter])
  work <- as.data.frame(weightsList[,iter])
  formula <- y ~ z
  
  est_tmp <- lapply(work, simEstimate, formula = formula, data = data, tau = tau, type = type)
  
  est_out <- do.call(c, est_tmp)
  names(est_out) <- names(work)
  
  return(est_out)
  
}

# finds average causal effect estimates
simEstimate <- function(formula, data, weights, tau, type = c("point", "se", "coverage")) {
  
  design <- svydesign(ids = ~ 1, weights = ~ weights, data = data.frame(weights, data))
  fit <- svyglm(formula, design = design, family = gaussian)
  
  if (type == "se")
    est <- SE(fit)[2]
  else if (type == "coverage")
    est <- as.numeric(confint(fit)[2,1] <= tau & confint(fit)[2,2] >= tau)
  else # type = "point"
    est <- coef(fit)[2]
  
  return(est)
  
}
