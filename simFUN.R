##################################################################
## TITLE: simFUN.R                                              ##     
## PURPOSE: Functions for simulating and fitting causal models  ##
##################################################################

require(cbal)
require(CBPS)
require(ATE)
require(survey)

# simulate scenarios
ks_data <- function(tau, n, sig2, rho, y_scen = c("a", "b"), z_scen = c("a", "b")) {
  
  # covariates
  x1 <- stats::rnorm(n, 0, 1)
  x2 <- stats::rnorm(n, 0, 1)
  x3 <- stats::rnorm(n, 0, 1)
  x4 <- stats::rnorm(n, 0, 1)
  
  # transformed predictors
  u1 <- as.numeric(scale(exp(x1/2)))
  u2 <- as.numeric(scale(x2/(1 + exp(x1)) + 10))
  u3 <- as.numeric(scale((x1*x3/25 + 0.6)^3))
  u4 <- as.numeric(scale((x2 + x4 + 20)^2))
  
  # treatment probabilities
  if (z_scen == "b")
    e_X <- 1/(1 + exp( -(-u1 + 0.5*u2 - 0.25*u3 - 0.1*u4) ) )
  else
    e_X <- 1/(1 + exp( -(-x1 + 0.5*x2 - 0.25*x3 - 0.1*x4) ) )
  
  r_exp <- stats::runif(n)
  z <- ifelse(r_exp < e_X, 1, 0)
  
  # error variance
  R <- matrix(rho, nrow = 2, ncol = 2)
  diag(R) <- 1
  V <- diag(sqrt(sig2), nrow = 2, ncol = 2)
  Sig <- V %*% R %*% V
  
  if (y_scen == "b")
    mu <- 210 + 27.4*u1 + 13.7*u2 + 13.7*u3 + 13.7*u4
  else
    mu <- 210 + 27.4*x1 + 13.7*x2 + 13.7*x3 + 13.7*x4
  
  eval <- eigen(Sig, symmetric = TRUE)
  y_init <- matrix(stats::rnorm(n*2, 0, 1), nrow = n, ncol = 2) # iid potential outcomes
  y_tmp <- t(eval$vectors %*% diag(sqrt(eval$values), nrow = 2) %*% t(y_init)) # SVD
  y_pot <- y_tmp + cbind(mu, mu + tau) # include causal effect
  
  # observed outcome
  y <- z*y_pot[,2] + (1 - z)*y_pot[,1]
  
  # create simulation dataset
  sim <- as.data.frame(cbind(y, z, x1, x2, x3, x4, u1, u2, u3, u4))
  
  return(sim)
  
}

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
simFit <- function(idx = 1, simDat) {
  
  dat <- simDat[,idx]
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
simPerf <- function(idx = 1, weightsList, simDat, tau, type = c("point", "se", "coverage")) {
  
  data <- as.data.frame(simDat[,idx])
  work <- as.data.frame(weightsList[,idx])
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
