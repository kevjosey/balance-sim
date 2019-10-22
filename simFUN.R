##################################################################
## TITLE: simFUN.R                                              ##     
## PURPOSE: Functions for simulating and fitting causal models  ##
##################################################################

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
  sim <- as.data.frame(cbind(y, z, x1, x2, x3, x4))
  
  return(sim)
  
}

hte_data <- function(n, sig2, rho, y_scen = c("a", "b"), z_scen = c("a", "b")){
  
  # error variance
  R <- matrix(rho, nrow = 2, ncol = 2)
  diag(R) <- 1
  V <- diag(sqrt(sig2), nrow = 2, ncol = 2)
  Sig <- V %*% R %*% V
  
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
  
  # effect coefficients
  beta <- c(210, 27.4, 13.7, 13.7, 13.7)
  gamma <- c(20, -13.7, 0, 0, 13.7)
  
  # propensity score
  if (z_scen == "b") {
    e_X <- 1/(1 + exp( -(-u1 + 0.5*u2 - 0.25*u3 - 0.1*u4 ) ) )
  } else { # z_scen == "a"
    e_X <- 1/(1 + exp( -(-x1 + 0.5*x2 - 0.25*x3 - 0.1*x4 ) ) )
  }
  
  z <- rbinom(n, 1, e_X)
  
  if (y_scen == "b") {
    X <- cbind(rep(1, times = n), u1, u2, u3, u4)
  } else { # y_scen == "b"
    X <- cbind(rep(1, times = n), x1, x2, x3, x4)
  }
  
  # outcome mean
  mu_0 <- X%*%beta
  mu_1 <- X%*%(beta + gamma)
  
  tau <- t(apply(X[z == 1,], 2, mean)) %*% gamma
  
  # potential outcomes
  eval <- eigen(Sig, symmetric = TRUE)
  y_init <- matrix(stats::rnorm(n*2, 0, 1), nrow = n, ncol = 2) # iid potential outcomes
  y_tmp <- t(eval$vectors %*% diag(sqrt(eval$values), nrow = 2) %*% t(y_init)) # SVD
  y_pot <- y_tmp + cbind(mu_0, mu_1) # include causal effect
  
  # observed outcome
  y <- z*y_pot[,2] + (1 - z)*y_pot[,1]
  
  # create simulation dataset
  sim <- list(y = y, z = z, x1 = x1, x2 = x2, x3 = x3, x4 = x4)
  
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
simFit_ATE <- function(idx = 1, simDat, tau) {
  
  dat <- simDat[,idx]
  formula <- as.formula(z ~ x1 + x2 + x3 + x4, env = environment(dat))
  Y <- dat$y
  Z <- dat$z
  
  # ipw
  fit_glm <- glm(formula, data = dat, family = binomial(link = "logit"))
  ps <- predict(fit_glm, type = "response")
  wts_glm <- ipw(ps, treat = Z, estimand = "ATE", standardize = FALSE)
  design_glm <- svydesign(ids = ~ 1, weights = ~ wts_glm, data = data.frame(wts_glm, dat))
  mod_glm <- svyglm(y ~ z, design = design_glm, family = gaussian)
  tau_glm <- coef(mod_glm)[2]
  se_glm <- SE(mod_glm)[2]
  cp_glm <- as.numeric(confint(mod_glm)[2,1] <= tau & confint(mod_glm)[2,2] >= tau)
  
  # cbps
  fit_cbps <- CBPS(formula, data = dat, ATT = 0, method = "exact", verbose = FALSE)
  wts_cbps <- fit_cbps$weights
  design_cbps <- svydesign(ids = ~ 1, weights = ~ wts_cbps, data = data.frame(wts_cbps, dat))
  mod_cbps <- svyglm(y ~ z, design = design_cbps, family = gaussian)
  tau_cbps <- coef(mod_cbps)[2]
  se_cbps <- SE(mod_cbps)[2]
  cp_cbps <- as.numeric(confint(mod_cbps)[2,1] <= tau & confint(mod_cbps)[2,2] >= tau)
  
  # sent
  fit_sent <- cbalance(formula, data = dat, distance = "shifted")
  est_sent <- cestimate(fit_sent, Y = Y, method = "sandwich")
  tau_sent <- est_sent$tau
  se_sent <- sqrt(est_sent$variance)
  cp_sent <- as.numeric(tau_sent - se_sent*1.96 <= tau & tau_sent + se_sent*1.96 >= tau)
  
  # bent
  fit_bent <- cbalance(formula, data = dat, distance = "binary")
  est_bent <- cestimate(fit_bent, Y = Y, method = "sandwich")
  tau_bent <- est_bent$tau
  se_bent <- sqrt(est_bent$variance)
  cp_bent <- as.numeric(tau_bent - se_bent*1.96 <= tau & tau_bent + se_bent*1.96 >= tau)

  # results
  tauh <- c(tau_glm, tau_cbps, tau_sent, tau_bent)
  seh <- c(se_glm, se_cbps, se_sent, se_bent)
  cph <- c(cp_glm, cp_cbps, cp_sent, cp_bent)
  
  out <- list(tauh = tauh, seh = seh, cph = cph)
    
  return(out)
  
}

simFit_HTE <- function(idx = 1, simDat, tau) {
  
  dat <- as.data.frame(simDat[,idx])
  formula <- as.formula(z ~ x1 + x2 + x3 + x4, env = environment(dat))
  outform <- as.formula(y ~ x1 + x2 + x3 + x4, env = environment(dat))
  cov_dat <- data.frame(dat[c("x1", "x2", "x3", "x4")])
  X <- as.matrix(cov_dat)
  Y <- dat$y
  Z <- dat$z
  inform <- as.formula(~ X)
  
  # aipw
  out <- ps <- matrix(nrow = nrow(X), ncol = 2)
  fit_0 <- lm(outform, data = dat[Z == 0,])
  fit_1 <- lm(outform, data = dat[Z == 1,])
  out[,1] <- predict(fit_0, newdata = dat)
  out[,2] <- predict(fit_1, newdata = dat)
  fit_ps <- glm(formula, data = dat, family = binomial(link = "logit"))
  ps[,1] <- 1 - predict(fit_ps, type = "response")
  ps[,2] <- predict(fit_ps, type = "response")
  fit_aipw <- ate.aipw(y = Y, tr = Z, mfp = ps, mfo = out)
  tau_aipw <- fit_aipw$diff.est[2]
  se_aipw <- sqrt(fit_aipw$diff.var[2])
  cp_aipw <- tau_aipw - se_aipw*1.96 <= tau & tau_aipw + se_aipw*1.96 >= tau
  
  # ate
  fit_ate <- ATE(Y = Y, Ti = Z, X = cov_dat, theta = 0, ATT = FALSE)
  sate <- summary(fit_ate)$Estimate
  tau_ate <- sate[3,1]
  se_ate <- sate[3,2]
  cp_ate <- as.numeric(sate[3,3] <= tau & sate[3,4] >= tau)
  
  # icbps
  fit_icbps <- CBPS(formula, data = dat, ATT = FALSE, method = "exact", verbose = FALSE, 
                    diff.formula = inform, baseline.formula = inform)
  wts_icbps <- fit_icbps$weights
  design <- svydesign(ids = ~ 1, weights = ~ wts_icbps, data = data.frame(wts_icbps, dat))
  mod_icbps <- svyglm(y ~ z, design = design, family = gaussian)
  tau_icbps <- coef(mod_icbps)[2]
  se_icbps <- SE(mod_icbps)[2]
  cp_icbps <- as.numeric(confint(mod_icbps)[2,1] <= tau & confint(mod_icbps)[2,2] >= tau)
  
  # hdcbps
  fit_hdcbps <- hdCBPS(formula, data = dat, y = Y, ATT = 0, method = "linear")
  tau_hdcbps <- fit_hdcbps$ATE
  se_hdcbps <- fit_hdcbps$s
  cp_hdcbps <- as.numeric(tau_hdcbps - se_hdcbps*1.96 <= tau & tau_hdcbps + se_hdcbps*1.96 >= tau)
  
  # sent
  fit_sent <- cbalance(formula, data = dat, distance = "shifted")
  est_sent <- cestimate(fit_sent, Y = Y, method = "sandwich")
  tau_sent <- est_sent$tau
  se_sent <- sqrt(est_sent$variance)
  cp_sent <- as.numeric(tau_sent - se_sent*1.96 <= tau & tau_sent + se_sent*1.96 >= tau)
  
  # results
  tauh <- c(tau_aipw, tau_ate, tau_icbps, tau_hdcbps, tau_sent)
  seh <- c(se_aipw, se_ate, se_icbps, se_hdcbps, se_sent)
  cph <- c(cp_aipw, cp_ate, cp_icbps, cp_hdcbps, cp_sent)
  
  out <- list(tauh = tauh, seh = seh, cph = cph)
  
  return(out)
  
}
