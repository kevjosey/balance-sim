###############################################################
## TITLE: unitTest.R                                         ##     
## PURPOSE: Test run simulations before using simScen.R      ##
###############################################################

set.seed(10)

# dependencies
library(survey)

# functions for fitting propensity scores
library(CBPS)
library(ATE)
library(cbal)

# additional functions
source("~/Github/cbal-sim/simFUN.R")

iter <- 20
n <- 200
rho <- -0.3
sig2 <- 5
y_scen <- "a"
z_scen <- "a"
tau <- 20

# Haggstrom or Setoguchi

simDat <- replicate(iter, ks_data(n = n, rho = rho, tau = tau, sig2 = sig2, y_scen = y_scen, z_scen = z_scen))

start <- Sys.time()

idx <- 1:iter
estList <- sapply(idx, simFit, tau = tau, simDat = simDat)

stop <- Sys.time()

tauHat_tmp <- do.call(rbind, estList[1,])
tauHat <- apply(tauHat_tmp, 2, mean)

tauSE_tmp <- do.call(rbind, estList[2,])
tauSE <- apply(tauSE_tmp, 2, mean)

coverageProb_tmp <- do.call(rbind, estList[3,])
coverageProb <- apply(coverageProb_tmp, 2, mean)

stop - start
