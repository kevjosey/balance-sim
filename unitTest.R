###############################################################
## TITLE: unitTest.R                                         ##     
## PURPOSE: Test run simulations before using simScen.R      ##
###############################################################

set.seed(10)

# dependencies
library(survey)

# functions for fitting propensity scores
library(twang)
library(CBPS)
library(ATE)
library(cbal)

# additional functions
source("E:/Github/cov-bal-sim/simFUN.R")

iter <- 20
n <- 200
rho <- -0.3
sig2 <- 5
y_scen <- "a"
z_scen <- "a"
tau <- 20

# Haggstrom or Setoguchi

simDat <- replicate(iter, gen_data(n = n, rho = rho, tau = tau, sig2 = sig2,
                                   y_scen = y_scen, z_scen = z_scen))

start <- Sys.time()

idx <- 1:iter
weightsList <- sapply(idx, simFit, simDat = simDat)

stop <- Sys.time()

tauHat_tmp <- t(sapply(idx, simPerf, weightsList = weightsList, simDat = simDat, tau = tau, type = "point"))
tauHat <- apply(tauHat_tmp, 2, mean)

tauSE_tmp <- t(sapply(idx, simPerf, weightsList = weightsList, simDat = simDat, tau, type = "se"))
tauSE <- apply(tauSE_tmp, 2, mean)

coverageProb_tmp <- t(sapply(idx, simPerf, weightsList = weightsList, simDat = simDat, tau = tau, type = "coverage"))
coverageProb <- apply(coverageProb_tmp, 2, mean)

stop - start
