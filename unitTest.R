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

iter <- 100
n <- 200
rho <- -0.3
sig2 <- 5
y_scen <- "a"
z_scen <- "a"
tau <- 20

# Kang and Schafer

simDat <- replicate(iter, ks_data(n = n, rho = rho, tau = tau, sig2 = sig2, y_scen = y_scen, z_scen = z_scen))

start <- Sys.time()

idx <- 1:iter
estList <- sapply(idx, simFit_ATE, tau = tau, simDat = simDat)

stop <- Sys.time()

apply(do.call(rbind, estList[1,]), 2, mean)
apply(do.call(rbind, estList[2,]), 2, mean)
apply(do.call(rbind, estList[3,]), 2, mean)

stop - start

# HTE

simDat <- replicate(iter, hte_data(n = n, rho = rho, sig2 = sig2, y_scen = y_scen, z_scen = z_scen))

start <- Sys.time()

idx <- 1:iter
estList <- sapply(idx, simFit_HTE, tau = tau, simDat = simDat)

stop <- Sys.time()

apply(do.call(rbind, estList[1,]), 2, mean)
apply(do.call(rbind, estList[2,]), 2, mean)
apply(do.call(rbind, estList[3,]), 2, mean)

stop - start
