###############################################################
## TITLE: simScen.R                                          ##     
## PURPOSE: Generates and runs simulation scenarios          ##
###############################################################

# dependencies
library(snow)

set.seed(1992)

### Haggstrom scenarios

n <- c(200, 1000)
rho <- c(-0.3, 0, 0.5)
sig2 <- c(2, 5, 10)
y_scen <- c("a", "b")
z_scen <- c("a", "b")

# simulation scenarios
simConditions <- expand.grid(n, sig2, rho, y_scen, z_scen, stringsAsFactors = FALSE)
names(simConditions) <- c("n", "sig2", "rho", "y_scen", "z_scen")

# eliminate scenarios that have already run
index <- 1:nrow(simConditions)
remove <- c()
index <- index[!(index %in% remove)]

# number of replications
iter <- 1000

## multicore simulation
cl <- makeCluster(6, type = "SOCK")

clusterEvalQ(cl, {
  
  # functions for fitting propensity scores
  library(cbal)
  library(CBPS)
  library(ATE)
  library(survey)
  
  # additional functions
  source("~/Github/cbal-sim/simFUN.R")

})

clusterExport(cl = cl, list = list("simConditions", "iter"), envir = environment())

start <- Sys.time()

clusterApply(cl, index, function(i,...) {
  
  dat <- simConditions[i,]
  
  n <- dat$n
  sig2 <- dat$sig2
  rho <- dat$rho
  y_scen <- dat$y_scen
  z_scen <- dat$z_scen
  tau <- 20
  
  simDat <- replicate(iter, ks_data(n = n, tau = tau, sig2 = sig2, rho = rho, 
                                    y_scen = y_scen, z_scen = z_scen))
  
  datFilename <- paste("~/Dropbox (ColoradoTeam)/JoseyDissertation/Data/cbal/simData/",
                       i, n, sig2, rho, y_scen, z_scen, ".RData", sep = "_")
  save(simDat, file = datFilename)
  
  idx <- 1:iter # simulation iteration index
  
  estList <- sapply(idx, simFit, tau = tau, simDat = simDat)
  
  misc_out <- data.frame(tau = rep(tau, times = iter),
                         n = rep(n, times = iter),
                         sig2 = rep(sig2, times = iter),
                         rho = rep(rho, times = iter), 
                         y_scen = rep(y_scen, times = iter),
                         z_scen = rep(z_scen, times = iter),
                         stringsAsFactors = FALSE)

  tauHat_tmp <- do.call(rbind, estList[1,])
  tauSE_tmp <- do.call(rbind, estList[2,])
  coverageProb_tmp <- do.call(rbind, estList[3,])
  colnames(tauHat_tmp) <- colnames(tauSE_tmp) <- colnames(coverageProb_tmp) <- c("CBPS", "SENT", "CAL", "ENT", "BENT")
  
  tauHat <- data.frame(misc_out, tauHat_tmp, stringsAsFactors = FALSE)
  tauSE <- data.frame(misc_out, tauSE_tmp, stringsAsFactors = FALSE)
  coverageProb <- data.frame(misc_out, coverageProb_tmp, stringsAsFactors = FALSE)
  
  tauFilename <- paste("~/Dropbox (ColoradoTeam)/JoseyDissertation/Data/cbal/tauHat/",
                       n, sig2, rho, y_scen, z_scen, ".RData", sep = "_")
  seFilename <- paste("~/Dropbox (ColoradoTeam)/JoseyDissertation/Data/cbal/tauSE/",
                      n, sig2, rho, y_scen, z_scen, ".RData", sep = "_")
  coverageFilename <- paste("~/Dropbox (ColoradoTeam)/JoseyDissertation/Data/cbal/coverageProb/",
                            n, sig2, rho, y_scen, z_scen, ".RData", sep = "_")
  
  save(tauHat, file = tauFilename)
  save(tauSE, file = seFilename)
  save(coverageProb, file = coverageFilename)
  
} )

stopCluster(cl)

stop <- Sys.time()
stop - start
