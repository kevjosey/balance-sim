Simulation Study Comparing Covariate Balance Methods
====================================================

This repository contains code for reproducing the simulation study examined in 'A framework for covariate balance using Bregman distances' (Josey et al., 2021). In this numerical experiment, we test and compare different methods for balancing covariate moments between treatment groups to estimate the average treatment effect. These methods include calibration approaches using code found in the github repo [`cbal`](https://github.com/kevjosey/balance-sim/tree/master/simFun.R), the GMM based covariate balance propensity scores contained in the [`cbps`](https://cran.r-project.org/web/packages/CBPS/index.html) package available on CRAN, among others. 

Below is an itemized list of the different R scripts within this repository:

- [`main.R`](https://github.com/kevjosey/balance-sim/tree/master/simFun.R): Main script for executing simulations. The primary methods we test include calibration weights such as entropy balancing, generalized method of moments estimators like covariate balance propensity score methods (CBPS), along with the classical doubly-robust estimators of the average treatment effect
- [`output.R`](https://github.com/kevjosey/balance-sim/tree/master/simFun.R): Code for reproducing the plots found in the main manuscript.
- [`simFun.R`](https://github.com/kevjosey/balance-sim/tree/master/simFun.R): Contains functions for generating data and for fitting the different estimates. Data constructs include both homogeneous treatment effects from Kang and Schaeffer (2013) as well as heterogeneous data (which includes effect modifiers).
- [unitTest.R]((https://github.com/kevjosey/balance-sim/tree/master/unitTest.R): Code for testing the accuracy and functionality of the scripts in simFun.R