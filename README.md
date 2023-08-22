Simulation Study Comparing Covariate Balance Methods
====================================================

This repository contains code for reproducing the simulation study examined in 'A framework for covariate balance using Bregman distances' (Josey et al., 2021). In this numerical experiment, we test and compare different methods for balancing covariate moments between treatment groups to estimate the average treatment effect. These methods include calibration approaches using code found in the github repo [`cbal`](https://github.com/kevjosey/cbal/) and the GMM based covariate balance propensity score approaches contained in the [`cbps`](https://cran.r-project.org/web/packages/CBPS/index.html) package available on CRAN, among others. 

## R Scripts

- [`main.R`](https://github.com/kevjosey/balance-sim/tree/master/main.R): Main script for executing simulations. The primary methods we test include calibration weights such as entropy balancing (Haimueller, 2012), generalized method of moments estimators like covariate balance propensity score methods (Imai & Ratkovic, 2014), along with the classical doubly-robust estimators of the average treatment effect (Bang & Robins, 2005).
- [`output.R`](https://github.com/kevjosey/balance-sim/tree/master/output.R): Code for reproducing the plots and tables (i.e. the simulation results) found in Josey et al. (2021).
- [`simFun.R`](https://github.com/kevjosey/balance-sim/tree/master/simFun.R): Contains functions for generating data and for fitting the different estimates. Data constructs include both homogeneous treatment effects from Kang and Schaeffer (2007) as well as heterogeneous data (which includes effect modifiers).
- [`unitTest.R`](https://github.com/kevjosey/balance-sim/tree/master/unitTest.R): Unit test the accuracy and functionality of the scripts in simFun.R

## References

Bang, H., & Robins, J. M. (2005). Doubly robust estimation in missing data and causal inference models. Biometrics, 61(4), 962-973.

Hainmueller, J. (2012). Entropy balancing for causal effects: A multivariate reweighting method to produce balanced samples in observational studies. Political analysis, 20(1), 25-46.

Imai, K., & Ratkovic, M. (2014). Covariate balancing propensity score. Journal of the Royal Statistical Society Series B: Statistical Methodology, 76(1), 243-263.

Josey, K. P., Juarez‐Colunga, E., Yang, F., & Ghosh, D. (2021). A framework for covariate balance using Bregman distances. Scandinavian Journal of Statistics, 48(3), 790-816.

Kang, J. D., & Schafer, J. L. (2007). Demystifying double robustness: A comparison of alternative strategies for estimating a population mean from incomplete data.