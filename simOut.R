###############################################################
## TITLE: simOut.R                                           ##     
## PURPOSE: Construction of tables and plots from outputs    ##
###############################################################

library(ggplot2)
library(gridExtra)

### Kang and Schafer results

dir_1 <- "E:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/tauHat/"
dir_2 <- "E:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/tauSE/"
dir_3 <- "E:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/coverageProb/"

files <- list.files(dir_1)
out_1 <- out_2 <- out_3 <- matrix("", nrow = length(files), ncol = 11)

colnames(out_1) <- colnames(out_2) <- colnames(out_3) <- 
  c("tau", "n", "variance", "correlation", "y_scen", "z_scen",
     "CBPS", "SENT", "CAL", "ENT", "BENT")

j <- 1

for (fdx in files) {
  
  file_est <- paste0(dir_1, fdx)
  file_se <- paste0(dir_2, fdx)
  file_coverage <- paste0(dir_3, fdx)
  load(file_est)
  load(file_se)
  load(file_coverage)
  
  lbl <- do.call(c, lapply(tauHat[1,1:6], as.character))
  
  out_1[j,1:6] <- out_2[j,1:6] <- out_3[j,1:6] <- lbl
  
  est <- apply(tauHat[,7:ncol(tauHat)], 2, mean)
  se <- apply(tauHat[,7:ncol(tauHat)], 2, sd)
  est_se <- sapply(1:length(est), function(i,...) 
    paste(round(est[i], 2), " (", round(se[i], 2), ")", sep = ""))
  
  mse <- apply(tauHat[,7:ncol(tauHat)], 2,
               function(x, tau) mean((x - tau)^2), 
               tau = tauHat$tau)
  bias <- apply(tauHat[,7:ncol(tauHat)], 2,
                function(x, tau) mean((x - tau)), 
                tau = tauHat$tau)
  mse_bias <- sapply(1:length(mse), function(i,...) 
    paste(round(mse[i], 2), " (", round(bias[i], 2), ")", sep = ""))
  
  coverage <- apply(coverageProb[,7:ncol(coverageProb)], 2, mean)
  
  out_1[j,7:ncol(out_1)] <- est_se
  out_2[j,7:ncol(out_2)] <- coverage
  out_3[j,7:ncol(out_3)] <- mse_bias
  
  j <- j + 1
  
}

# plot outcomes

load("E:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/tauHat/_200_10_0_a_a_.RData")
dat1 <- stack(as.data.frame(tauHat[,7:ncol(tauHat)]))
dat1$ind <- factor(dat1$ind, labels = c("CBPS", "SENT", "CAL", "ENT", "BENT"))
p1 <- ggplot(dat1) + 
  geom_boxplot(aes(x = ind, y = values, fill = ind)) + 
  ylab("tau") +
  ylim(0, 40) +
  xlab("") +
  ggtitle("outcome scenario: a, treatment scenario: a") +
  guides(fill =  FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

load("E:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/tauHat/_200_10_0_a_b_.RData")
dat2 <- stack(as.data.frame(tauHat[,7:ncol(tauHat)]))
dat2$ind <- factor(dat2$ind, labels = c("CBPS", "SENT", "CAL", "ENT", "BENT"))
p2 <- ggplot(dat2) + 
  geom_boxplot(aes(x = ind, y = values, fill = ind)) + 
  ylab("tau") +
  ylim(0, 40) +
  xlab("") +
  ggtitle("outcome scenario: a, treatment scenario: b") +
  guides(fill =  FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

load("E:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/tauHat/_200_10_0_b_a_.RData")
dat3 <- stack(as.data.frame(tauHat[,7:ncol(tauHat)]))
dat3$ind <- factor(dat3$ind, labels = c("CBPS", "SENT", "CAL", "ENT", "BENT"))
p3 <- ggplot(dat3) + 
  geom_boxplot(aes(x = ind, y = values, fill = ind)) + 
  ylab("tau") +
  ylim(0, 40) +
  xlab("") +
  ggtitle("outcome scenario: b, treatment scenario: a") +
  guides(fill =  FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

load("E:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/tauHat/_200_10_0_b_b_.RData")
dat4 <- stack(as.data.frame(tauHat[,7:ncol(tauHat)]))
dat4$ind <- factor(dat4$ind, labels = c("CBPS", "SENT", "CAL", "ENT", "BENT"))
p4 <- ggplot(dat4) + 
  geom_boxplot(aes(x = ind, y = values, fill = ind)) + 
  ylab("tau") +
  ylim(0, 40) +
  xlab("") +
  ggtitle("outcome scenario: b, treatment scenario: b") +
  guides(fill =  FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
  
png("E:/Dropbox (ColoradoTeam)/JoseyDissertation/Output/Figures/ATE_plot.png", 
    width = 1000, 
    height = 1000,
    res = 100, 
    units = "px")

grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
  
dev.off()

# Generate better tables

out_1.tmp <- subset(data.frame(out_1, stringsAsFactors = FALSE)) # n == "200" & variance == "10" & correlation == "0", 
                   # select = -c(tau, n, variance, correlation))
out_1.tmp$n <- as.numeric(out_1.tmp$n)
out_1.tmp$variance <- as.numeric(out_1.tmp$variance)
out_1.tmp$correlation <- as.numeric(out_1.tmp$correlation)
out_1.tmp <- out_1.tmp[order(out_1.tmp$n, out_1.tmp$variance, out_1.tmp$correlation, out_1.tmp$y_scen, out_1.tmp$z_scen),]

out_2.tmp <- subset(data.frame(out_2, stringsAsFactors = FALSE)) # n == "200" & variance == "10" & correlation == "0", 
                    # select = -c(tau, n, variance, correlation))
out_2.tmp$n <- as.numeric(out_2.tmp$n)
out_2.tmp$variance <- as.numeric(out_2.tmp$variance)
out_2.tmp$correlation <- as.numeric(out_2.tmp$correlation)
out_2.tmp <- out_2.tmp[order(out_2.tmp$n, out_2.tmp$variance, out_2.tmp$correlation, out_2.tmp$y_scen, out_2.tmp$z_scen),]

out_3.tmp <- subset(data.frame(out_3, stringsAsFactors = FALSE)) # n == "200" & variance == "10" & correlation == "0", 
                    # select = -c(tau, n, variance, correlation))
out_3.tmp$n <- as.numeric(out_3.tmp$n)
out_3.tmp$variance <- as.numeric(out_3.tmp$variance)
out_3.tmp$correlation <- as.numeric(out_3.tmp$correlation)
out_3.tmp <- out_3.tmp[order(out_3.tmp$n, out_3.tmp$variance, out_3.tmp$correlation, out_3.tmp$y_scen, out_3.tmp$z_scen),]
  
filename1 <- "E:/Dropbox (ColoradoTeam)/JoseyDissertation/Output/Tables/estimates.csv"
filename2 <- "E:/Dropbox (ColoradoTeam)/JoseyDissertation/Output/Tables/coverageProbs.csv"
filename3 <- "E:/Dropbox (ColoradoTeam)/JoseyDissertation/Output/Tables/mse_bias.csv"

write.csv(out_1.tmp, file = filename1)
write.csv(out_2.tmp, file = filename2)
write.csv(out_3.tmp, file = filename3)
