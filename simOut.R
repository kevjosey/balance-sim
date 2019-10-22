###############################################################
## TITLE: simOut.R                                           ##     
## PURPOSE: Construction of tables and plots from outputs    ##
###############################################################

library(ggplot2)
library(gridExtra)

### Kang and Schafer results

dir_1 <- "D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/cbal/ATE/tauHat/"
dir_2 <- "D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/cbal/ATE/tauSE/"
dir_3 <- "D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/cbal/ATE/coverageProb/"

files <- list.files(dir_1)
out_1 <- out_2 <- out_3 <- matrix("", nrow = length(files), ncol = 10)

colnames(out_1) <- colnames(out_2) <- colnames(out_3) <- 
  c("tau", "n", "variance", "correlation", "y_scen", "z_scen", "IPW", "CBPS", "SENT", "BENT")

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

load("D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/cbal/ATE/tauHat/_200_10_0_a_a_.RData")
dat1 <- stack(as.data.frame(tauHat[,7:ncol(tauHat)]))
dat1$ind <- factor(dat1$ind, labels = c("IPW", "CBPS", "SENT", "BENT"))
p1 <- ggplot(dat1) + 
  geom_boxplot(aes(x = ind, y = values, fill = ind)) + 
  ylab("tau") +
  ylim(0, 40) +
  xlab("") +
  ggtitle("outcome scenario: a, treatment scenario: a") +
  guides(fill =  FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

load("D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/cbal/ATE/tauHat/_200_10_0_a_b_.RData")
dat2 <- stack(as.data.frame(tauHat[,7:ncol(tauHat)]))
dat2$ind <- factor(dat2$ind, labels = c("IPW", "CBPS", "SENT", "BENT"))
p2 <- ggplot(dat2) + 
  geom_boxplot(aes(x = ind, y = values, fill = ind)) + 
  ylab("tau") +
  ylim(0, 40) +
  xlab("") +
  ggtitle("outcome scenario: a, treatment scenario: b") +
  guides(fill =  FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

load("D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/cbal/ATE/tauHat/_200_10_0_b_a_.RData")
dat3 <- stack(as.data.frame(tauHat[,7:ncol(tauHat)]))
dat3$ind <- factor(dat3$ind, labels = c("IPW", "CBPS", "SENT", "BENT"))
p3 <- ggplot(dat3) + 
  geom_boxplot(aes(x = ind, y = values, fill = ind)) + 
  ylab("tau") +
  ylim(0, 40) +
  xlab("") +
  ggtitle("outcome scenario: b, treatment scenario: a") +
  guides(fill =  FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

load("D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/cbal/ATE/tauHat/_200_10_0_b_b_.RData")
dat4 <- stack(as.data.frame(tauHat[,7:ncol(tauHat)]))
dat4$ind <- factor(dat4$ind, labels = c("IPW", "CBPS", "SENT", "BENT"))
p4 <- ggplot(dat4) + 
  geom_boxplot(aes(x = ind, y = values, fill = ind)) + 
  ylab("tau") +
  ylim(0, 40) +
  xlab("") +
  ggtitle("outcome scenario: b, treatment scenario: b") +
  guides(fill =  FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
  
png("D:/Dropbox (ColoradoTeam)/JoseyDissertation/Output/cbal/Figures/ATE_plot.png", 
    width = 1000, 
    height = 1000,
    res = 100, 
    units = "px")

grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
  
dev.off()

# Generate better tables

out_1_clean <- as.data.frame(out_1)
out_1_clean$n <- as.numeric(out_1_clean$n)
out_1_clean$variance <- as.numeric(out_1_clean$variance)
out_1_clean$correlation <- as.numeric(out_1_clean$correlation)
out_1_clean <- out_1_clean[with(out_1_clean, order(n, variance, correlation, y_scen, z_scen)),]

out_2_clean <- as.data.frame(out_2)
out_2_clean$n <- as.numeric(out_2_clean$n)
out_2_clean$variance <- as.numeric(out_2_clean$variance)
out_2_clean$correlation <- as.numeric(out_2_clean$correlation)
out_2_clean <- out_2_clean[with(out_2_clean, order(n, variance, correlation, y_scen, z_scen)),]

out_3_clean <- as.data.frame(out_3)
out_3_clean$n <- as.numeric(out_3_clean$n)
out_3_clean$variance <- as.numeric(out_3_clean$variance)
out_3_clean$correlation <- as.numeric(out_3_clean$correlation)
out_3_clean <- out_3_clean[with(out_3_clean, order(n, variance, correlation, y_scen, z_scen)),]
  
filename1 <- "D:/Dropbox (ColoradoTeam)/JoseyDissertation/Output/cbal/Tables/estimates_ATE.csv"
filename2 <- "D:/Dropbox (ColoradoTeam)/JoseyDissertation/Output/cbal/Tables/coverageProbs_ATE.csv"
filename3 <- "D:/Dropbox (ColoradoTeam)/JoseyDissertation/Output/cbal/Tables/mse_bias_ATE.csv"

write.csv(out_1_clean, file = filename1)
write.csv(out_2_clean, file = filename2)
write.csv(out_3_clean, file = filename3)

### HTE results

dir_1 <- "D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/cbal/HTE/tauHat/"
dir_2 <- "D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/cbal/HTE/tauSE/"
dir_3 <- "D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/cbal/HTE/coverageProb/"

files <- list.files(dir_1)
out_1 <- out_2 <- out_3 <- matrix("", nrow = length(files), ncol = 11)

colnames(out_1) <- colnames(out_2) <- colnames(out_3) <- 
  c("tau", "n", "variance", "correlation", "y_scen", "z_scen", "AIPW", "CAL", "iCBPS", "hdCBPS", "SENT")

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

load("D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/cbal/HTE/tauHat/_200_10_0_a_a_.RData")
dat1 <- stack(as.data.frame(tauHat[,7:ncol(tauHat)]))
dat1$ind <- factor(dat1$ind, labels = c("AIPW", "CAL", "iCBPS", "hdCBPS", "SENT"))
p1 <- ggplot(dat1) + 
  geom_boxplot(aes(x = ind, y = values, fill = ind)) + 
  ylab("tau") +
  ylim(0, 40) +
  xlab("") +
  ggtitle("outcome scenario: a, treatment scenario: a") +
  guides(fill =  FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

load("D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/cbal/HTE/tauHat/_200_10_0_a_b_.RData")
dat2 <- stack(as.data.frame(tauHat[,7:ncol(tauHat)]))
dat2$ind <- factor(dat2$ind, labels = c("AIPW", "CAL", "iCBPS", "hdCBPS", "SENT"))
p2 <- ggplot(dat2) + 
  geom_boxplot(aes(x = ind, y = values, fill = ind)) + 
  ylab("tau") +
  ylim(0, 40) +
  xlab("") +
  ggtitle("outcome scenario: a, treatment scenario: b") +
  guides(fill =  FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

load("D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/cbal/HTE/tauHat/_200_10_0_b_a_.RData")
dat3 <- stack(as.data.frame(tauHat[,7:ncol(tauHat)]))
dat3$ind <- factor(dat3$ind, labels = c("AIPW", "CAL", "iCBPS", "hdCBPS", "SENT"))
p3 <- ggplot(dat3) + 
  geom_boxplot(aes(x = ind, y = values, fill = ind)) + 
  ylab("tau") +
  ylim(0, 40) +
  xlab("") +
  ggtitle("outcome scenario: b, treatment scenario: a") +
  guides(fill =  FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

load("D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/cbal/HTE/tauHat/_200_10_0_b_b_.RData")
dat4 <- stack(as.data.frame(tauHat[,7:ncol(tauHat)]))
dat4$ind <- factor(dat4$ind, labels = c("AIPW", "CAL", "iCBPS", "hdCBPS", "SENT"))
p4 <- ggplot(dat4) + 
  geom_boxplot(aes(x = ind, y = values, fill = ind)) + 
  ylab("tau") +
  ylim(0, 40) +
  xlab("") +
  ggtitle("outcome scenario: b, treatment scenario: b") +
  guides(fill =  FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

png("D:/Dropbox (ColoradoTeam)/JoseyDissertation/Output/cbal/Figures/HTE_plot.png", 
    width = 1000, 
    height = 1000,
    res = 100, 
    units = "px")

grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

dev.off()

# Generate better tables

out_1_clean <- as.data.frame(out_1)
out_1_clean$n <- as.numeric(out_1_clean$n)
out_1_clean$variance <- as.numeric(out_1_clean$variance)
out_1_clean$correlation <- as.numeric(out_1_clean$correlation)
out_1_clean <- out_1_clean[with(out_1_clean, order(n, variance, correlation, y_scen, z_scen)),]

out_2_clean <- as.data.frame(out_2)
out_2_clean$n <- as.numeric(out_2_clean$n)
out_2_clean$variance <- as.numeric(out_2_clean$variance)
out_2_clean$correlation <- as.numeric(out_2_clean$correlation)
out_2_clean <- out_2_clean[with(out_2_clean, order(n, variance, correlation, y_scen, z_scen)),]

out_3_clean <- as.data.frame(out_3)
out_3_clean$n <- as.numeric(out_3_clean$n)
out_3_clean$variance <- as.numeric(out_3_clean$variance)
out_3_clean$correlation <- as.numeric(out_3_clean$correlation)
out_3_clean <- out_3_clean[with(out_3_clean, order(n, variance, correlation, y_scen, z_scen)),]

filename1 <- "D:/Dropbox (ColoradoTeam)/JoseyDissertation/Output/cbal/Tables/estimates_HTE.csv"
filename2 <- "D:/Dropbox (ColoradoTeam)/JoseyDissertation/Output/cbal/Tables/coverageProbs_HTE.csv"
filename3 <- "D:/Dropbox (ColoradoTeam)/JoseyDissertation/Output/cbal/Tables/mse_bias_HTE.csv"

write.csv(out_1_clean, file = filename1)
write.csv(out_2_clean, file = filename2)
write.csv(out_3_clean, file = filename3)
