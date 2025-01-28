################################################################################
## Simulatie sample size and power PRECISe-ENERGY study
## Sander van Kuijk
##
## Start: 28-11-2024
################################################################################

library(MASS)
library(reshape2)
library(nlme)
library(matrixcalc)
library(Matrix)
library(lqmm)


# Sample characteristics based on PRECISe study
setwd("L:/SPEC/ICU/RESEARCH/PRECISe/Analyses")
load("PRECISeobjects_V1.0.RData")
rm(list=setdiff(ls(), "p"))

# Correlation between 30 and 90 days estimated from PRECISe, correlation between
# 30 and 60 days, and 60 and 90 days assuming an AR(1) model.

# Simulation parameters
iters <- 500
c03 <- round(cor(p$EQ5D.HUS.proxy.0, p$EQ5D.HUS.imp0.30, use = "complete.obs"), 2)
c06 <- round(cor(p$EQ5D.HUS.proxy.0, p$EQ5D.HUS.imp0.90, use = "complete.obs"), 2)
c09 <- round(cor(p$EQ5D.HUS.proxy.0, p$EQ5D.HUS.imp0.90, use = "complete.obs"), 2)

c36 <- round(sqrt(cor(p$EQ5D.HUS.imp0.30, p$EQ5D.HUS.imp0.90, use = "complete.obs")), 2)
c39 <- round(cor(p$EQ5D.HUS.imp0.30, p$EQ5D.HUS.imp0.90, use = "complete.obs"), 2)

c69 <- round(sqrt(cor(p$EQ5D.HUS.imp0.30, p$EQ5D.HUS.imp0.90, use = "complete.obs")), 2)

mu00 <- 0.78 # Ref PRECISe, Lancet
mu30 <- 0.33 # Ref PRECISe, Lancet
mu90 <- 0.38 # Ref PRECISe, Lancet
mu60 <- (mu30 + mu90)/2 # Interpolation

sd00 <- 0.25
sd30 <- 0.33
sd90 <- 0.38
sd60 <- (sd30 + sd90)/2

mu    <- c(mu00, mu30, mu60, mu90)
sd    <- c(sd00, sd30, sd60, sd90)
es    <- 0.06

# To vary:
n <- 950

# Compute covariance as correlation*sd1*sd2
sig   <- matrix(c(sd[1]^2,       c03*sd00*sd30, c06*sd00*sd60, c09*sd00*sd90,
                  c03*sd00*sd30, sd[2]^2,       c36*sd30*sd60, c39*sd30*sd90,
                  c06*sd00*sd60, c36*sd30*sd60, sd[3]^2,       c69*sd60*sd90,
                  c09*sd00*sd90, c39*sd30*sd90, c69*sd60*sd90, sd[4]^2),
                nrow = 4, byrow = TRUE)
is.positive.definite(sig)

# Dubbelcheck simulatie data
mean30 <- rep(NA, iters)
mean60 <- rep(NA, iters)
mean90 <- rep(NA, iters)
std30  <- rep(NA, iters)
std60  <- rep(NA, iters)
std90  <- rep(NA, iters)
c3060  <- rep(NA, iters)
c3090  <- rep(NA, iters)
c6090  <- rep(NA, iters)

# Resultaten opslaan
beta    <- rep(NA, iters)
pval    <- rep(NA, iters)
  
# Reproduceerbaarheid
set.seed(7181)

# Simulatie
for(i in 1:iters){
  
  Xi <- mvrnorm(n, mu = mu, Sigma = sig, empirical = TRUE)
  Xi[, 1] <- Xi[, 1]
  Xi[, 2] <- Xi[, 2]
  Xi[, 3] <- Xi[, 3]
  Xi[, 4] <- Xi[, 4]
  
  d  <- data.frame(ID = seq(1:n), Xi)
  dl <- melt(d, id.vars = c("ID", "X1"), measure.vars = c("X2", "X3", "X4"),
           variable.name = "fu_moment", value.name = "qol")
  names(dl)[2] <- "baseline"
  dl <- dl[order(dl$ID), ]
  dl$fu_t <- as.numeric(rep(c(30, 60, 90), n))
  dl$rand <- c(rep(0, 0.5*length(dl$ID)), rep(1, 0.5*length(dl$ID)))
  dl$qol  <- ifelse(dl$rand == 1, dl$qol + 0.5*es, dl$qol - 0.5*es)
  
  mean30[i] <- mean(dl$qol[dl$fu_t == 30])
  mean60[i] <- mean(dl$qol[dl$fu_t == 60])
  mean90[i] <- mean(dl$qol[dl$fu_t == 90])
  std30[i]  <- sd(dl$qol[dl$fu_t == 30])
  std60[i]  <- sd(dl$qol[dl$fu_t == 60])
  std90[i]  <- sd(dl$qol[dl$fu_t == 90])
  c3060[i]  <- cor(dl$qol[dl$fu_t == 30], dl$qol[dl$fu_t == 60])
  c3090[i]  <- cor(dl$qol[dl$fu_t == 30], dl$qol[dl$fu_t == 90])
  c6090[i]  <- cor(dl$qol[dl$fu_t == 60], dl$qol[dl$fu_t == 90])

  fit <- lme(qol ~ rand + baseline, data = dl, rand = ~ 1|ID,
             correlation = corAR1(form = ~ fu_t | ID),
             control = list(opt = "optim"))

  beta[i] <- fixef(fit)[2]
  pval[i] <- summary(fit)$tTable[2, 5]
  
}

# Simulatie check
mean(mean30)
mean(mean60)
mean(mean90)
mean(std30)
mean(std60)
mean(std90)
mean(c3060)
mean(c3090)
mean(c6090)

# Resultaten
noquote(paste("Power to detect a difference of ", es, " points is ",
               (sum(pval <= 0.05)/iters)*100, "%.", sep = ""))
round(mean(beta), 3)

### Einde file.