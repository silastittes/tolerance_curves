#prior predictive checks
library(truncnorm)
library(gdata)
library(dplyr)

#simulate data set for tolerance model
source("bayes/tolerance_functions.R")
emery <- load_emery()

prior_data <- read.xls(xls = "bayes/EcoLettData_PriorPredChecks.xls") %>%
  filter(Treatment == "NR")


## simulate parameters and data ----------------------------
nSpp <- length(unique(emery$Species))
ntreat <- length(unique(emery$treat))
nobs <- nSpp * ntreat * 7
x_o <- seq(1, ntreat, length.out = ntreat) #data scale
nSppTreat <- nobs/nSpp/ntreat

simreps <- 50
ydat <- matrix(NA, nrow = nobs, ncol = simreps)
for (z in 1:simreps){
  
  nu <- rgamma(1, 20, 0.2)
  a <- rtruncnorm(n = nSpp, mean = 4, sd = 1, a = 1)
  b <- rtruncnorm(n = nSpp, mean = 3, sd = 1, a = 1)
  
  #hyper paramters for c
  mean_c <- rtruncnorm(n = nSpp, mean = 0, sd = 10, a = 0)
  var_c <- rtruncnorm(n = nSpp, mean = 0, sd = 10, a = 0)
  
  c <- rtruncnorm(n = nSpp, mean = mean_c, sd = var_c, a = 0)
  d <- rtruncnorm(n = nSpp, mean = min(x_o), sd = 2, b = min(x_o))
  e <- rtruncnorm(n = nSpp, mean = max(x_o), sd = 2, a = max(x_o))
  e1 <- e - d
  beta_0 <- rnorm(n = nSpp, mean = 0, sd = 2)
  beta_1 <- rtruncnorm(n = nSpp, mean = 0, sd = 2, b = 0)
  tseq <- seq(0,20, length.out = 1000)
  
  
  #convert data scale to supported kimurswamy scale using d and e
  #produce mean y values
  x <- array(NA, c(ntreat, nSpp))
  mus <- array(NA, dim = c(ntreat,nSpp))
  p_zero <- array(dim = c(ntreat, nSpp))
  for (i in 1:ntreat) {
    x[i,] <- (x_o[i] - d) / e1
    mus[i,] <- stretch.kumara(x[i,], a = a, b = b, c = c)
    p_zero[i, ] <- plogis(beta_0 + beta_1 * mus[i, ])
  }


  #generate gamma distributed variation around mean y values
  #store data
  dimDat <- list(treat = rep(NA, nobs), spp = rep(NA, nobs), y = rep(NA,nobs))
  
  drow <- 1
  for (i in 1:nSpp) {
    for (j in 1:ntreat) {
      for (k in 1:nSppTreat) {
        dimDat$treat[drow] <- x_o[j]
        dimDat$spp[drow] <- i
        is_zero <- rbinom(1, 1, p_zero[j, i])
        if (is_zero) {
          dimDat$y[drow] <- 0
        } else {
          dimDat$y[drow] <- rgamma(n = 1, shape = nu, rate = (1-p_zero[j, i]) * nu / mus[j, i])
        }
        drow <- drow + 1
      }
    }
  #add to ydat matrix
  ydat[,z] <- dimDat$y
  }
}

plot(NA,NA, xlim=range(dimDat$treat), ylim = range(ydat), 
     cex = 0.2,xlab = 'Treatment', ylab = 'Response')

apply(ydat, 2, 
      function(x){
        points(jitter(dimDat$treat), x,
             #col = as.color(dimDat$spp, alpha = 1), 
             pch = 19, cex=0.2)
      })


plot(density(ydat), lwd=2)
lines(density(prior_data$Infl_totwt, na.rm = T), lty=2, lwd=2)

qseq <- seq(0, 1, length.out = 100)
qsim <- quantile(ydat, qseq)
qprior <- quantile(prior_data$Infl_totwt, qseq, na.rm = T)
plot( qsim, qprior, xlim=range(qprior), ylim=range(qprior))
plot( qsim, qprior)
abline(0,1, lty=2)

quantile(ydat)
quantile(prior_data$Infl_totwt, na.rm = T)

