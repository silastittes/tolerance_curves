#prior predictive checks

library(rstan)
library(truncnorm)

#simulate data set for tolerance model
setwd("~/Documents/Projects/toleranceCurves/bayes/")
source("tolerance_functions.R")
emery <- load_emery()

## simulate parameters and data ----------------------------
nSpp <- length(unique(emery$Species))
ntreat <- length(unique(emery$treat))
nobs <- nSpp * ntreat * 7
x_o <- seq(1, ntreat, length.out = ntreat) #data scale
nSppTreat <- nobs/nSpp/ntreat

simreps <- 20
ydat <- matrix(NA, nrow = nobs, ncol = simreps)
for (z in 1:simreps){
  
  nu <- rgamma(1, 10, 0.2)
  a <- rtruncnorm(n = nSpp, mean = 2, sd = 2, a = 1)
  b <- rtruncnorm(n = nSpp, mean = 3, sd = 2, a = 1)
  c <- rtruncnorm(n = nSpp, mean = 2, sd = 10, a = 0)
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
             col = as.color(dimDat$spp, alpha = 1), pch = 19, cex=0.2)
      })


quantile(emery$Inflor_biomass)
quantile(ydat)
hist(ydat[,10])
hist(emery$Inflor_biomass)
