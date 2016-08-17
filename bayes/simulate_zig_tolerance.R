library(rstan)
library(truncnorm)

#simulate data set for tolerance model

# Silas heads up - I removed the call to setwd()
setwd("~/Documents/Projects/toleranceCurves/bayes/")

## FUNCTIONS FOR TOLERANCE CURVE -------------------------------
stretch.kumara <- function(x, a, b, c){
  return(c*((a*b*x^(a-1) ) * (1-x^a)^(b-1)))
}

kumara_dydx <- function(x, a, b, c){
  -a * b * c * x^(-2 + a) *
    (1 - x^a)^(-2 + b) *
    (1 - x^a + a * (-1 + b * x^a))
}

HDI <- function(values, percent=0.95){
  sorted <- sort(values)
  index <- floor(percent * length(sorted))
  nCI <- length(sorted) - index
  width <- rep(0, nCI)
  for (i in 1:nCI){
    width[i] <- sorted[i + index] - sorted[i]
    }
  HDImin <- sorted[which.min(width)]
  HDImax <- sorted[which.min(width) + index]
  HDIlim <- c(HDImin, HDImax)
  return(HDIlim)
}

## simulate parameters and data ----------------------------
nSpp <- 2
nobs <- 400
ntreat <- 10
x_o <- seq(1, ntreat, length.out = ntreat) #data scale
nSppTreat <- nobs/nSpp/ntreat

simreps <- 1
ydat <- matrix(NA, nrow = nobs, ncol = simreps)
for (z in 1:simreps) {
  nu <- rgamma(1, 100, 100 / 50)
  a <- rtruncnorm(n = nSpp, mean = 5, sd = 2, a = 2)
  b <- rtruncnorm(n = nSpp, mean = 5, sd = 2, a = 2)
  c <- rtruncnorm(n = nSpp, mean = 1, sd = 2, a = 2)
  d <- rtruncnorm(n = nSpp, mean = min(x_o), sd = 2, b = min(x_o))
  e <- rtruncnorm(n = nSpp, mean = max(x_o), sd = 2, a = max(x_o))
  e1 <- e - d
  beta_0 <- rtruncnorm(1, mean = 0, sd = 1, a = -10)
  beta_1 <- rtruncnorm(1, mean = -2, sd = 0.1, b = 0)
  tseq <- seq(0,20, length.out = 1000)
  plot(tseq, plogis(beta_0 + beta_1*tseq))

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
          dimDat$y[drow] <- rgamma(n = 1, shape = nu, rate = nu / mus[j, i])
        }
        drow <- drow + 1
      }
    }
  }
  #add to ydat matrix
  ydat[,z] <- dimDat$y
}

range(ydat)

fmu <- 0.1
rgamma(n = 1, shape = nu, rate = nu / fmu)
rgamma(n = 1, shape = nu, rate = nu / fmu/0.9)


dimDat$yScale <- by(data = dimDat$y, INDICES = dimDat$spp, FUN = function(x) x/max(x))
dimDat$yScale <- c(dimDat$yScale$`1` , dimDat$yScale$`2`)
by(data = dimDat$y, INDICES = dimDat$spp, FUN = max)

#visualize simulated data
plot(jitter(dimDat$treat), ydat[,1],
     col = dimDat$spp, pch = 19,
     ylim = range(ydat), cex = 0.2,
     xlab = 'Treatment',
     ylab = 'Response')


#visualize simulated data
plot(jitter(dimDat$treat), dimDat$y,
     col = dimDat$spp, pch = 19,
     cex = 0.2,
     xlab = 'Treatment',
     ylab = 'Response')


# bundle data for stan
dimDat <- data.frame(dimDat)
dataList <- list(N = length(dimDat$spp),
               x = dimDat$treat,
               #y = dimDat$yScale, #!!!!!!!!!! SCALED !!!!!!!!
               y = dimDat$y,
               numSpp = length(unique(dimDat$spp)),
               sppint = dimDat$spp)
# estimate parameters
watch <- c("a", "b", "c", "d", "e", "e1", "nu", "beta_0", "beta_1", "mu")
stan.fit.sim <- stan(file = "tolerance_zig_v3_nohier.stan",
                     data = dataList, 
                     #iter = 200, chains = 2,
                pars = watch)

#stan.fit.sim
#rstan::traceplot(stan.fit.sim, pars = watch)


watch2 <- c("a", "b", "c", "d", "e", "e1", "nu")
#stan.fit.sim2 <- stan(file = "~/Documents/GradSchool/courses/hier_reg/project/Max_help/tolerance_gamma.stan",
#                     data = dataList, 
#                     #iter = 200, chains = 2,
#                     pars = watch2)

#stan.fit.sim2
#stan.fit.sim

# explore estimated relationship between mu and p_zero
post <- rstan::extract(stan.fit.sim)
mu_vals <- seq(min(mus), max(mus), length.out = 100)

# add line for every draw
ndraw <- length(post$lp__)

plot(x = NULL, y = NULL, xlim = range(mus), ylim = c(0, 1),
     xlab = "mu", ylab = 'Pr(zero)')
for (i in 1:ndraw) {
  lines(mu_vals, plogis(post$beta_0[i] + post$beta_1[i] * mu_vals))
}
# add true
lines(mu_vals, plogis(beta_0 + beta_1 * mu_vals), col = 'red', lwd = 3)

#plot(jitter(dimDat$treat), ydat,
#     col = dimDat$spp, pch = 19,
#     ylim = range(ydat), cex = 0.2,
#     xlab = 'Treatment',
#     ylab = 'Response')


plot(jitter(dimDat$treat), ydat,
     col = dimDat$spp, pch = 19,
     ylim = c(0,2), xlim=c(5,8),
     cex = 0.2,
     xlab = 'Treatment',
     ylab = 'Response')

#plot(jitter(dimDat$treat), dimDat$yScale,
plot(jitter(dimDat$treat), dimDat$y,
     col = dimDat$spp, pch = 19,
     cex = 0.2,
     type="n",
     xlab = 'Treatment',
     ylab = 'Response', xlim=range(dimDat$treat)*c(1,1.2))

xseq <- seq(0,1, length.out=500)
for(i in 1:ndraw){
  for(j in 1:ncol(post$e1)){
    
    xdat1 <- xseq*post$e1[i,j] + post$d[i,j]
    ydat1 <- stretch.kumara(xseq, post$a[i,j], post$b[i,j], post$c[i,j])
    yzig <- ydat1 * (1 - plogis(post$beta_0[i] + post$beta_1[i] * ydat1))
    
    lines(xdat1, ydat1, col=j)
    #lines(xdat1,yzig, col=j)
    
  }
}

slen <- unique(dimDat$spp)
for(i in slen){
  subDat <- subset(dimDat, spp == i)
  #spln <- smooth.spline(subDat$treat, subDat$yScale)
  spln <- smooth.spline(subDat$treat, subDat$y)
  lines(spln$x, spln$y, col="cyan", lwd=3)
}

points(jitter(dimDat$treat), dimDat$y,
       bg = dimDat$spp, pch = 21)

plotInt <- function(param, d, value){
  
  di <- dim(param)
  
  if( length(di) < 2){
    pDen <- density(param)
    plot(pDen)
    high <- max(pDen$y)*0.2
    qs <- quantile(param, c(0.025, 0.5, 0.975))
    pMean <- mean(param)
  } else {
    pDen <- density(param[,d])
    plot(pDen)
    high <- max(pDen$y)*0.2
    qs <- quantile(param[,d], c(0.025, 0.5, 0.975))
    pMean <- mean(param[,d])
  }
  
  segments(x0 = value[d], y0 = 0, x1 = value[d], y1 = high, col="blue", lwd=4)
  segments(x0 = qs, y0 = 0, x1 = qs, y1 = high, col="red", lwd=4)
  segments(x0 = pMean, y0 = 0, x1 = pMean, y1 = high, col="green", lwd=4)
  
  legend("topleft", c("true", "95% cred int", "mean"), 
         lwd=3, col=c("blue","red", "green"), bty="n", cex = 1)
}

plotInt(param = post$a, d = 1, value = a)
plotInt(param = post$a, d = 2, value = a)

plotInt(param = post$b, d = 1, value = b)
plotInt(param = post$b, d = 2, value = b)

plotInt(param = post$c, d = 1, value = c)
plotInt(param = post$c, d = 2, value = c)

plotInt(param = post$d, d = 1, value = d)
plotInt(param = post$d, d = 2, value = d)

plotInt(param = post$e1, d = 1, value = e1)
plotInt(param = post$e1, d = 2, value = e1)

plotInt(param = post$nu, d = 1, value = nu)

plotInt(param = post$beta_0, d = 1, value = beta_0)
plotInt(param = post$beta_1, d = 1, value = beta_1)


meanDraws <- apply(X = post$mu, 2, mean)
#plot(emery$Inflor_biomass, meanDraws)
#abline(0,1, lty=2, col="red")
plot(meanDraws, ydat - meanDraws,
     pch=as.integer(emery$Species))

#abline(h=0)

smth <- smooth.spline(meanDraws, ydat - meanDraws)
lines(smth$x, smth$y)

