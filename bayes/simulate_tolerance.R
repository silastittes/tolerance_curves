#simulate data set for tolerance model

#load functions and data
source("bayes/tolerance_functions.R")

emery <- load_emery()

## simulate parameters and data ----------------------------
nSpp <- 4
nobs <- 400
ntreat <- 10
x_o <- seq(1, ntreat, length.out = ntreat) #data scale
nSppTreat <- nobs/nSpp/ntreat

simreps <- 1
ydat <- matrix(NA, nrow = nobs, ncol = simreps)
for (z in 1:simreps) {
  nu <- rgamma(1, 20, 0.2)
  a <- rtruncnorm(n = nSpp, mean = 4, sd = 1, a = 2)
  b <- rtruncnorm(n = nSpp, mean = 3, sd = 1, a = 2)
  c <- rtruncnorm(n = nSpp, mean = 0, sd = 2, a = 0)
  d <- rtruncnorm(n = nSpp, mean = min(x_o), sd = 2, b = min(x_o))
  e <- rtruncnorm(n = nSpp, mean = max(x_o), sd = 2, a = max(x_o))
  e1 <- e - d
  
  
  #beta_0 <- rnorm(1, mean = 0, sd = 2)
  #beta_1 <- rtruncnorm(1, mean = 0, sd = 2, b = 0)

  
  beta_0 <- rnorm(1, mean = 2, sd = 0.2)
  beta_1 <- rtruncnorm(1, mean = 0, sd = 0.2, b = 0)
  
  
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

  
#approximate zero gamma sampling
  drow <- 1
  for (i in 1:nSpp) {
    for (j in 1:ntreat) {
      for (k in 1:nSppTreat) {
        dimDat$treat[drow] <- x_o[j]
        dimDat$spp[drow] <- i
        dimDat$y[drow] <- rgamma(n = 1, shape = nu, 
                                 rate = (p_zero[j, i]*nu) / mus[j, i])
        drow <- drow + 1
      }
    }
  }
  
  #add to ydat matrix
  ydat[,z] <- dimDat$y

  
}

  
  
#zero inflated sampling 
#  drow <- 1
#  for (i in 1:nSpp) {
#    for (j in 1:ntreat) {
#      for (k in 1:nSppTreat) {
#        dimDat$treat[drow] <- x_o[j]
#        dimDat$spp[drow] <- i
#        is_zero <- rbinom(1, 1, p_zero[j, i])
#        if (is_zero) {
#          dimDat$y[drow] <- 0
#        } else {
#          dimDat$y[drow] <- rgamma(n = 1, shape = nu, rate = nu / mus[j, i])
#        }
#        drow <- drow + 1
#      }
#    }
#  }
  #add to ydat matrix
# ydat[,z] <- dimDat$y
#}

range(ydat)


dimDat$yScale <- by(data = dimDat$y, INDICES = dimDat$spp, FUN = function(x) x/max(x))
dimDat$yScale <- c(dimDat$yScale$`1` , dimDat$yScale$`2`)
by(data = dimDat$y, INDICES = dimDat$spp, FUN = max)

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
stan.fit.sim <- stan(file = "bayes/tolerance_v3.stan",
                     data = dataList, 
                     iter = 200, chains = 4,
                pars = watch)


#stan.fit.sim
#rstan::traceplot(stan.fit.sim, pars = watch)


watch2 <- c("a", "b", "c", "d", "e", "e1", "nu")
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

#plot(jitter(dimDat$treat), dimDat$yScale,
plot(jitter(dimDat$treat), dimDat$y,
     col = dimDat$spp, pch = 19,
     cex = 0.2,
     type="n",
     xlab = 'Treatment',
     ylab = 'Response')#, xlim=range(dimDat$treat)*c(1,1.2))

xseq <- seq(0,1, length.out=500)
for(i in 1:ndraw){
  for(j in 1:ncol(post$e1)){
    xdat1 <- xseq*post$e1[i,j] + post$d[i,j]
    ydat1 <- stretch.kumara(xseq, post$a[i,j], post$b[i,j], post$c[i,j])
    
    lines(xdat1, ydat1, col=j)
    
  }
}


## ADD "TRUE" CURVES 
for(i in 1:nSpp){
  xdat1 <- xseq*e1[i] + d[i]
  ydat1 <- stretch.kumara(xseq, a[i], b[i], c[i])
  lines(xdat1, ydat1, col="cyan", lwd=2)
}


slen <- unique(dimDat$spp)
for(i in slen){
  subDat <- subset(dimDat, spp == i)
  #spln <- smooth.spline(subDat$treat, subDat$yScale)
  spln <- smooth.spline(subDat$treat, subDat$y)
  lines(spln$x, spln$y, col="yellow", lwd=3)
}

points(jitter(dimDat$treat), dimDat$y,
       bg = dimDat$spp, pch = 21)

plotInt <- function(param, d, value){
  
  di <- dim(param)
  par_name <- deparse(substitute(param))
  only_name <- unlist(strsplit(x = par_name, split = "\\$"))[2]
  main_par <- paste(only_name, d)
  
  if( length(di) < 2){
    pDen <- density(param)
    plot(pDen, main = main_par)
    high <- max(pDen$y)*0.2
    qs <- quantile(param, c(0.025, 0.5, 0.975))
    pMean <- mean(param)
  } else {
    pDen <- density(param[,d])
    plot(pDen, main = main_par)
    high <- max(pDen$y)*0.2
    qs <- quantile(param[,d], c(0.025, 0.5, 0.975))
    pMean <- mean(param[,d])
  }
  

  segments(x0 = qs, y0 = 0, x1 = qs, y1 = high, col="red", lwd=4)
  segments(x0 = pMean, y0 = 0, x1 = pMean, y1 = high, col="green", lwd=4)
  segments(x0 = value[d], y0 = 0, x1 = value[d], y1 = high, col="blue", lwd=4)
  
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
plot(meanDraws, ydat - meanDraws,
     pch=as.integer(dimDat$spp))
smth <- smooth.spline(meanDraws, ydat - meanDraws)
lines(smth$x, smth$y)

by(data = dimDat$y, 
   INDICES = list(species = dimDat$spp, trt = dimDat$treat), 
   FUN = function(x) max(x) > 0 )



cl <- 1
mcmc_p <- post$b[, cl]
t_p <- a[cl]
mcmc_q <- quantile(mcmc_p, c(0.025, 0.975))
t_p
mcmc_q[1] < t_p & t_p < mcmc_q[2]
