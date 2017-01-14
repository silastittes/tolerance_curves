#simulate data set for tolerance model

#load functions and data comment once loaded to allow seed setting here

source("tolerance_functions.R")
emery <- load_emery()
mod <- rstan::stan_model("bayes/tolerance_v3.stan")


set.seed(890567703)

nsims <- 100
#storage for credible interval
cred_mat <- array(NA, dim = c(4, 8, nsims))
prop_zero <- array(NA, dim = c(4, nsims))


for(sim in 1:nsims){
## simulate parameters and data ----------------------------
#nSppTreat <- nobs/nSpp/ntreat
nSppTreat <- 10
ntreat <- 5
nSpp <- 4
#nobs <- 50*nSpp
nobs <- nSpp*ntreat*nSppTreat
x_o <- seq(1, ntreat, length.out = ntreat) #data scale

simreps <- 1

ydat <- matrix(NA, nrow = nobs, ncol = simreps)
for (z in 1:simreps) {
  nu <- rep(rgamma(1, 20, 0.2), nSpp)
  a <- rtruncnorm(n = nSpp, mean = 4, sd = 1, a = 2)
  b <- rtruncnorm(n = nSpp, mean = 3, sd = 1, a = 2)
  c_hyper <- rtruncnorm(n = 1, mean = 0, sd = 2, a = 0)
  #c <- rtruncnorm(n = nSpp, mean = c_hyper, sd = 1, a = 0)
  c <- rep(c_hyper, nSpp) #for demonstration
  d <- rtruncnorm(n = nSpp, mean = min(x_o), sd = 1, b = min(x_o))
  e <- rtruncnorm(n = nSpp, mean = max(x_o), sd = 1, a = max(x_o))
  
  d[c(3, 4)] <- rtruncnorm(n = 1, mean = min(x_o)-4, sd = 2, b = min(x_o))
  e[c(3, 4)] <- rtruncnorm(n = 1, mean = max(x_o)+7, sd = 2, a = max(x_o))
  e1 <- e - d


  beta_0 <- rnorm(nSpp, mean = -5, sd = 0.2)
  beta_1 <- rtruncnorm(nSpp, mean = -0.1, sd = 0.2, b = 0)
  
  #for illustrative purposes
  #c[2] <- c[2]/2 #ensure plently of zeros for species 2
  beta_0 <- rep(-5, nSpp)
  beta_1 <- rep(-1, nSpp)
  beta_0[c(2,4)] <- 1
  beta_1[c(2,4)] <- -0.2
  
  #tseq <- seq(0,20, length.out = 1000)
  #plot(tseq, plogis(beta_0[2] + beta_1[2]*tseq))

  #convert data scale   to supported kimurswamy scale using d and e
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
        
        #dimDat$y[drow] <- rgamma(n = 1, shape = nu, 
        #                         rate = (p_zero[j, i]*nu) / mus[j, i])
        
        dimDat$y[drow] <- rgamma(n = 1, shape = nu, 
                                 rate = (nu) / mus[j, i])
        
        if(rbernoulli(n = 1, p = p_zero[j, i])){
          dimDat$y[drow] <- 0
        }
        drow <- drow + 1
      }
    }
  }
  
  #add to ydat matrix
  ydat[,z] <- dimDat$y

}

#dimDat$yScale <- by(data = dimDat$y, INDICES = dimDat$spp, FUN = function(x) x/max(x))
#dimDat$yScale <- c(dimDat$yScale$`1` , dimDat$yScale$`2`)
#by(data = dimDat$y, INDICES = dimDat$spp, FUN = max)

# bundle data for stan
dimDat <- data.frame(dimDat)

#colset <- c("blue", "black", "green", "red")
#visualize simulated data
#plot(jitter(dimDat$treat), dimDat$y,
#     #col = dimDat$spp, pch = 19,
#     col = colset, pch = 19,
#     cex = 0.5,
#     xlab = 'Treatment',
#     ylab = 'Response')

#smoothdata <- dimDat %>% group_by(spp) %>%
#  do(smooth = smooth.spline(.$treat, .$y))

#invisible(
#  lapply(smoothdata$spp, 
#         function(i) lines(smoothdata$smooth[[i]]$x, 
#                           smoothdata$smooth[[i]]$y, col = colset[i]
#         )
#  )
#)


dataList <- list(N = length(dimDat$spp),
               x = dimDat$treat,
               #y = dimDat$yScale, #!!!!!!!!!! SCALED !!!!!!!!
               y = dimDat$y,
               numSpp = length(unique(dimDat$spp)),
               sppint = dimDat$spp)


# estimate parameters
watch <- c("a", "b", "c", "d", "e1", "nu", "beta_0", "beta_1", "mu")
stan.fit.sim <- sampling(mod,
                     data = dataList, 
                     iter = 500, chains = 4,
                     pars = watch)

post <- rstan::extract(stan.fit.sim)
mu_vals <- seq(min(mus), max(mus), length.out = 100)


cred_match <- sapply(as.list(watch[-9]), function(p){
  sapply(as.list(1:nSpp), function(s){
    param <- eval(parse(text = p))
    
    dimtest <- dim(post[[p]])
    if( length(dimtest) > 1 ){
      if(dimtest[2] < s){
        NULL
      } else {
      qn <- quantile(post[[p]][,s], c(0.025, 0.975))
      param[s] > qn[1] & param[s] < qn[2]  
      } 
    } else {
      qn <- quantile(post[[p]], c(0.025, 0.975))
      param[s] > qn[1] & param[s] < qn[2]  
    }
      })
  })


  cred_mat[1:4,1:8,sim] <- cred_match
  prop_zero[,sim] <- by(dimDat$y, dimDat$spp, function(x) mean(x ==0))
  print(paste("loop", sim))
} #sim loop starts at very top of file


dump("cred_mat", file = "bayes/simulate100.R")
dump("prop_zero", file = "bayes/simulate100_prop_zero.R")

# add line for every draw
ndraw <- length(post$lp__)

#plot(x = NULL, y = NULL, xlim = range(mus), ylim = c(0, 1),
#     xlab = "mu", ylab = 'Pr(zero)')
#for (i in 1:ndraw) {
#  lines(mu_vals, plogis(post$beta_0[i] + post$beta_1[i] * mu_vals))
#}
# add true
#lines(mu_vals, plogis(beta_0 + beta_1 * mu_vals), col = 'red', lwd = 3)

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
    
    lines(xdat1, ydat1, col=alpha(j, 0.05))
    
  }
}

## ADD "TRUE" CURVES 
for(i in 1:nSpp){
  xdat1 <- xseq*e1[i] + d[i]
  ydat1 <- stretch.kumara(xseq, a[i],
                          b[i], c[i]) #!!!
  lines(xdat1, ydat1, col="gold", lwd=2)
}

slen <- unique(dimDat$spp)
for(i in slen){
  subDat <- subset(dimDat, spp == i)
  #spln <- smooth.spline(subDat$treat, subDat$yScale)
  spln <- smooth.spline(subDat$treat, subDat$y)
  lines(spln$x, spln$y, col = alpha("yellow", 0.8), lwd=3, lty = 3)
}

points(jitter(dimDat$treat), dimDat$y,
       col = dimDat$spp, pch = 21, bg = "grey")

plotInt <- function(param, d, value, ...){
  
  di <- dim(param)
  #par_name <- deparse(substitute(param))
  #only_name <- unlist(strsplit(x = par_name, split = "\\$"))[2]
  #main_par <- paste(only_name, d)
  
  if( length(di) < 2){
    pDen <- density(param)
    plot(pDen, ...)
    high <- max(pDen$y)*0.2
    qs <- quantile(param, c(0.025, 0.5, 0.975))
    pMean <- mean(param)
  } else {
    pDen <- density(param[,d])
    plot(pDen, ...)
    high <- max(pDen$y)*0.2
    qs <- quantile(param[,d], c(0.025, 0.5, 0.975))
    pMean <- mean(param[,d])
  }
  

  segments(x0 = qs, y0 = 0, x1 = qs, y1 = high, col="red", lwd=4)
  segments(x0 = pMean, y0 = 0, x1 = pMean, y1 = high, col="green", lwd=4)
  
  #mean(t_p / (p_zero[,cl]))
  segments(x0 = value[d],
           y0 = 0, 
           x1 = value[d],
           y1 = high, col="blue", lwd=4)
  
  legend("topleft", c("true", "95% cred int", "mean"), 
         lwd=3, col=c("blue","red", "green"), bty="n", cex = 1)
}


#plot all post densitites with cred int and true vals
yy <- as.list(watch[-9])
invisible(lapply(as.list(yy), function(y){
  cPost <- eval(parse(text=paste("post$", y , sep = "")))
  lapply(as.list(1:nSpp), function(x){
    plotInt(param = cPost, d = x, value = eval(parse(text = y)), main = paste(y, x))
  })
}))


meanDraws <- apply(X = post$mu, 2, mean)
plot(meanDraws, ydat - meanDraws,
     pch=as.integer(dimDat$spp))
smth <- smooth.spline(meanDraws, ydat - meanDraws)
lines(smth$x, smth$y)

by(data = dimDat$y, 
   INDICES = list(species = dimDat$spp, trt = dimDat$treat), 
   FUN = function(x) max(x) > 0 )


range(dimDat$y)
cl <- 2
mcmc_p <- post$c[,cl]
pp <- plogis(mean(post$beta_0[,cl]) + mean(post$beta_1[,cl])* median(dimDat$treat))
t_p <- c[cl]
quantile(mcmc_p, c(0.0225, 0.975))
t_p * (1-pp)


source("bayes/simulate100.R")
source("bayes/simulate100_prop_zero.R")
watch[-9]


1:8 %>% sapply( function(x) mean(cred_mat[1,x,]))
1:8 %>% sapply( function(x) mean(cred_mat[2,x,]))
1:8 %>% sapply( function(x) mean(cred_mat[3,x,]))
1:8 %>% sapply( function(x) mean(cred_mat[4,x,]))
apply(prop_zero, 1, mean)
