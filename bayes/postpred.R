#POSTERIOR PREDICTIVE CHECKS FOR LASTHENIA TOLERANCE CURVE MODEL
source("tolerance_functions.R")
source("load_data.R")

emery <- load_emery() %>% filter(Inflor_biomass != 0)

ndraws <- nrow(posts$lp__)
ydat <- 1:ndraws %>% mclapply(function(z){
  nu <- posts$nu[z]
  a <- posts$a[z,]
  b <- posts$b[z,]
  c <- posts$c[z,]
  d <- posts$d[z,]
  e1 <-posts$e1[z,]
  beta_0 <- posts$beta_0[z,]
  beta_1 <- posts$beta_1[z,]
  
  
  1:nrow(emery) %>% sapply(function(i){
    x <- (emery$treat[i] - d[emery$sppint[i]]) / e1[emery$sppint[i]]
    mus <- stretch.kumara(x, a = a[[emery$sppint[i]]], 
                             b = b[[emery$sppint[i]]], 
                             c = c[[emery$sppint[i]]])
    p_zero <- plogis(beta_0[[emery$sppint[i]]] + 
                          beta_1[[emery$sppint[i]]] * mus)
    
    rgamma(n = 1, shape = nu, 
           rate = (1-p_zero) * nu / mus)
    
  })
})


#ssq_sim <- 1:ndraws %>% lapply(function(z){
ssq_sim <- 1:ndraws %>% sapply(function(z){
  nu <- posts$nu[z]
  a <- posts$a[z,]
  b <- posts$b[z,]
  c <- posts$c[z,]
  d <- posts$d[z,]
  e1 <-posts$e1[z,]
  beta_0 <- posts$beta_0[z,]
  beta_1 <- posts$beta_1[z,]
  
  squares <- 1:nrow(emery) %>% sapply(function(i){
    x <- (emery$treat[i] - d[emery$sppint[i]]) / e1[emery$sppint[i]]
    mus <- stretch.kumara(x, a = a[[emery$sppint[i]]], 
                          b = b[[emery$sppint[i]]], 
                          c = c[[emery$sppint[i]]])
    p_zero <- plogis(beta_0[[emery$sppint[i]]] + 
                       beta_1[[emery$sppint[i]]] * mus)
    
    pred <- rgamma(n = 1, shape = nu, 
                   rate = (1-p_zero) * nu / mus)
    
    rbind((emery$Inflor_biomass[i] - mus)^2, (pred - mus)^2)
  })
  
  apply(squares, 1, sum)
})

plot(ssq_sim[1,], ssq_sim[2,], pch = ".")
abline(0,1, lwd=2, col = "grey")

ggplot(data.frame(t(ssq_sim)), aes(x = X1, y = X2)) +
  geom_hex() +
  geom_abline(slope = 1, intercept = 0, colour = "grey") +
  theme_bw()

#bayesian p-value for sum of squares discrepency
mean(ssq_sim[1,] < ssq_sim[2,])

#create logical vector to remove zeros
trt <- emery$Inflor_biomass != 0
plot(ydat[[1]][trt], emery$Inflor_biomass[trt])
val <- 0.5
quant_rep_y <- sapply(ydat, function(x) quantile(x[trt], probs = val))
hist(quant_rep_y)
abline(v=quantile(emery$Inflor_biomass[trt], val))
mean(quant_rep_y > quantile(emery$Inflor_biomass[trt], val))

min_rep_y <- sapply(ydat, function(x) min(x[trt]))
hist(min_rep_y)
abline(v=min(emery$Inflor_biomass[trt]))
mean(min_rep_y > min(emery$Inflor_biomass[trt]))

max_rep_y <- sapply(ydat, function(x) max(x[trt]))
hist(max_rep_y)
abline(v=max(emery$Inflor_biomass[trt]))
mean(max_rep_y > max(emery$Inflor_biomass[trt]))


Q <- 3
10^(-Q/10)
