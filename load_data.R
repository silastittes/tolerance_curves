source("tolerance_functions.R")
#LOAD DATA
emery <-load_emery()
stanDat <- load_stanDat()
posts <- rstan::extract(stanDat)
ndraws <- nrow(posts$lp__)
lasth <- load_lasth()
