source("tolerance_functions.R")
#LOAD DATA
emery <-load_emery()
stanDat <- load_stanDat()
posts <- extract(stanDat)
ndraws <- nrow(posts$lp__)
summs <- summary(stanDat)$summary
lasth <- load_lasth()
integraldf <- load_integral()
maximadf <- load_maxima()
