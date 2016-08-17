setwd("~/Documents/Projects/toleranceCurves/")
source("bayes/tolerance_functions.R")
#LOAD DATA
emery <-load_emery()
maximadf <- load_maximadf()
stanDat <- load_stanDat()
posts <- extract(stanDat)
ndraws <- nrow(posts$lp__)
summs <- summary(stanDat)$summary
lasth <- load_lasth()