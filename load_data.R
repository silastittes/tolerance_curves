source("tolerance_functions.R")
source("derived_files/lasth_100_post.R")

#LOAD DATA
emery <-load_emery()
stanDat <- load_stanDat()
posts <- rstan::extract(stanDat)
ndraws <- nrow(posts$lp__)
lasth <- load_lasth()
grad <- load_gradient()
draws <- read_csv("bayes/stan_par1_df.csv") #draws for penalized zero model
