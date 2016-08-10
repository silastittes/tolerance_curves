setwd("~/Documents/Projects/toleranceCurves/")

################
################
### packages ###
################
################

#MARKDOWN AND GRAPHICS
library(rmarkdown)
library(knitr)
library(scales)

#Hadley
library(dplyr)
library(ggplot2)
library(purrr)

library(smoothmest)
library(as.color)
library(formatR)
library(dtw) #?
library(GPfit) #?
library(truncnorm)

#phylo comp tools
library(ape)
library(picante)
library(phytools)
library(OUwie)
library(bayou)
library(geiger)

#rstan
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

################################
################################
######### PAR SETTINGS #########
################################
################################

options(scipen=99) #avoid scientific notation when printing
opts_chunk$set(tidy=TRUE, tidy.opts=list(width.cutoff=70, message=FALSE, results='hide'))
options(bitmapType="cairo")

opar <- par()
opar$cin <- NULL
opar$cra <- NULL
opar$csi <- NULL
opar$cxy <- NULL
opar$din <- NULL
opar$page <- NULL



####################
####################
#### TRAIT DATA ####
####################
####################


emery <- read.csv("bayes/NEmeryData.csv", header=T)
emery <- emery[is.na(emery$Inflor_biomass) == 0,]

###DELETE THIS###
#emery <- emery[emery$Inflor_biomass != 0,]

emery$treat <- rep(NA, nrow(emery))

#make new column with treatments quant instead of nominal
for(i in 1:nrow(emery)){
  if( emery$Treatment[i] == "F"){
    emery$treat[i]  <- 5
  } else if(emery$Treatment[i] == "MF"){
    emery$treat[i]  <- 4
  } else if(emery$Treatment[i] == "B"){
    emery$treat[i]  <- 3
  } else if (emery$Treatment[i] == "MD"){
    emery$treat[i]  <- 2
  } else if (emery$Treatment[i] == "D"){
    emery$treat[i]  <- 1
  } else {
    emery$treat[i]  <- NA
  }
}

#make another new column with species as integers instead of nominal for stan
species <- unique(emery$Species)
sppint <- 1:length(species)
emery$sppint <- rep(NA, nrow(emery))
for(i in 1:length(species)){
  emery$sppint[emery$Species %in% species[i]]  <- sppint[i]
}


#emery <- arrange(emery, sppint)
#emery <- arrange(emery, treat)

#by(data = emery$treat, INDICES = emery$Species, FUN = function(x) sum(x == min(x)))
#emery$Inflor_biomass[emery$Species == "minor"]
#unique(names(biomass_max))


#########################################
#########################################
##### FUNCTIONS FOR TOLERANCE CURVE #####
#########################################
#########################################

stretch.kumara <- function(x, a, b, c){
  return(c*((a*b*x^(a-1) ) * (1-x^a)^(b-1)))
}

stretch.kumara2 <- function(xs, a, b, c, d, e1){
  x <- (xs - d)/e1
  mod.fit <- c*((a*b*x^(a-1) ) * (1-x^a)^(b-1))
  return(mod.fit)
}

plot.kumara <- function(xs, a, b, c, d, e1){
  x <- xs * e1 + d
  mod.fit <- c*((a*b*xs^(a-1) ) * (1-xs^a)^(b-1))
  list(x, mod.fit)
}

scale.kumara <- function(xs, a, b, c, d, e1){
  x <- (xs - d)/e1
  mod.fit <- c*((a*b*x^(a-1) ) * (1-x^a)^(b-1))
  return(mod.fit)
}

kumara_dydx <- function(x, a, b, c){
  -a * b * c * x^(-2 + a) *
    (1 - x^a)^(-2 + b) *
    (1 - x^a + a * (-1 + b * x^a))
}


opt.kumara <- function(xs, a, b, c, d, e1){
  x <- xs * e1 + d
  mod.fit <- c*((a*b*xs^(a-1) ) * (1-xs^a)^(b-1))
  x[which.max(mod.fit)]
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

#################################
#################################
######## LOAD HMCMC DATA ########
#################################
#################################


#load posterior draws of all 5 parameters plus the derived maxima parameter
maximadf <- read.table(file = "bayes/maxima_draws.txt", header = T)

stan_samples <- list.files("bayes/")[grep ("tolerance_zig_v3_nohier.samples",
                                              list.files("bayes/"))]

stanDat <- rstan::read_stan_csv(paste0("bayes/",stan_samples))

posts <- extract(stanDat)
ndraws <- nrow(posts$lp__) #*0.1
#summs <- rstan::summary(stanMod)$summary
posts <- extract(stanDat)
summs <- summary(stanDat)$summary

#warnings()

#alternate (normal) model
#stan_samples_Normal <- list.files("bayes/")[grep ("tolerance_zig_vNormal_nohier.samples",
#                                           list.files("bayes/"))]

#stanDat_Normal <- rstan::read_stan_csv(paste0("bayes/",stan_samples_Normal))

#posts_Normal <- extract(stanDat)
#ndraws_Normal <- nrow(posts_Normal$lp__) #*0.1
#summs <- rstan::summary(stanMod)$summary
#posts_Normal <- extract(stanDat_Normal)
#summs <- summary(stanDat_Normal)$summary


#########################################
#########################################
############ LOAD PHYLO DATA ############
#########################################
#########################################
#Load Emery Tree
#lasth1 <- read.tree("LastheniaBayesian.tre")
lasth <- read.tree("ASR/LastheniaBayesian.tre")
lasth <- root(phy = lasth, outgroup = c("eriophyllum", "amblyopappus"))
lasth$node.label <- round(as.numeric(lasth$node.label), 2)
plot.phylo(lasth, show.node.label = T)

lasth_tips <- unlist(strsplit(lasth$tip.label, split = "L."))
lasth_tips <- unlist(strsplit(lasth_tips , split = "'"))
lasth_tips <- lasth_tips[lasth_tips != ""]
lasth$tip.label <- lasth_tips

holoAdd <- which(lasth$tip.label =="sect.Hologymne")
newtips <- c("coulteri", "glabrata", "ferrisiae", "chrysantha")

addbr <- 0.0001
addtree <- rtree(n = length(newtips), rooted = T, 
                 tip.label = newtips, 
                 br = rtruncnorm(n = length(newtips), a = addbr, mean = addbr, sd = 0.0001)
                 )

lasth <- bind.tree(x = lasth, y = addtree, where = holoAdd)


#add polytomies

#drop
drop <- lasth$tip.label[!lasth$tip.label %in% unique(emery$Species)]
lasth <- drop.tip(phy = lasth, tip = drop)
lasth <- chronopl(lasth, lambda = 1)

#######################
#######################
##### MAXIMA DATA #####
#######################
#######################

maximadf <- read.table(file = "bayes/maxima_draws.txt", header = T)
