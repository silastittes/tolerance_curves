################################
################################
### install missing packages ###
################################
################################

pkgs <- c("scales", "dplyr", "parallel", 
          "truncnorm", "geomorph", "ape", 
          "phytools", "OUwie", "rstan", "purrr",
          "TeachingDemos", "gdata"
          )

needed <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]

if(length(needed)) install.packages(needed, dependencies = TRUE)


################
################
### packages ###
################
################

#MARKDOWN AND GRAPHICS
#library(rmarkdown)
#library(knitr)
library(scales)

#tidy
library(dplyr)
library(purrr)
#library(ggplot2)


#other
library(parallel)
#library(smoothmest)
#library(as.color)
#library(formatR)
#library(dtw) 
#library(GPfit)
library(truncnorm)
library(gdata) #read.xls


#phylo comp tools
library(geomorph)
library(ape)
library(phytools)
library(OUwie)

#rstan
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

################################
################################
######### PAR SETTINGS #########
################################
################################

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


load_emery <- function(){
  
  emery <- read.xls("data/Inundation_compiled_FINAL.xlsx") %>%
    mutate(Inflor_biomass = ifelse( is.na(Inflor_biomass) & 
                                      ifEmerge.Y.N. == 1, yes =  0, 
                                    no = Inflor_biomass)) %>%
    mutate(treat = ifelse(Treatment == "F", yes = 5, 
                          no = ifelse(Treatment == "MF", yes = 4, 
                                      no = ifelse(Treatment == "B", yes = 3, 
                                                  no = ifelse(Treatment == "MD", yes = 2, 
                                                              no = 1))))) %>%
    filter(!is.na(Inflor_biomass))
  
  species <- unique(emery$Species)
  sppint <- 1:length(species)
  emery$sppint <- rep(NA, nrow(emery))
  for(i in 1:length(species)){
    emery$sppint[emery$Species %in% species[i]]  <- sppint[i]
  }
  return(emery)
}

#########################################
#########################################
##### FUNCTIONS FOR TOLERANCE CURVE #####
#########################################
#########################################

#input (0,1) vector, output corresponding values from model described in text.
stretch.kumara <- function(x, a, b, c){
  return(c*((a*b*x^(a-1) ) * (1-x^a)^(b-1)))
}

#convert data scale values then return values from model
stretch.kumara2 <- function(xs, a, b, c, d, e1){
  x <- (xs - d)/e1
  mod.fit <- c*((a*b*x^(a-1) ) * (1-x^a)^(b-1))
  return(mod.fit)
}

#for returning rescaled x axis and corresponding output 
plot.kumara <- function(xs, a, b, c, d, e1){
  x <- xs * e1 + d
  mod.fit <- c*((a*b*xs^(a-1) ) * (1-xs^a)^(b-1))
  list(x, mod.fit)
}

#for returning output after rescaling to supported (0,1)
scale.kumara <- function(xs, a, b, c, d, e1){
  x <- (xs - d)/e1
  mod.fit <- c*((a*b*x^(a-1) ) * (1-x^a)^(b-1))
  return(mod.fit)
}

#for returning rescaled x axis and corresponding output
unscale.kumara <- function(x, a, b, c, d, e1){
  xs <- x*e1 + d
  mod.fit <- c*((a*b*x^(a-1) ) * (1-x^a)^(b-1))
  return(list(xs = xs, mod.fit=mod.fit))
}

#for getting appoximate integral
int_kumara <- function(xs, a, b, c, d, e1){
  x <- xs * e1 + d
  dx <- diff(range(x))/length(x)
  mod.fit <- c*((a*b*xs^(a-1) ) * (1-xs^a)^(b-1))
  sum(mod.fit)*dx
}

#for getting approximate optimum
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
load_maxima <- function(){
  if(file.exists("derived_files/maxima_draws.txt")){
    read.table(file = "derived_files/maxima_draws.txt", header = T)
  } else {
    warning("file does not exist yet. Make using generate_derived.R")
  }
}

load_integral <- function(){
  if(file.exists("derived_files/integral_draws.txt")){
    read.table(file = "derived_files/integral_draws.txt", header = T)
  } else{
    warning("file does not exist yet. Make using generate_derived.R")
  }
}

load_stanDat <- function(){
  stan_samples <- list.files("bayes/samples/")[grep ("tolerance_v3.samples",
                                             list.files("bayes/samples/"))]
 rstan::read_stan_csv(paste0("bayes/samples/", stan_samples))
}

#########################################
#########################################
############ LOAD PHYLO DATA ############
#########################################
#########################################
#Load Tree
load_lasth <- function(){
  set.seed(07142015)
  lasth <- read.tree("data/LastheniaBayesian.tre")
  lasth <- root(phy = lasth, outgroup = c("eriophyllum", "amblyopappus"))
  lasth$node.label <- round(as.numeric(lasth$node.label), 2)

  lasth_tips <- unlist(strsplit(lasth$tip.label, split = "L."))
  lasth_tips <- unlist(strsplit(lasth_tips , split = "'"))
  lasth_tips <- lasth_tips[lasth_tips != ""]
  lasth$tip.label <- lasth_tips
  
  holoAdd <- which(lasth$tip.label =="sect.Hologymne")
  newtips <- c("coulteri", "glabrata", "ferrisiae", "chrysantha")
  
  
  addbr <- 0.0001
  addtree <- rtree(n = length(newtips), rooted = T, 
                   tip.label = newtips, 
                   br = rtruncnorm(n = length(newtips), a = addbr, 
                                   mean = addbr, sd = 0.0001)
                   )
  
  lasth <- bind.tree(x = lasth, y = addtree, where = holoAdd)
  
  #drop
  drop <- lasth$tip.label[!lasth$tip.label %in% unique(emery$Species)]
  lasth <- drop.tip(phy = lasth, tip = drop)
  lasth <- chronos(lasth, lambda = 1)
  return(lasth)
}

?ape::chronos

load_lasth()
