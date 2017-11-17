################################
################################
### install missing packages ###
################################
################################

pkgs <- c("scales", "tidyverse", "parallel", 
          "magrittr", "truncnorm", "geomorph", 
          "ape", "phytools", "rstan", 
          "gdata", "xtable", "abind",
          "glmnet", "nlme", "ggrepel",
          "ggjoy", "rlang", "stringr", "rethinking"
          )

needed <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]

if(length(needed)) install.packages(needed, dependencies = TRUE)


################
################
### packages ###
################
################

#MARKDOWN AND GRAPHICS
library(rmarkdown)
library(knitr)
library(scales)

#tidy
library(tidyverse)
library(stringr)
library(magrittr)
#library(rlang)
#library(ggjoy)
#library(ggrepel)
#library(GGally)


#library(purrr)
#library(ggplot2)


#other
library(truncnorm)
library(gdata) #read.xls
library(xtable)
library(abind)
library(glmnet)


#phylo comp tools
library(ape)
library(geomorph)
library(phytools)
library(nlme)

#rstan and bayesian stuff
library(rethinking)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
extract <- rstan::extract

#overwrite to ensure right function namespace
summarise <- dplyr::summarise
map <- purrr::map

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
  holo <- c("ferrisiae", "glabrata", "coulteri", "chrysantha")
  emery <- read.xls("data/Inundation_compiled_FINAL.xlsx") %>%
    mutate(Inflor_biomass = ifelse(
      is.na(Inflor_biomass) &
        ifEmerge.Y.N. == 1, yes =  0,
      no = Inflor_biomass)) %>%
    mutate(
      treat = ifelse(
        Treatment == "F", yes = 5, 
        no = ifelse(
          Treatment == "MF", yes = 4, 
          no = ifelse(
            Treatment == "B", yes = 3, 
            no = ifelse(
              Treatment == "MD", yes = 2, 
              no = 1
            )
          )
        )
      )
    ) %>%
    filter(!is.na(Inflor_biomass)) %>%
    group_by(Species) %>%
    mutate(sppint = as.integer(Species)) %>% 
    group_by(Species, treat) %>%
    filter(sum(Inflor_biomass > 0) > 0) %>%
    ungroup() %>%
    mutate(
      Species_h = ifelse(
        as.character(Species) %in% holo, 
        "hologymne", as.character(Species)
      ), 
      sppint_h = as.integer(as.factor(Species_h))
    ) 
  
  return(emery)
}

#########################################
#########################################
##### FUNCTIONS FOR TOLERANCE CURVE #####
#########################################
#########################################

#input (0,1) vector, output corresponding values from model described in text.
stretch_kumara <- function(x, a, b, c){
  return(c*((a*b*x^(a-1) ) * (1-x^a)^(b-1)))
}

#convert data scale values then return values from model
stretch_kumara2 <- function(xs, a, b, c, d, e){
  x <- (xs - d)/(e-d)
  mod.fit <- c*((a*b*x^(a-1) ) * (1-x^a)^(b-1))
  return(mod.fit)
}

#for returning rescaled x axis and corresponding output 
plot_kumara <- function(xs, a, b, c, d, e){
  x <- xs * (e - d) + d
  mod.fit <- c*((a*b*xs^(a-1) ) * (1-xs^a)^(b-1))
  list(x, mod.fit)
}

#for returning output after rescaling to supported (0,1)
scale_kumara <- function(xs, a, b, c, d, e){
  x <- (xs - d)/(e-d)
  mod.fit <- c*((a*b*x^(a-1) ) * (1-x^a)^(b-1))
  return(mod.fit)
}

#for returning rescaled x axis and corresponding output
unscale_kumara <- function(x, a, b, c, d, e){
  x <- xs * (e - d) + d
  mod.fit <- c*((a*b*x^(a-1) ) * (1-x^a)^(b-1))
  return(list(xs = xs, mod.fit=mod.fit))
}


#simulate data set for tolerance model
map_kumara <- function(xs, par_df){
  #generate model fits where zeros affect parameters directly
  1:nrow(par_df) %>% map_df(~{
    draw_x <- par_df$draw[.x]
    spp <- par_df$Species[.x]
    d <- par_df$d[.x]
    e <- par_df$e[.x]
    a <- par_df$a[.x]
    b <- par_df$b[.x]
    c <- par_df$c[.x]
    x <- xs*(e - d) + d
    mod_fit <- c*((a*b*xs^(a-1)) * (1-xs^a)^(b-1))
    data.frame(Species = spp, x = x, y = mod_fit, draw = draw_x)
  })
}


map_kumara2 <- function(xs, par_df){
  #generate model fits where zeros DO NOT affect parameters directly
  1:nrow(par_df) %>% map_df(~{
    draw_x <- par_df$draw[.x]
    spp <- par_df$Species[.x]
    d <- par_df$d[.x]
    e <- par_df$e[.x]
    a <- par_df$a[.x]
    b <- par_df$b[.x]
    c <- par_df$c[.x]
    b0 <- par_df$beta_0[.x]
    b1 <- par_df$beta_1[.x]
    x <- xs*(e - d) + d
    mod_fit <- c*((a*b*xs^(a-1)) * (1-xs^a)^(b-1))
    mod_fit <- (1 - plogis(b0 + b1*mod_fit))*mod_fit
    data.frame(Species = spp, x = x, y = mod_fit, draw = draw_x)
  })
}


#simulate data set for tolerance model
revmap_kumara <- function(x, par_df){
  #generate model fits where zeros affect parameters directly
  1:nrow(par_df) %>% map_df(~{
    draw_x <- par_df$draw[.x]
    spp <- par_df$Species[.x]
    d <- par_df$d[.x]
    e <- par_df$e[.x]
    a <- par_df$a[.x]
    b <- par_df$b[.x]
    c <- par_df$c[.x]
    xs <- (xseq_comm - d)/(e - d)
    mod_fit <- c*((a*b*xs^(a-1)) * (1-xs^a)^(b-1))
    data.frame(Species = spp, x = xseq_comm, y = mod_fit, draw = draw_x)
  })
}





#From Max Joseph's course
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


load_stanDat <- function(){
  stan_samples <- list.files("bayes/samples/")[grep ("tolerance_v3_zig_0.samples",
                                             list.files("bayes/samples/"))]
 rstan::read_stan_csv(paste0("bayes/samples/", stan_samples))
}

#########################################
#########################################
############ LOAD PHYLO DATA ############
#########################################
#########################################
#Load Tree
load_lasth <- function(addbr = 0.001){
  #set.seed(seed)
  lasth <- read.tree("data/LastheniaBayesian.tre")
  lasth <- root(phy = lasth, outgroup = c("eriophyllum", "amblyopappus"))
  lasth$node.label <- round(as.numeric(lasth$node.label), 2)

  lasth_tips <- unlist(strsplit(lasth$tip.label, split = "L."))
  lasth_tips <- unlist(strsplit(lasth_tips , split = "'"))
  lasth_tips <- lasth_tips[lasth_tips != ""]
  lasth$tip.label <- lasth_tips
  
  holoAdd <- which(lasth$tip.label =="sect.Hologymne")
  newtips <- c("coulteri", "glabrata", "ferrisiae", "chrysantha")
  
  addtree <- rtree(
    n = length(newtips), rooted = T, 
    tip.label = newtips, 
    #br = rexp(
    #  n = length(newtips), 
    #  rate = 0)
    
    br = rep(0.000001, length(newtips))
    )

  
  lasth <- bind.tree(x = lasth, y = addtree, where = holoAdd)
  
  #drop
  drop <- lasth$tip.label[!lasth$tip.label %in% unique(emery$Species)]
  lasth <- drop.tip(phy = lasth, tip = drop)
  lasth <- di2multi(lasth, 0.000001)
  lasth <- chronos(lasth, lambda = 1)
  return(lasth)
}
