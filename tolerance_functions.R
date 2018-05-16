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
          "ggjoy", "rlang", "stringr", "rethinking", 
          "patchwork"
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
library(patchwork)
#library(rlang)
#library(ggjoy)
#library(GGally)


#library(purrr)
#library(ggplot2)


#other
library(truncnorm)
library(gdata) #read.xls
library(xtable)
library(abind)
library(glmnet)
library(ggrepel)


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

#the proper version
tolerance_mu <- function(xs, a, b, c, d, e){
  x <- (xs - d)/(e-d)
  x %>% map_dbl(~ {
    if(.x > 0 & .x < 1){
      c*((a*b*.x^(a-1) ) * (1-.x^a)^(b-1))
    } else{
      0
    }
  })
}


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



#HIGHEST POSTERIOR DENSITY INTERVAL

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




##############################################
##############################################
### CREATE HABITAT DATA FOR LASTHENIA TAXA ###
##############################################
##############################################

#THIS ISN'T IDEAL, BUT INFORMATION WAS PUT TOGETHER FROM PREVIOUS PUBLISHED PAPERS,
#SO HARD TO THINK OF AN ALTERNATIVE BESIDES HARD-CODING IT.

Genus_species <- c("debilis", "ferrisiae", "chrysantha", 
                   "fremontii", "coulteri", "microglossa",
                   "platycarpha", "conjugens", "gracilis", 
                   "minor", "glabrata", "burkei",
                   "californica", "glaberrima")

state_reg <- c("terrestrial", "vernal", "vernal",
               "vernal", "vernal", "terrestrial",
               "aqua_terr", "vernal", "aqua_terr",
               "terrestrial", "vernal", "vernal",
               "aqua_terr", "vernal")
names(state_reg) <- Genus_species

state_reg_aqua_terr2terr <- c("terrestrial", "vernal", "vernal",
                              "vernal", "vernal", "terrestrial",
                              "terrestrial", "vernal", "terrestrial",
                              "terrestrial", "vernal", "vernal",
                              "terrestrial", "vernal")
names(state_reg_aqua_terr2terr) <- Genus_species

state_reg_aqua_terr2vernal <- c("terrestrial", "vernal", "vernal",
                                "vernal", "vernal", "terrestrial",
                                "vernal", "vernal", "vernal",
                                "terrestrial", "vernal", "vernal",
                                "vernal", "vernal")
names(state_reg_aqua_terr2vernal) <- Genus_species


reg_df <- tibble(
  habit = state_reg, 
  aqua_terr2terr = state_reg_aqua_terr2terr, 
  aqua_terr2vernal = state_reg_aqua_terr2vernal, 
  Species = names(state_reg)
) 







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
    
    br = rep(1e-6, length(newtips))
    )

  
  lasth <- bind.tree(x = lasth, y = addtree, where = holoAdd)
  
  #drop
  drop <- lasth$tip.label[!lasth$tip.label %in% unique(emery$Species)]
  lasth <- drop.tip(phy = lasth, tip = drop)
  lasth <- di2multi(lasth, 1e-5)
  lasth <- chronos(lasth, lambda = 1)
  lasth <- ladderize(lasth, right = T)
  return(lasth)
}

#plot(load_lasth())


################################
################################
#### RUN STAN ON INPUT DATA ####
################################
################################

#MANY OPTIONS ARE HARD CODED, BUT CAN AND SHOULD BE CHANGE FOR NEW DATA SETS

gen_tolerance <- function(df, response = Inflor_biomass, 
                          treatment = treat, group_ids = sppint, 
                          file_id = NULL, rseed = 2323,
                          treed = 10,
                          iters = 5000,
                            thins = 10){
  
  response <- enquo(response)
  treatment <- enquo(treatment)
  group_ids <- enquo(group_ids)
  
  mod <- rstan::stan_model("bayes/tolerance_v3_alt.stan")
  
  min_max_df <- df %>%
    filter(!! response > 0) %>%
    group_by( !! group_ids) %>%
    summarise(minx = min(!! treatment),
              maxx = max(!! treatment))
  
  sel_pull <- function(column){
    df %>% select(!! column) %>% pull()
  }
  
  stan_out <- c(0,1) %>% map(~{
    stan_in <- list(
      N = sel_pull(response) %>% length(), 
      y = sel_pull(response), 
      x = sel_pull(treatment),
      minx = c(min_max_df$minx),
      maxx = c(min_max_df$maxx),
      numSpp = sel_pull(group_ids) %>% unique() %>% length(),
      sppint = sel_pull(group_ids),
      zig = .x,
      # a_pr_mu = -1, a_pr_sig = 0.2,
      # b_pr_mu = -1, b_pr_sig = 0.2,
      # c_pr_mu = 0, c_pr_sig = 0.4,
      # d_pr_mu = 0, d_pr_sig = 0.4,
      # e_pr_mu = 0, e_pr_sig = 0.4
      
      #!!!
      
      a_pr_mu = -1, a_pr_sig = .5,
      b_pr_mu = -1, b_pr_sig = .5,
      c_pr_mu = 0, c_pr_sig = .5,
      d_pr_mu = 0, d_pr_sig = .5,
      e_pr_mu = 0, e_pr_sig = .5
    )
    
    adelt <- 0.8
    
    chain <- 4
    
    
    if(is.null(file_id)){
      sampling(mod,
               data = stan_in,
               control = list(adapt_delta = adelt, max_treedepth = treed), 
               iter = iters, chains = chain, thin = thins,
               seed = rseed
      )
    } else {
      sampling(mod,
               data = stan_in,
               control = list(adapt_delta = adelt, max_treedepth = treed), 
               iter = iters, chains = chain, thin = thins,
               sample_file = paste0("bayes/samples/tolerance_v3_zig_", .x, eval(substitute(file_id)), ".samples"),
              seed = rseed
              )
      
    }
    
    
  })
  
  stan_out
}




gen_stan_df_2 <- function(stan_out, species_order){
  
  #stan_out output from gen_tolerance()
  #species_order vector of desired labels matching group_ids used in gen_tolerance()
  
  new_post <- stan_out %>% extract()
  
  params <- c("d", "e", "a", "b", "c", "nu")
  
  par_df <- params %>% map( ~{
    new_post[[.x]] %>%
      data.frame() %>%
      set_colnames(species_order) %>%
      gather("Species", x) %>%
      set_colnames(c("Species", .x)) %>%
      group_by(Species) %>%
      mutate(draw = 1:n())
  }) %>% 
    do.call(cbind, .) %>%
    select(-starts_with("Species")) %>%
    select(-matches("draw[0-9]")) %>%
    mutate(
      maxima = (((a - 1)/(a*b - 1))^(1/a) * (e - d) + d),
      breadth = (e-d),
      area = c * breadth,
      special = c/breadth
    )
  
  par_df
}  




sim_tolerance_data <- function(n_spp = 3){
  
  
  #priors
  mu_sig <- 1
  
  a_pr_mu <- 4
  a_pr_sig <- mu_sig
  
  b_pr_mu <- 7
  b_pr_sig <- mu_sig
  
  c_pr_mu <- 0
  c_pr_sig <- mu_sig
  
  d_pr_mu <- -5
  d_pr_sig <- mu_sig
  
  e_pr_mu <- 5
  e_pr_sig <- mu_sig
  
  #generate parameters
  pr_sig <- 1
  mu_a <- rnorm(1, a_pr_mu, pr_sig)
  a <- truncnorm::rtruncnorm(n_spp, a = 2, mean = mu_a, sd = a_pr_sig)
  
  mu_b <- rnorm(1, b_pr_mu, pr_sig)
  b <- truncnorm::rtruncnorm(n_spp, a = 2, mean = mu_b, sd = b_pr_sig)
  
  mu_c <- rnorm(1, c_pr_mu, pr_sig)
  c <- truncnorm::rtruncnorm(n_spp, a = 0, mean = mu_c, sd = c_pr_sig)
  
  mu_d <- rnorm(1, d_pr_mu, pr_sig)
  d <- rnorm(n_spp, mu_d, d_pr_sig)
  
  mu_e <- rnorm(1, e_pr_mu, pr_sig)
  e <- rnorm(n_spp, mu_e, e_pr_sig)
  
  ed_test <- 1:n_spp %>% map_lgl(~d[.x] > e[.x]) %>% sum()
  if(ed_test > 0){
    stop("parameter d is greater than parameter e")
  }
  
  mu_nu <- truncnorm::rtruncnorm(1, a = 0, mean = 0, sd = 1)
  nu <- rgamma(n_spp, shape = 10, scale = mu_nu)
  #nu <- truncnorm::rtruncnorm(n_spp, a = 0, mean = mu_nu, sd = 1)
  
  
  #assume same nu for each group
  # nu <- rgamma(1, shape = 10, scale = 1) %>%
  #   rep(n_spp)
  # 
  list(a = a, b = b, c = c, d = d, e = e, nu = nu, 
       mu_a = mu_a, mu_b = mu_b, mu_c = mu_c, 
       mu_d = mu_d, mu_e = mu_e, mu_nu = mu_nu
  )
}



#create predictions from parameters input -- assumes parameter are scalars
gen_tolerance_data <- function(xs, a, b, c, d, e, nu, species_id){
  mu <- tolerance_mu(xs, a, b, c, d, e)
  zero_idx <- xs < d | xs > e
  mu_spp <- mu %>% 
    map_dbl(function(x){
      rnorm(n = 1, mean = x, sd = (1+x)*1/nu) %>%
        (function(z) ifelse(z < 0, 0, z)) 
    }) %>%
    replace(zero_idx, 0)
  
  tibble(x = xs,
         trait = mu_spp, 
         mu = mu, 
         spp = rep(species_id, length(mu))
  )
}



spp_order <- load_emery() %>% 
 select(Species, sppint) %>%
 unique() %>%
 arrange(sppint) %>% 
 select(1) %>%
 as_vector()  

gen_stan_df <- function(stan_out, species_order){
  
  #stan_out output from gen_tolerance()
  #species_order vector of desired labels matching group_ids used in gen_tolerance()
  
  new_post <- stan_out %>% map(~extract(.x))
  
  params <- c("d", "e", "a", "b", "c", "beta_0", "beta_1", "nu")
  
  par_df <- new_post %>% map(function(y){
    params %>% map( ~{
      y[[.x]] %>%
        data.frame() %>%
        set_colnames(species_order) %>%
        gather("Species", x) %>%
        set_colnames(c("Species", .x)) %>%
        group_by(Species) %>%
        mutate(draw = 1:n())
    }) %>% 
      do.call(cbind, .) %>%
      select(-starts_with("Species")) %>%
      select(-matches("draw[0-9]")) %>%
      mutate(
        maxima = (((a - 1)/(a*b - 1))^(1/a) * (e - d) + d),
        breadth = (e-d),
        area = c * breadth,
        special = c/breadth
      )
  })
  
  par_df
}  


######################################
######################################
#### SIMULATE POSITION ALONG AXIS ####
######################################
######################################

rand_axis <- function(reps = 1){
  emery <- load_emery()
  s1 <- rgamma(n = 1, shape = 2, scale = 10)
  s2 <- rgamma(n = 1, shape = 2, scale = 10)
  treat_map <- rbeta(n = 5, shape1 = s1, shape2 = s2) %>% 
    sort() %>% (function(x) 1 + (x - min(x)) * ( 5 - 1 ) / (max(x) - min(x))) %>%
    tibble(treat = 1:5,
           rand_treat = .)
  left_join(emery, treat_map, by = "treat")
}



########################################
########################################
### LOAD AND PREP POOL GRADIENT DATA ###
########################################
########################################

load_gradient <- function(){
  emery <- load_emery()
  gradient <- read.xls("data/Pool depths_FINAL summary_REVISED.xls", 
                       header = T, skip = 1, stringsAsFactors = F) %>%
    mutate(taxa = strsplit(X, "_") %>% map_chr(~ .x[length(.x)]),
           taxa = ifelse(taxa == "deblilis", "debilis", taxa)) %>% 
    filter(taxa %in% as.character(unique(emery$Species)))
  
  droppers <- which(gradient %>% select(-taxa, -X) %>% apply(1, function(x) mean(is.na(x))) == 1)
  gradient[-droppers,] %>% 
    select(Mean, SE, taxa) %>%
    rename(Species = taxa)
}

#################################################
#################################################
### CONSTRUCT STAN AND GRADIENT COMBINED DATA ###
#################################################
#################################################

gen_gradient_df <- function(tolerance_df){
  
  grad <- load_gradient()
  #reg_df is constructed above
  wide_params <- full_join(x = tolerance_df, y = grad, by = "Species") %>%
    full_join(., reg_df, by = "Species")
  
  wide_params <- wide_params %>%
    ungroup() %>%
    select(-c(Species, draw, Mean, SE, habit, aqua_terr2terr, aqua_terr2vernal)) %>% 
    scale %>% 
    as_tibble %>%
    set_colnames( paste0(names(.), "_sc" )) %>%
    bind_cols(., wide_params) %>%
    mutate(
      aqua_terr2terr_bin = ifelse(
        aqua_terr2terr == "vernal", 0, 1
      ),
      aqua_terr2vernal_bin = ifelse(
        aqua_terr2vernal == "vernal", 0, 1
      )
    )
  
  wide_params
}
