#GENERATE DERIVED PARAMETERS AND ANALYSES FOR LASTHENIA TOLERANCE CURVE DATA
source("load_data.R")
source("derived_files/state_reg.R")

post_trees <- source("derived_files/lasth_100_post.R")$value

draws <- read.csv("bayes/stan_par1_df.csv") #draws for penalized zero model
draw_fits <- read.csv("bayes/fitted_points_mod1.csv")

n_draws <- posts$lp__ %>% nrow


###############################################################
##CURVE SIGNAL-------------------------------------------------
###############################################################


xseq_comm <- seq(1,5,length.out=100)
signal_set1 <- draws %>% revmap_kumara(xseq_comm, .)

curve_sig <- 1:n_draws %>%
  map(~{
    d2 <- signal_set1 %>% 
      filter(draw == .x) %>%
      mutate(y = ifelse(is.na(y), 0, y))
    out <- abind(split(d2[,2:3], d2$Species), along=3)
    post_trees %>% map_dbl(function(lasth_i){
      physignal(A = out, phy = lasth_i, iter = 1)$phy.signal  
    })
  }) %>% 
  flatten %>% 
  as_vector %>% 
  as_tibble %>%
  mutate(param = rep("curve", n()),
         tree = rep(1:length(post_trees), n_draws)) %>%
  rename(signal = value)


###############################################################
##PARAM SIGNAL-------------------------------------------------
###############################################################


param_names <- c("area", "d", "maxima", "e", "c", "breadth")
signal_df_pre <- param_names %>% map(function(var){
  1:n_draws %>% map(~{
    param_df <- draws %>% 
      group_by(Species) %>%
      filter(draw == .x) %>% 
      select(Species, one_of(var))
    param <- param_df[[2]]
    names(param) <- param_df[[1]]
    post_trees %>% map_dbl(function(lasth_i){
      physignal(A = param, phy = lasth_i, iter = 1)$phy.signal  
    })
  })
}) %>% do.call(cbind, .) %>% 
  as_tibble %>% 
  set_colnames(param_names) %>%
  gather("param", "signal")

signal_df <- signal_df_pre$signal %>%
  flatten %>% 
  as_vector %>% 
  as_tibble %>%
  mutate( param = rep(signal_df_pre$param, each = length(trees_post)),
          tree = rep(1:length(post_trees), length(param_names)*n_draws)) %>%
  rename(signal = value)


all_sig <- bind_rows(curve_sig, signal_df) %>%
  mutate(param = factor(param))

levels(all_sig$param) <- c("curve", param_names)

write_csv(all_sig, "derived_files/curve_K.csv")
#write.table(x = curveK_draws, file = "derived_files/curve_K.txt")
