source("tolerance_functions.R")

emery <- load_emery()

#RUN STAN
stan_out <- gen_tolerance(df = emery, file_id = "")

summary(stan_out[[1]])$summary

stan_out %>% map(~ print(.x, pars = c("e", "nu")))

traceplot(stan_out[[1]], pars = "e_t")

par_df <- gen_stan_df(stan_out, spp_order)

par_set1 <- par_df[[1]] %>% group_by(Species) %>%
  map_kumara(seq(0,1, length.out = 100), .)

par_set2 <- par_df[[2]] %>% group_by(Species) %>%
  map_kumara2(seq(0,1, length.out = 100), .)

qs <- seq(0,1, length.out = 1000)

q1 <- quantile(par_set1$y, qs)
q2 <- quantile(par_set2$y, qs)
plot(q1, q2)
abline(0,1, lwd = 2, col = "red")
cor(q1, q2)

write.csv(par_df[[1]], file = "bayes/stan_par1_df.csv", quote = F, row.names = F)
write.csv(par_set1, file = "bayes/fitted_points_mod1.csv", quote = F, row.names = F)
write.csv(par_df[[2]], file = "bayes/stan_par2_df.csv", quote = F, row.names = F)
write.csv(par_set2, file = "bayes/fitted_points_mod2.csv", quote = F, row.names = F)


###Delete here down?

# emery2 <- emery %>% select(treat, Inflor_biomass, Species)
# ggplot() +
#  geom_line(data = par_set1, aes(x=x, y=y, group=as.character(draw)), alpha = 0.05, colour = "dodgerblue") +
#  geom_line(data = par_set2, aes(x=x, y=y, group=as.character(draw)), alpha = 0.05, colour = "yellow") +
#  geom_jitter(data = emery, mapping = aes(x = treat, y = Inflor_biomass),
#              height = 0, width = 0.1, alpha = 0.5) +
#  stat_smooth(data = emery, aes(x=treat, y=Inflor_biomass), se = F, colour = "green") +
#  facet_wrap(~Species, scales = "free") +
#  theme_minimal()


#new_post[[1]]$b %>%
#  data.frame() %>%
#  set_colnames(spp_order) %>% 
#  gather("Species", "par") %>%
#  ggplot(aes(y = Species, x = par)) +
#  geom_joy(alpha = 0.5, colour = "white") + 
#  theme_joy(grid = FALSE) +
#  scale_x_continuous(expand = c(0.01, 0)) +
#  scale_y_discrete(expand = c(0.01, 0))

#par_df[[1]] %>% 
#  ggplot(aes(colour = Species, x = e)) +
#  geom_density()
  

#min_max_df <- emery %>%
#  group_by(sppint_h) %>%
#  summarise(minx = min(treat),
#            maxx = max(treat))

#stan_out_hologymne <- c(0,1) %>% map(~{
  
#  stan_in <- list(
#    N = length(emery$Inflor_biomass), 
#    y = emery$Inflor_biomass, 
#    x = emery$treat,
#    minx = c(min_max_df$minx),
#    maxx = c(min_max_df$maxx),
#    numSpp = length(unique(emery$sppint_h)),
#    sppint = emery$sppint_h,
#    zig = .x
#  )
  
#  sampling(mod,
#           data = stan_in,
#           control = list(adapt_delta = 0.8, max_treedepth = 10), 
#           iter = 2000, chains = 4, 
#           sample_file = paste0("bayes/samples/tolerance_v3_zig_",.x,".samples"),
#           seed = 23265
#  )
#})


#library(loo)
#log_lik_1 <- extract_log_lik(stan_out[[2]])
#loo_1 <- loo(log_lik_1)
#waic_1 <- waic(log_lik_1)
#print(loo_1)
#emery[loo_1$pareto_k > 0.5,]

#log_lik_h <- extract_log_lik(stan_out_hologymne[[1]])
#loo_h <- loo(log_lik_h)
#waic_h <- waic(log_lik_h)
#print(loo_h)
#emery[loo_h$pareto_k > 0.5,]


#compare(loo_h, loo_1)
#compare(waic_h, waic_1)
