#POSTERIOR PREDICTIVE CHECKS FOR LASTHENIA TOLERANCE CURVE MODEL
source("tolerance_functions.R")
#source("load_data.R")

emery <- load_emery() %>%
  select(Species, sppint, Inflor_biomass, treat)

zig_draws <- read.csv("bayes/stan_par1_df.csv") #draws for zero weighted model (generative)

post_pred <- 1:nrow(zig_draws) %>% 
  map_df(~{
    #post_pred <- 1:100 %>% map_df(~{
    draw_i <- zig_draws[.x, ]
    zig_spp <- zig_draws$Species[.x]
    emery_i <- emery %>% filter(Species == zig_spp)
  
    x <- (emery_i$treat - draw_i$d) / (draw_i$e - draw_i$d)
  
    mus <- stretch_kumara(
      x, a = draw_i$a, 
      b = draw_i$b, 
      c = draw_i$c
      )
  
  p_zero <- plogis(
    draw_i$beta_0 + draw_i$beta_1 * mus
    )
  
  #pseudo <- rgamma(n = length(mus), 
  #                 shape = draw_i$nu, 
  #                 rate =  draw_i$nu / mus)
  #zero <- rbinom(n = length(mus), size = 1, p_zero)
  #pseudo[as.logical(zero)] <- 0
  
  
  pseudo <- rgamma(
    n = length(mus), shape = draw_i$nu, 
    rate =  (draw_i$nu * (1 - p_zero)) / mus
    )
  
  data.frame(
    ssq_obs = sum((mus - emery_i$Inflor_biomass)^2),
    ssq_pseudo = sum((mus - pseudo)^2),
    Species = zig_spp
    #pseudo = pseudo,
    #Inflor_biomass = emery_i$Inflor_biomass,
    #Species = zig_spp,
    #draw = .x
  )

})

mean(post_pred$ssq_obs > post_pred$ssq_pseudo)

pred_spp <- post_pred %>% group_by(Species) %>%
  summarise(lab = as.character(mean(ssq_obs > ssq_pseudo) %>% round(2)),
            ssq_obs = ssq_obs %>% max*.9,
            ssq_pseudo = ssq_pseudo %>% max*.9)

write_csv(pred_spp, "bayes/postpred_zig1_out.csv")

post_plot <- post_pred %>% ggplot(aes(x = ssq_obs, y = ssq_pseudo)) +
  geom_point(alpha = 0.1, aes(colour = Species)) +
  geom_text(data = pred_spp, label = pred_spp$lab, size = 2) +
  facet_wrap(~Species, scales = "free") +
  geom_abline(intercept = 0, slope = 1, lty = 2, col = "blue") +
  ylab("Posterior predictive ssq") +
  xlab("Observed ssq") +
  theme_minimal() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  theme(text = element_text(size = 4))

ggsave("figures/figB12.pdf", plot = post_plot, device = "pdf")
#ggsave("analyses_and_viz/postpred.png", plot = post_plot, device = "png")



# Modeling checking instead using sampled posterior p-values

#Grab a posterior parameter vector at random
df_draw1 <- zig_draws %>% filter(
  draw == sample(1:max(zig_draws$draw), 1)
  )

#simulate 100 data sets from the random post vector
post_pred <- 
  1:500 %>% map_df(~{
    1:nrow(df_draw1) %>% 
      map_df(~{
        draw_i <- df_draw1[.x, ]
        zig_spp <- df_draw1$Species[.x]
        emery_i <- emery %>% filter(Species == zig_spp)
        
        x <- (emery_i$treat - draw_i$d) / (draw_i$e - draw_i$d)
        
        mus <- stretch_kumara(
          x, a = draw_i$a, 
          b = draw_i$b, 
          c = draw_i$c
        )
        
        p_zero <- plogis(
          draw_i$beta_0 + draw_i$beta_1 * mus
        )
        
        pseudo <- rgamma(
          n = length(mus), 
          shape = draw_i$nu, 
          rate =  (draw_i$nu * (1 - p_zero)) / mus
        )
        
        #pseudo <- rgamma(
        #  n = length(mus), 
        #  shape = draw_i$nu, 
        #  rate =  (draw_i$nu) / mus
        #)
        
        #pseudo[rbinom(length(mus), 1, (1 - p_zero))] <- 0
        
        
        data.frame(
          
          ssq_obs = sum((emery_i$Inflor_biomass - mus)^2)/mus,
          ssq_pseudo = sum((pseudo - mus)^2)/mus,
          
          Species = zig_spp
        )
        
      })
  })

#compare ssqs by species
post_pred %>% group_by(Species) %>%
  summarise(lab = as.character(mean(ssq_obs > ssq_pseudo) %>% round(2)),
            ssq_obs = ssq_obs %>% max*.9,
            ssq_pseudo = ssq_pseudo %>% max*.9)

#compare whole ssqs
mean(post_pred$ssq_obs > post_pred$ssq_pseudo)

