#POSTERIOR PREDICTIVE CHECKS FOR LASTHENIA TOLERANCE CURVE MODEL
source("tolerance_functions.R")
#source("load_data.R")

emery <- load_emery() %>%
  select(Species, sppint, Inflor_biomass, treat)

zig_draws <- read.csv("bayes/stan_par1_df.csv") #draws for zero weighted model (generative)

post_pred <- 1:nrow(zig_draws) %>% map_df(~{
#post_pred <- 1:100 %>% map_df(~{
  draw_i <- zig_draws[.x, ]
  zig_spp <- zig_draws$Species[.x]
  emery_i <- emery %>% filter(Species == zig_spp)
  
  x <- (emery_i$treat - draw_i$d) / (draw_i$e - draw_i$d)
  
  mus <- stretch_kumara(x, a = draw_i$a, 
                        b = draw_i$b, 
                        c = draw_i$c)
  
  p_zero <- plogis(draw_i$beta_0 + 
                     draw_i$beta_1 * mus)
  
  #pseudo <- rgamma(n = length(mus), 
  #                 shape = draw_i$nu, 
  #                 rate =  draw_i$nu / mus)
  #zero <- rbinom(n = length(mus), size = 1, p_zero)
  #pseudo[as.logical(zero)] <- 0
  
  
  pseudo <- rgamma(n = length(mus), shape = draw_i$nu, 
                   rate =  (draw_i$nu * (1 - p_zero)) / mus)
  
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
            ssq_obs = ssq_obs %>% max*.75,
            ssq_pseudo = ssq_pseudo %>% max*.75)

write_csv(pred_spp, "bayes/postpred_zig1_out.csv")

post_plot <- post_pred %>% ggplot(aes(x = ssq_obs, y = ssq_pseudo)) +
  geom_point(alpha = 0.1, aes(colour = Species)) +
  geom_text(data = pred_spp, label = pred_spp$lab) +
  facet_wrap(~Species, scales = "free") +
  geom_abline(intercept = 0, slope = 1, lty = 2, col = "blue") +
  ylab("Posterior predictive ssq") +
  xlab("Observed ssq") +
  theme_minimal() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  theme(text = element_text(size = 5))

ggsave("analyses_and_viz/postpred.pdf", plot = post_plot, device = "pdf")
ggsave("analyses_and_viz/postpred.png", plot = post_plot, device = "png")