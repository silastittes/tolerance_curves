#---------
source("load_data.R")

select <- dplyr::select
rename <- dplyr::rename


#nothing missing in either direction
tree_miss <- as.character(unique(emery$Species))[!as.character(unique(emery$Species)) %in% unique(grad$Species)]
unique(grad$Species)[!unique(grad$Species) %in% as.character(unique(emery$Species))]

col_match <- match(grad$Species, as.character(unique(emery$Species)))
col_match <- col_match[!is.na(col_match)]

wide_params <- gen_gradient_df(tolerance_df = draws)


### master dataframe with all parameters, gradient measures, and habitats
#--------------------------------------

mv_plgs <- wide_params %>% 
  group_by(draw) %>%
  do({
    temp_df <- .
    trees_post %>% map_df(~{
      0:1 %>% map_df(function(x){
        mv_pgls_tree <- drop.tip(.x, tree_miss)
        
        mod_full <- temp_df %>%
          filter(Species %in% grad$Species) %>% 
          as.data.frame %>%
          set_rownames(.$Species) %>%
          gls(Mean ~ c_sc + d_sc + e_sc +  maxima_sc, #keep maxima_sc? 
              data = .,
              correlation = corPagel(x, phy = mv_pgls_tree, fixed = T),
              method = "ML"
              )
        
        mod_intercept <- temp_df %>%
          filter(Species %in% grad$Species) %>% 
          as.data.frame %>%
          set_rownames(.$Species) %>%
          gls(Mean ~ 1, 
              data = .,
              correlation = corPagel(x, phy = mv_pgls_tree, fixed = T),
              method = "ML"
          )
        
          aic_f <- AIC(mod_full)        
          aic_i <- AIC(mod_intercept)
          names(aic_f) <- "AIC_full"
          names(aic_i) <- "AIC_intercept"
          intercept <- coef(mod_intercept)
          pvals <- summary(mod_full)$tTable[,4]
          names(pvals) <- paste0(names(pvals), "_pvals")
          names(intercept) <- "intercept_1"
          c(coef(mod_full), pvals, intercept, aic_f, aic_i) %>% 
          rbind %>% 
          cbind(., x) %>% 
          as_tibble %>% 
          rename(pagel_lambda = x)
        })
      })
    })


mv_plgs %>%
  group_by(pagel_lambda) %>%
  summarise(mean(AIC_full < AIC_intercept))


mv_plgs_long <- mv_plgs %>% ungroup %>% 
  gather(
    param, slopes, 
    `(Intercept)`, c_sc, d_sc, e_sc, maxima_sc
  ) %>% 
  gather(params, pvals,
         `(Intercept)_pvals`,
         c_sc_pvals,
         d_sc_pvals,
         e_sc_pvals,
         maxima_sc_pvals
         )

mv_plgs_long %>% 
  group_by(pagel_lambda, params) %>%
  summarise(q_25 = quantile(pvals, 0.25) %>% round(2),
            q_50 = quantile(pvals, 0.5) %>% round(2),
            q_75 = quantile(pvals, 0.75) %>% round(2),
            harmonic_mean = round(1/mean(1/pvals), 4)) %>%
  arrange(params)

mv_plgs_long %>% 
  group_by(pagel_lambda, params) %>%
  summarise(hpdi_low = HPDI(pvals, 0.95)  %>% .[1] %>% round(2),
            q_50 = quantile(pvals, 0.5) %>% round(2),
            hpdi_high = HPDI(pvals, 0.95)  %>% .[2] %>% round(2))

pval_hm <- mv_plgs_long %>% 
  group_by(pagel_lambda, params) %>%
  summarise(harmonic_mean = round(1/mean(1/pvals), 5)) %>%
  arrange(params) %>%
  ungroup()

cairo_pdf(filename = "figures/B13.pdf")

mv_plgs_long %>% 
  mutate(params = str_replace_all(params, 
                                  c("^c_sc_pvals$" = "ζ", 
                                    "\\(Intercept\\)_pvals" = "Intercept",
                                    "^maxima_sc_pvals$" = "optima",
                                    "^d_sc_pvals$" = "δ", 
                                    "^e_sc_pvals$" = "ε")
                                  ),
         pagel_lambda = str_replace_all(pagel_lambda, 
                                 c("0" = "λ = 0",
                                   "1" = "λ = 1"))
  ) %>%
  ggplot(aes(x = pvals)) +
  facet_grid(params ~ pagel_lambda, scales = "free_y") +
  geom_histogram() +
  theme_minimal()
dev.off()

options(xtable.sanitize.colnames.function=identity,
        xtable.sanitize.rownames.function=identity)

mv_plgs_long %>% 
  group_by(param, pagel_lambda) %>%
  summarise(
    means = mean(slopes),
    prob_zero = mean(slopes > 0)) %>%
  ungroup() %>%
  mutate(
    pagel_lambda = as.integer(pagel_lambda),
    pval_hm = pval_hm$harmonic_mean
    ) %>%
  set_colnames(c(
    "parameter", 
    "Pagel's $\\lambda$", 
    "coefficient mean",
    "proportion of draws $>$ 0", 
    "harmonic mean of p values")
    ) %>%
  as.data.frame %>% 
  filter(parameter != "(Intercept)") %>%
  set_rownames(c("$\\zeta $", "$\\zeta$", 
                 "$\\delta $", "$\\delta$", 
                 "$\\epsilon $", "$\\epsilon$",
                 "$optima $", "$optima$")) %>%
  select(-parameter) %>%
  xtable(caption = "something here", digits = 4) %>% 
  print.xtable(file = "derived_files/pgls.tex")


#how do aic values and predictions compare for full models, but alternate pagel lambdas?
pg_0 <- mv_plgs_long %>% 
  filter(pagel_lambda == 0)

pg_1 <- mv_plgs_long %>% 
  filter(pagel_lambda == 1)

c("AIC_full", "AIC_intercept", "slopes")   %>%
  map(~ mean(pg_0[[.x]] < pg_1[[.x]]))

mean(pg_0$AIC_full < pg_1$AIC_full)
mean(pg_0$slopes > pg_1$slopes)

#library(ggjoy)
#mv_plgs_long %>%
#  select(AIC_full, AIC_intercept, pagel_lambda) %>%
#  gather("model", "aic", -pagel_lambda) %>%
#  mutate(pagel_lambda = paste0("λ = ", pagel_lambda),
#         model = gsub("_", " ", model)) %>%
#  ggplot(aes(x = aic, y = model)) +
#  facet_wrap(~factor(pagel_lambda), ncol = 1) +
#  geom_joy(col = "white") +
#  theme_bw() +
#  theme(text = element_text(size=16))


cairo_pdf(filename = "figures/B14.pdf")
layout(
  mat = matrix(c(1,2,3,4), byrow = F, nrow = 2)
)

rdraw <- 155
rtree_int <- 24
for(i in 0:1){
  temp_df <- wide_params %>% 
    filter(draw == rdraw)
  
  mv_pgls_tree <- drop.tip(trees_post[[rtree_int]], tree_miss)
  
  mod_full <- temp_df %>%
    filter(Species %in% grad$Species) %>% 
    as.data.frame %>%
    set_rownames(.$Species) %>%
    gls(Mean ~ c_sc + d_sc + e_sc +  maxima_sc, #keep maxima_sc? 
        data = .,
        correlation = corPagel(i, phy = mv_pgls_tree, fixed = T),
        method = "ML"
    )
  print(summary(mod_full))
  
  st_resid <- mod_full$residuals/attr(mod_full$residuals, "std")[1]
  plot(
    mod_full$fitted,
    st_resid,
    main = paste("λ = ", i),
    xlab = "fitted values",
    ylab = ifelse(i == 0, "standardized residuals", "" )
    )
  
  lines(
    loess.smooth(
      mod_full$fitted,
      st_resid
    )
  )
  abline(h = 0, lty = 2)
  qqnorm(
    st_resid, 
    main = "",
    ylab = ifelse(i == 0, "Sample Quantiles", "" )
    )
  abline(0,1, lty = 2)
}

dev.off()


#NON PHYLOGENETIC LASSO CLASSIFICATION BY HABITAT
lasso_all <- wide_params %>% group_by(draw) %>%
  do({
    test_df <- . 
    test_df %>% as.data.frame %>% set_rownames(.$Species)
    test_modmat <- test_df %>% select(contains("_sc")) %>%
      select(c_sc, maxima_sc, e_sc, d_sc) %>%
      model.matrix(test_df$aqua_terr2terr_bin ~ . -1, data = .)

    cv_coef <- cv.glmnet(test_modmat, 
                         test_df$aqua_terr2terr_bin,
                         family = "binomial", 
                         nfolds = nrow(test_modmat),
                         grouped = F) %>%
      coef(s = "lambda.min") %>% 
      as.matrix %>% 
      t %>% 
      as_tibble
  })


lasso_all %>% 
  ungroup %>%
  select(-draw, -`(Intercept)`) %>%
  summarise_all(
    funs(
      mean = mean(., na.rm = T),
      prop = mean(. > 0, na.rm = T)
      )
    ) %>% 
  gather(variable, value) %>%
  separate(variable, c("var", "stat"), sep = "_sc\\_") %>%
  spread(var, value) %>%
  set_colnames(c(" ", "$\\zeta$", "$\\delta$", "$\\epsilon$", "$optima$")) %>%
  mutate(" " = c("coefficient mean", "proportion of draws > 0")) %>%
  xtable(caption = "something here", digits = 3) %>% 
  print.xtable(include.rownames=FALSE, file = "derived_files/lasso.tex")




#r_draw <- 353
r_draw <- sample(1:400, 1)
test_df <- wide_params %>% 
  filter(draw == r_draw) %>% 
  as.data.frame %>% 
  set_rownames(.$Species)
test_modmat <- test_df %>% select(contains("_sc")) %>%
  select(c_sc, maxima_sc, e_sc, d_sc) %>%
  model.matrix(test_df$aqua_terr2terr_bin ~ . -1, data = .)

cv_coef <- cv.glmnet(test_modmat, 
                     test_df$aqua_terr2terr_bin,
                     family = "binomial", 
                     nfolds = nrow(test_modmat),
                     grouped = F, 
                     type.measure = "deviance"
                     )

#cairo_pdf(filename = "figures/B15.pdf")
plot(cv_coef)
#dev.off()

as.matrix(coef(cv_coef, s = "lambda.min"))

#glm(aqua_terr2terr_bin ~  c_sc + breadth_sc, 
#    data = test_df, family = binomial)

### Phylogenetic Independent Contrasts 
#--------------------------------------
pic_cor <- function(var1, var2){
  wide_params %>% 
    group_by(draw) %>%
    do({
      temp_df <- .
      trees_post %>% map_df(~{
        
        lasth_pic <- .x

        pic_cor <- cor(
          pic(temp_df[[var1]], lasth_pic),
          pic(temp_df[[var2]], lasth_pic), 
          method = "pearson"
        ) %>% as_tibble %>%
          rename(pic_cor = value)
        
        cor_cor <- cor(
          temp_df[[var1]],
          temp_df[[var2]], 
          method = "pearson"
        ) %>% as_tibble %>%
          rename(cor_cor = value)
        
        bind_cols(pic_cor, cor_cor)
      })
    })
} 
  

cor_c_e <- pic_cor("c", "e")

cor_c_e %>% 
  gather(corr_type, corr, -draw) %>%
  group_by(corr_type) %>%
  summarise(mean(corr),
            mean(corr > 0))

hist(cor_c_e$pic_cor, breaks = 500)


#old stuff, can probably delete
# 
# draws_join <- full_join(x = draws, y = grad, by = "Species") %>%
#   full_join(., reg_df, by = "Species") %>% 
#   gather(param_name, param_value, 
#          -c(draw, Species, Mean, habit, aqua_terr2terr, aqua_terr2vernal))
# 
# param_names <- c("area", "d", "maxima", "e", "c", "special", "breadth")
# 
# max_draws <- draws_join %>% 
#   filter(param_name %in% param_names) %>%
#   group_by(Species, param_name) %>%
#   summarise_all(.funs = max)

#wide_params %>% 
#  ggplot(aes(x = .$maxima_sc, y = .$e_sc, colour = Species)) +
#  geom_point()

#wide_params %>% 
#  select(Species, maxima, area, d, e, breadth) %>%
#  ggpairs(aes(colour = Species, alpha = 0.1))
