### Set-up
#---------

setwd("~/Documents/Projects/tolerance-curve2/")
source("load_data.R")
source("derived_files/lasth_100_post.R")

source("derived_files/state_reg.R") #for reg_df
draws <- read_csv("bayes/stan_par1_df.csv") #draws for penalized zero model

select <- dplyr::select

gradient <- read.xls("data/Pool depths_FINAL summary_REVISED.xls", 
                     header = T, skip = 1, stringsAsFactors = F) %>%
  mutate(taxa = strsplit(X, "_") %>% map_chr(~ .x[length(.x)]),
         taxa = ifelse(taxa == "deblilis", "debilis", taxa)) %>% 
  filter(taxa %in% as.character(unique(emery$Species)))

droppers <- which(gradient %>% select(-taxa, -X) %>% apply(1, function(x) mean(is.na(x))) == 1)
gradient <- gradient[-droppers,]


#nothing missing in either direction
tree_miss <- as.character(unique(emery$Species))[!as.character(unique(emery$Species)) %in% unique(gradient$taxa)]
unique(gradient$taxa)[!unique(gradient$taxa) %in% as.character(unique(emery$Species))]

col_match <- match(gradient$taxa, as.character(unique(emery$Species)))
col_match <- col_match[!is.na(col_match)]



### master dataframe with all parameters, gradient measures, and habitats
#--------------------------------------

grad <- gradient %>% select(Mean, taxa) %>%
  rename(Species = taxa)

draws_join <- full_join(x = draws, y = grad, by = "Species") %>%
  full_join(., reg_df, by = "Species") %>% 
  gather(param_name, param_value, 
         -c(draw, Species, Mean, habit, aqua_terr2terr, aqua_terr2vernal))

param_names <- c("area", "d", "maxima", "e", "c", "special", "breadth")

max_draws <- draws_join %>% 
  filter(param_name %in% param_names) %>%
  group_by(Species, param_name) %>%
  summarise_all(.funs = max)

#multivariate prediction of position along vernal pool
wide_params <- full_join(x = draws, y = grad, by = "Species") %>%
  full_join(., reg_df, by = "Species")

wide_params <- wide_params %>%
  select(-Species, -draw, -Mean, -habit, -aqua_terr2terr, -aqua_terr2vernal) %>% 
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


wide_params %>% 
  ggplot(aes(x = .$maxima_sc, y = .$e_sc, colour = Species)) +
  geom_point()

wide_params %>% 
  select(Species, maxima, area, d, e, breadth) %>%
  ggpairs(aes(colour = Species, alpha = 0.1))

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
          gls(Mean ~ c_sc + d_sc + e_sc, 
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
          names(intercept) <- "intercept_1"
          c(coef(mod_full), intercept, aic_f, aic_i) %>% 
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
    -draw, -pagel_lambda, 
    -`(Intercept)`, -intercept_1, 
    -AIC_full, -AIC_intercept
  )


options(xtable.sanitize.colnames.function=identity,
        xtable.sanitize.rownames.function=identity)

mv_plgs_long %>% 
  group_by(param, pagel_lambda) %>%
  summarise(means = mean(slopes),
            prob_zero = mean(slopes > 0)) %>%
  mutate(pagel_lambda = as.integer(pagel_lambda)) %>% 
  set_colnames(c("parameter", "Pagel's $\\lambda$", "coefficient mean", "proportion of draws $>$ 0")) %>%
  as.data.frame %>% 
  set_rownames(c("$\\zeta $", "$\\zeta$", 
                 "$\\delta $", "$\\delta$", 
                 "$\\epsilon $", "$\\epsilon$")) %>%
  select(-parameter) %>%
  xtable(caption = "something here", digits = 3) %>% 
  print.xtable()

  

#how do aic values and predictions compare for full models, but alternate pagel lambdas?
pg_0 <- mv_plgs_long %>% 
  filter(pagel_lambda == 0)

pg_1 <- mv_plgs_long %>% 
  filter(pagel_lambda == 1)

c("AIC_full", "AIC_intercept", "slopes")   %>%
  map(~ mean(pg_0[[.x]] < pg_1[[.x]]))

mean(pg_0$AIC_full < pg_1$AIC_full)
mean(pg_0$slopes > pg_1$slopes)


mv_plgs_long %>%
  select(AIC_full, AIC_intercept, pagel_lambda) %>%
  gather("model", "aic", -pagel_lambda) %>%
  mutate(pagel_lambda = paste0("λ = ", pagel_lambda),
         model = gsub("_", " ", model)) %>%
  ggplot(aes(x = aic, y = model)) +
  facet_wrap(~factor(pagel_lambda), ncol = 1) +
  geom_joy(col = "white") +
  theme_bw() +
  theme(text = element_text(size=16))


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
      coef(s = "lambda.min")
    
    lasso_df <- rep(0, length(cv_coef@Dimnames[[1]])) %>% rbind %>% 
      as_tibble %>%
      set_colnames(cv_coef@Dimnames[[1]])
    
    covs <- cv_coef@Dimnames[[1]]
    cov_idx <- match(covs[cv_coef@i+1],
          covs
          )
    lasso_df[cov_idx] <- cv_coef@x
    lasso_df
  })


lasso_all %>% 
  ungroup %>%
  select(-draw, -`(Intercept)`) %>%
  summarise_each(
    funs(
      mean = mean(., na.rm = T),
      prop = mean(. > 0, na.rm = T)
      )
    ) %>% 
  gather(variable, value) %>%
  separate(variable, c("var", "stat"), sep = "_sc\\_") %>%
  spread(var, value) %>%
  set_colnames(c(" ", "$\\zeta$", "$\\delta$", "$\\epsilon$", "maxima")) %>%
  mutate(" " = c("coefficient mean", "proportion of draws > 0")) %>%
  xtable(caption = "something here", digits = 3) %>% 
  print.xtable(include.rownames=FALSE)


#glm(aqua_terr2terr_bin ~  c_sc + breadth_sc, 
#    data = test_df, family = binomial)

###
###
###

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
  
#cor_c_breadth <- pic_cor("c", "breadth")
#cor_d_breadth <- pic_cor("d", "breadth")
#cor_e_breadth <- pic_cor("e", "breadth")
#cor_e_maxima <- pic_cor("e", "maxima")
#cor_d_maxima <- pic_cor("d", "maxima")
#cor_c_maxima <- picklopioipi,,o,_cor("c", "maxima")
#cor_c_d <- pic_cor("c", "d")
#cor_d_e <- pic_cor("d", "e")

cor_c_e <- pic_cor("c", "e")

cor_c_e %>% 
  gather(corr_type, corr, -draw) %>%
  group_by(corr_type) %>%
  summarise(mean(corr),
            mean(corr > 0))
hist(cor_c_e$pic_cor, breaks = 1000)
