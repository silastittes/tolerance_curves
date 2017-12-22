#load functions and data comment once loaded to allow seed setting here
#setwd("~/Documents/Projects/tolerance-curve2/")
source("tolerance_functions.R")

mod <- rstan::stan_model("bayes/tolerance_v3_alt.stan")

#Function to generate data
sim_tolerance <- function(
  n_spp_treat = 10,
  n_treat = 5,
  n_spp = 16,
  
  hyper_means = list(
    a = -1,
    b = -1,
    c = -1,
    d = 0,
    e = 0,
    beta_0 = 0,
    beta_1 = -1.5,
    nu = 20
  ),
  
  hyper_sds = list(
    a = 0.2, 
    b = 0.2, 
    c = 0.2, 
    d = 0.2,
    e = 0.2,
    beta_0 = 0.1,
    beta_1 = 0.1,
    nu = 0.1
  ),
  
  random_axis = F,
  shift = 2
  
){
  
  exp_tr <- function(x, lim, upper = T){
    if(upper) lim + exp(-x) else lim - exp(-x)
  }
  
  #n_spp_treat = 2
  #n_treat = 20
  #n_spp = 3
  
  nobs <- n_spp*n_treat*n_spp_treat
  
  x_o <- 1:n_spp %>% map(~seq(1, n_treat, length.out = n_treat))
  
  if(random_axis){
    x_o <- 1:n_spp %>% map(~seq(1, n_treat, length.out = n_treat) + rnorm(1, 0, shift))
  }
  
  
  gen_val <- function(param_chr, lim, upper, n_spp){
    hyper <- rnorm(n = 1, mean = 0, sd = 1)
    #hyper_sd <- rtruncnorm(n = 1, mean = 0, sd = 1, a = 0)
    hyper_sd <- 1
    return(
      exp_tr(
        hyper_means[[param_chr]] +
          hyper_sds[[param_chr]] *
          rnorm(n = n_spp, mean = hyper, sd = hyper_sd),
        lim, upper
      ))
  }
  
  a <- gen_val(param_chr = "a", lim = 1, upper = T, n_spp)
  b <- gen_val(param_chr = "b", lim = 1, upper = T, n_spp)
  c <- gen_val(param_chr = "c", lim = 0, upper = T, n_spp)
  d <- x_o %>% map_dbl(~ gen_val(param_chr = "d", lim = min(.x), upper = F, n_spp = 1))
  e <- x_o %>% map_dbl(~ gen_val(param_chr = "e", lim = max(.x), upper = T, n_spp = 1))
  
  beta_0_hyper <- rnorm(n = 1, mean = 0, sd = 1)
  #beta_0_hyper_sd <- rtruncnorm(n = 1, mean = 0, sd = 1, a = 0)
  beta_0_hyper_sd <- 1
  beta_0 <- hyper_means$beta_0 + 
    hyper_sds$beta_0*
    rnorm(n = n_spp, mean = beta_0_hyper, sd = beta_0_hyper_sd)
  
  beta_1_hyper <- rtruncnorm(n = 1, mean = 0, sd = 1, b = 0)
  #beta_1_hyper_sd <- rtruncnorm(n = 1, mean = 0, sd = 1, a = 0)
  beta_1_hyper_sd <- 1
  beta_1 <- hyper_means$beta_1 + 
    hyper_sds$beta_1*
    rnorm(n = n_spp, mean = beta_1_hyper, sd = beta_1_hyper_sd)
  
  #nu_hyper <- rtruncnorm(n = 1, mean = hyper_means$nu, sd = hyper_sds$nu, a = 0)
  #nu_hyper_sd <- rtruncnorm(n = 1, mean = 0, sd = 1, a = 0)
  nu_hyper <- rtruncnorm(n = 1, mean = 1, sd = 1, a = 0)
  #nu_hyper_sd <- rtruncnorm(n = 1, mean = 0, sd = 1, a = 0)
  nu_hyper_sd <- 1
  nu <- rgamma(n = n_spp, shape =  nu_hyper, rate = nu_hyper_sd)
  #convert data scale to supported kimurswamy scale using d and e
  #produce mean y values
  
  x <- 1:n_spp %>% map(~(x_o[[.x]] - d[.x]) / (e[.x] - d[.x]))
  
  mus <- 1:n_spp %>% map(~stretch_kumara(x[[.x]], a[.x], b[.x], c[.x]))
  
  p_zero <- 1:n_spp %>% map(~plogis(beta_0[.x] + beta_1[.x] * mus[[.x]]))
  
  
  #generate gamma distributed variation around mean y values
  rgamma_v <- Vectorize(rgamma, "rate")
  rbern_v <- Vectorize(rbernoulli, "p")
  
  sim_data <- 1:n_spp %>% 
    map_df(~{
      #!!!
      #.x <- 1
      y <- rgamma_v(n = n_spp_treat, shape = nu[[.x]], rate = nu[[.x]] / mus[[.x]])
      zer <- rbern_v(n = n_spp_treat, p = p_zero[[.x]])
      y[zer] <- 0
      y %>% 
        data.frame() %>%
        gather(key = "treat_x", value = "y") %>% 
        mutate(treat = rep( x_o[[.x]], each = n_spp_treat),
               spp = rep(.x, nrow(.))) #%>%
      #select(-treat_x)
      
    })
  
  #ggplot(sim_data, aes(x = treat, y = y, colour = factor(spp))) +
  #  geom_point()
  
  return(
    list(
      sim_data = sim_data, 
      a = a,
      b = b,
      c = c,
      d = d,
      e = e,
      beta_0 = beta_0,
      beta_1 = beta_1,
      nu = nu
    )
  )
}


#function to parse and fit stan to simulated data
stan_tolerance <- function(sim_data){
  ## simulate parameters and data ----------------------------
  #sim_data <- simulate_tolerance(n_spp_treat, n_treat, n_spp)
  a <- sim_data$a
  b <- sim_data$b
  c <- sim_data$c
  d <- sim_data$d
  e <- sim_data$e
  beta_0 <- sim_data$beta_0
  beta_1 <- sim_data$beta_1
  nu <- sim_data$nu
  #truncate ends to exclude zero-only treatments  
  sim_data <- sim_data$sim_data
  sim_data <- sim_data %>% group_by(spp, treat) %>%
    filter(sum(y > 0) > 0)
  
  
  min_max_df <- sim_data %>% group_by(spp) %>%
    summarise(minx = min(treat),
              maxx = max(treat))
  
  # estimate parameters
  stan_sim <- zigs <- 0:1 %>% map(~{
    dataList <- list(N = length(sim_data$spp),
                     y = sim_data$y,
                     x = sim_data$treat,
                     minx = array(min_max_df$minx),
                     maxx = array(min_max_df$maxx),
                     numSpp = length(unique(sim_data$spp)),
                     sppint = sim_data$spp,
                     zig = .x,
                     
                     a_pr_mu = -1, a_pr_sig = 0.2,
                     b_pr_mu = -1, b_pr_sig = 0.2,
                     c_pr_mu = 0, c_pr_sig = 0.4,
                     d_pr_mu = 0, d_pr_sig = 0.4,
                     e_pr_mu = 0, e_pr_sig = 0.4)
    
    sampling(mod,
             data = dataList, 
             iter = 100, chains = 4, 
             control = list(adapt_delta = 0.8, max_treedepth = 10),
             seed = 123
    )
    
  })
  
  post <- stan_sim %>% map(~ rstan::extract(.x))
  
  params <- c("d", "e", "a", "b", "c", "beta_0", "beta_1", "nu")
  
  par_df <- post %>% map(function(y){
    params %>% map( ~{
      y[[.x]] %>%
        data.frame() %>%
        gather("Species", x) %>%
        set_colnames(c("Species", .x)) %>% 
        group_by(Species) %>%
        mutate(draw = 1:n())
    }) %>% 
      do.call(cbind, .) %>%
      subset(., select=which(!duplicated(names(.)))) 
  }) 
  
  par_set1 <- par_df[[1]] %>% group_by(Species) %>%
    map_kumara(seq(0,1, length.out = 100), .)
  
  par_set2 <- par_df[[2]] %>% group_by(Species) %>%
    map_kumara2(seq(0,1, length.out = 100), .)
  
  qq <- qqplot(par_set1$y, par_set2$y, plot.it = F)
  corr_models <-  cor(qq$x, qq$y)
  diff_models <- sum((qq$x - qq$y)^2) / (2*length(qq$y))
  
  #plot predictions of two models!
  true_ps <- c("a", "b", "c", 
               "d", "e", 
               "beta_0", "beta_1", 
               "nu")
  #sim_data %>% group_by(spp) %>% summarise(min(y))
  
  prop_true <- 1:length(unique(sim_data$spp)) %>% 
    map_df( ~ {
      seq_along(true_ps) %>% map_dbl(function(y){
        tp <- true_ps[y]
        t_vals <- eval(parse(text = tp))
        mean(t_vals[.x] > post[[2]][[tp]][,.x]) ##!!????!?!
      }) %>% 
        rbind() %>% 
        data.frame()
    })
  
  valid_df <- true_ps %>% map_df(~ {
    param <- .x
    param_true = eval(parse(text = param))
    true_df <- data.frame(
      param_true = eval(parse(text = param)), 
      spp = paste0("V", 1:length(param_true)), 
      stringsAsFactors = F
    )
    
    post[[2]][[param]] %>%
      as_tibble %>% 
      gather("spp", "param") %>%
      full_join(., true_df, by = "spp") %>%
      group_by(spp) %>%
      summarise(prop_true = mean(param_true > param),
                cred_in = quantile(param, 0.025) < unique(param_true) & 
                  quantile(param, 0.975) > unique(param_true)) %>%
      ungroup() %>%
      summarise(prop_true = mean(prop_true),
                cred_in = mean(cred_in))
  })
  
  
  prop_zero <- sim_data %>%
    group_by(spp) %>%
    summarise(prop_zero = mean(y == 0))
  #print("prop zero is:")
  #print(prop_zero)
  
  
  results <- list(
    stan_posts = post, #list of 2 model fits
    prop_zero = prop_zero$prop_zero,
    prop_true = prop_true,
    corr_models = corr_models,
    diff_models = diff_models,
    valid_df = valid_df
  )
  return(results)
}


#function to visualize data, fit, and gam fit
plot_sim_out <- function(out){
  
  #Make the "true" curves for each species
  true_tol <- seq_along(out$sim_data$a) %>% 
    map_df(~{
      plot_kumara(
        seq(0, 1, length.out = 100),
        out$sim_data$a[.x], 
        out$sim_data$b[.x], 
        out$sim_data$c[.x], 
        out$sim_data$d[.x], 
        out$sim_data$e[.x]
      ) %>% 
        do.call(cbind, .) %>%
        set_colnames(paste0(c("x", "y"))) %>%
        as_tibble %>%
        mutate(Species = paste0("X", .x))
    })
  
  
  params <- c("d", "e", "a", "b", "c", "beta_0", "beta_1", "nu")
  par_df <- out$sim_results$stan_posts %>% map(function(y){
    params %>% map( ~{
      y[[.x]] %>%
        data.frame() %>%
        gather("Species", x) %>%
        set_colnames(c("Species", .x)) %>%
        group_by(Species) %>%
        mutate(draw = 1:n())
    }) %>% 
      do.call(cbind, .) %>%
      subset(., select=which(!duplicated(names(.)))) 
  }) 
  
  par_set1 <- par_df[[1]] %>% group_by(Species) %>%
    map_kumara(seq(0,1, length.out = 100), .)
  
  par_set2 <- par_df[[2]] %>% group_by(Species) %>%
    map_kumara2(seq(0,1, length.out = 100), .) 
  
  data_plot <- out$sim_data$sim_data %>%
    mutate(Species = paste0("X",spp))
  
  ggplot() + 
    geom_line(data = par_set1, aes(x=x, y=y, group=as.character(draw)), alpha = 0.05, colour = "dodgerblue") +
    geom_line(data = par_set2, aes(x=x, y=y, group=as.character(draw)), alpha = 0.05, colour = "yellow") +
    geom_jitter(data = data_plot, mapping = aes(x = treat, y = y), 
                height = 0, width = 0.05, alpha = 0.5) +
    facet_wrap(~Species, scales = "free_y") +
    geom_smooth(data = data_plot, aes(x=treat, y=y), 
                se = F, colour = "green", alpha = 0.5, lwd = 0.5) +
    geom_line(data = true_tol, aes(x = x, y = y), 
              colour = "blue", lwd = 0.5, alpha = 0.5) +
    theme_minimal()
}


#produce latex table to summarize fit (using xtable)
sim_table <- function(out, file_name){
  options(xtable.sanitize.colnames.function=identity,
          xtable.sanitize.rownames.function=identity)
  
  sim_table <- out$sim_results$valid_df %>%
    data.frame() %>%
    set_colnames(c(
      "proportion simulated $>$ true parameter", 
      "Propoortion of true parameters \n within 95\\% credible interval")
    ) %>%
    set_rownames(c("$\\alpha$", "$\\beta$", 
                   "$\\zeta$", "$\\delta$", 
                   "$\\epsilon$", "$\\beta_0$", 
                   "$\\beta_1$", "$\\nu$")) %>%
    xtable()
  align( sim_table ) <- c( 'l', 'p{1.5in}', 'p{1.5in}' )
  
  if(missing(file_name)) print.xtable(x = sim_table) 
  else print.xtable(x = sim_table, file = file_name)
}



#fit idea data and challenging data
set.seed(123)
ideal_data <- sim_tolerance(
  n_spp_treat = 10, n_treat = 5, 
  n_spp = 14, random_axis = F, shift = 0
)

ideal_results <- stan_tolerance(ideal_data)
ideal_out <- list(sim_data=ideal_data, sim_results=ideal_results)

challenge_data <- sim_tolerance(
  n_spp_treat = 10, n_treat = 5, 
  n_spp = 14, random_axis = F, shift = 0,
  hyper_sds = list(
    a = 1, 
    b = 1, 
    c = 1, 
    d = 1,
    e = 1,
    beta_0 = 0.5,
    beta_1 = 0.5,
    nu = 1
  )
)



challenge_results <- stan_tolerance(challenge_data)
challenge_out <- list(sim_data=challenge_data, sim_results=challenge_results)

#sim_table(ideal_out)
#sim_table(challenge_out)
sim_table(ideal_out, "derived_files/idea_simtable.tex")
sim_table(challenge_out, "derived_files/idea_challengetable.tex")

ideal_out$sim_results$prop_zero
ideal_out$sim_results$corr_models
ideal_out$sim_results$diff_models

challenge_out$sim_results$prop_zero
challenge_out$sim_results$corr_models
challenge_out$sim_results$diff_models

cairo_pdf(filename = "figures/B1.pdf")
plot_sim_out(ideal_out)
dev.off()

cairo_pdf(filename = "figures/B2.pdf")
plot_sim_out(challenge_out)
dev.off()