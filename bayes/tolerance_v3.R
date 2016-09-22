#Run stan model for lasthenia tolerance curves

#load lasthenia data set
source("tolerance_functions.R")
emery <- load_emery()

#construct list to pass to stan
dataList <- list(N = nrow(emery), 
                 x = emery$treat, 
                 y = emery$Inflor_biomass,
                 numSpp = length(unique(emery$sppint)),
                 sppint = emery$sppint
)


#Run stan with default options:
  #4 chains
  #2000 iterations
  #save every sample
  #50% warmup
  #random initial parameter values
  #(saving 1000 draws per chain is PLENTY for this model in stan)
  #save samples and diagnostic data to files

#parameters to track
watch <- c("a", "b", "c","d", "e", "e1", "mu",
           "nu", "beta_0", "beta_1", 
           "var_c", "mean_c")
stan.fit <- stan(file = "bayes/tolerance_v3.stan",
                 data = dataList,
                 pars = watch,
                 sample_file = "bayes/samples/tolerance_v3.samples",
                 diagnostic_file = "bayes/samples/tolerance_v3.diag"
)

