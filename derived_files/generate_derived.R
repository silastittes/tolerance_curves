#GENERATE DERIVED PARAMETERS AND ANALYSES FOR LASTHENIA TOLERANCE CURVE DATA

source("load_data.R")


###############################################################
##OPTIMA-------------------------------------------------------
###############################################################

#Approximate location of water treatment optima
xseq_big <- seq(0,1, length.out=10000)
draw_list <- as.list(1:(ndraws*1))
sppint_list <- as.list(1:max(emery$sppint))
opt <- mclapply(sppint_list, function(sp){
  sapply(draw_list, function(pr){
    opt.kumara(xs = xseq_big,
               a = posts$a[pr,sp],
               b = posts$b[pr,sp],
               c = posts$c[pr,sp],
               d = posts$d[pr,sp],
               e1 = posts$e1[pr,sp]
    )
  }
  )
}
)

names(opt) <- unique(emery$Species)
maximadf <- data.frame(opt)
write.table(x = maximadf, file = "derived_files/maxima_draws.txt")
maximadf <- read.table("derived_files/maxima_draws.txt")


#SANITY CHECK
#mden <- apply(maximadf, 2, density)
#ylimit <- max(sapply(mden, function(x) max(x$y)))
#xlimit <- range(sapply(mden, function(x) range(x$x)))
#plot(NA,NA, xlim=xlimit, ylim=c(0,ylimit), xlab="", ylab="", main="")
#sapply(mden, FUN = function(x) lines(x = x$x, y = x$y) )


###############################################################
##INTEGRAL-----------------------------------------------------
###############################################################

#Approximate location of water treatment optima
xseq_big <- seq(0,1, length.out=10000)
draw_list <- as.list(1:(ndraws*1))
sppint_list <- as.list(1:max(emery$sppint))
int <- mclapply(sppint_list, function(sp){
  sapply(draw_list, function(pr){
    int_kumara(xs = xseq_big,
               a = posts$a[pr,sp],
               b = posts$b[pr,sp],
               c = posts$c[pr,sp],
               d = posts$d[pr,sp],
               e1 = posts$e1[pr,sp])
  })
})

names(int) <- unique(emery$Species)
integraldf <- data.frame(int)
write.table(x = integraldf, file = "derived_files/integral_draws.txt")
integraldf <- read.table("derived_files/integral_draws.txt")


#SANITY CHECK
#mden <- apply(integraldf, 2, density)
#ylimit <- max(sapply(mden, function(x) max(x$y)))
#xlimit <- range(sapply(mden, function(x) range(x$x)))
#plot(NA,NA, xlim=xlimit, ylim=c(0,ylimit), xlab="", ylab="", main="")
#sapply(mden, FUN = function(x) lines(x = x$x, y = x$y) )


###############################################################
##CURVE SIGNAL-------------------------------------------------
###############################################################

nspp <- nrow(unique(emery[,c(2,7)]))
xseq_comm <- seq(1,5, length.out=200)

sppMaxVal <- by(data = emery$Inflor_biomass, INDICES = emery$sppint, FUN = max)
names(sppMaxVal) <- unique(emery$Species)

draw_list <- as.list(1:ndraws) #posterior draws as list
sppint_list <- as.list(1:max(emery$sppint)) #species are columns
curveK_draws <- sapply(draw_list, 
                       function(pr){ 
                         spp_mat <- t(sapply(sppint_list, function(sp){
                           scale.kumara(xs = xseq_comm,
                                        a = posts$a[pr,sp],
                                        b = posts$b[pr,sp],
                                        c = posts$c[pr,sp],
                                        d = posts$d[pr,sp],
                                        e1 = posts$e1[pr,sp]
                                        )
                           })
                           )
                         rownames(spp_mat) <- unique(emery[,2])
                         physignal(A = spp_mat, phy = lasth, 
                                   iter = 1)$phy.signal
                       })



curveK_draws_scaled <- sapply(draw_list, function(pr){
  spp_mat <- t(sapply(sppint_list, function(sp){
    scale.kumara(xs = xseq_comm,
                 a = posts$a[pr,sp],
                 b = posts$b[pr,sp],
                 c = posts$c[pr,sp]/sppMaxVal[sp], #alt
                 d = posts$d[pr,sp],
                 e1 = posts$e1[pr,sp]
                 )
    })
  )
  rownames(spp_mat) <- unique(emery[,2])
  physignal(A = spp_mat, phy = lasth, iter = 1)$phy.signal
})

write.table(x = curveK_draws_scaled, file = "derived_files/curve_K_scaled.txt")
write.table(x = curveK_draws, file = "derived_files/curve_K.txt")


##############################################################
##OUwie Draws-------------------------------------------------
##############################################################


#avoid the annoying hit enter for next plot thing
par(ask=F)

#Set up data for OUwie
Genus_species <- as.character(unique(emery$Species))
Reg <- c("terrestrial", "vernal", "vernal",
         "vernal", "vernal", "vernal",
         "aqua_terr", "vernal", "vernal",
         "aqua_terr", "terrestrial", "aqua_terr",
         "vernal", "terrestrial"
)


state_reg <- Reg
names(state_reg) <- as.character(Genus_species)


#ml est for habitat reconstruction
regFit_ml <- rerootingMethod(lasth, state_reg, model = "ER", tips = T)
#stochastic mapping for anc. state. recon. of habitat

#number of stochastic maps of habitat
mapsims <- 100

#generate maps
reg_fit <- make.simmap(tree = lasth, state_reg, model = "ER", nsim = mapsims)


#analyze transitions
#tt <- 8
#colnames(countSimmap(reg_fit)$Tr)[tt]
#plot(table(countSimmap(reg_fit)$Tr[,tt]))
#quantile(countSimmap(reg_fit)$Tr[,tt], c(0.225, 0.975))

#prep for OUwie
states <- colnames(reg_fit[[1]]$mapped.edge)
names(states) <- c("red", "green", "blue")
colz <- names(states)
names(colz) <- states

#example map
#plotSimmap(reg_fit[[1]], colors = colz)
#reg_fit[[1]]$mapped.edge
#nodelabels()
#tiplabels()
#dev.off() #reset plotting options


#assign node habitat state as maximum stochastic map node state
simstates <- lapply(reg_fit, function(x){
  states[apply(x$mapped.edge, 1, which.max)]
  #apply(x$mapped.edge, 1, function(z){
  #sample(x = states, size = 1, prob = z, replace = F)})
})

#assign habitat to ml nodes
states_ml <- apply(regFit_ml$marginal.anc, 1, function(x) states[which.max(x)])

#subset to nodes only
anc_sims <- lapply(simstates, function(x) x[lasth$edge[,2] > 
                                              length(lasth$tip.label)])
hab_sims <- lapply(simstates, function(x) x)

#function to apply OUwie over posterior tolerance draws, 
#and stochastic map draws
#OUwie data frame = taxa, regime, trait value
ouwie_draws <- function(draws, df, mod, tree, sims){
  
  #process across stochastic map sims
  mapdraw <- lapply(sims, function(y){
    #combine true tips states to stochastic map anc states
    anc_states <- y[tree$edge[,2] >= length(tree$tip.label)]
    tree$node.label <- anc_states
    tip_states <- y[tree$edge[,2] <= length(tree$tip.label)]
    
    #run posterior draws against each stochastic habitat reconstruction
    draws <- mclapply(X = as.list(1:draws), FUN = function(x){
      X <- matrix(df[x,])
      #X <- apply(maximadf, 2, mean)
      ouwieDat <- data.frame(Genus_species=Genus_species,
                             Reg=tip_states, X=unlist(X))
      ouwieOUMVA <- OUwie(tree, data = ouwieDat,
                          model = mod, root.station = T)
      #message(x)
    })
  })
  return(mapdraw)
}

ouwie_draws_ml <- function(draws, df, mod, tree, anc){
  
  anc_states <- anc[tree$edge[,2] > length(tree$tip.label)]
  tree$node.label <- anc_states
  tip_states <- anc[tree$edge[,2] <= length(tree$tip.label)]
  
  #run posterior draws against each stochastic habitat reconstruction
  draws <- mclapply(X = as.list(1:draws), FUN = function(x){
    X <- matrix(df[x,])
    #X <- apply(maximadf, 2, mean)
    ouwieDat <- data.frame(Genus_species=Genus_species,
                           Reg=tip_states, X=unlist(X))
    ouwieOUMVA <- OUwie(tree, data = ouwieDat,
                        model = mod, root.station = T)
  })
  return(draws)
}

#set up parameter list
post_params <- list(integ = integraldf, optma = maximadf,
                    c = posts$c, d = posts$d, e1 = posts$e1)
par_draws <- 100

### WITH maximum likelihood

#OU1 ML
for(x in names(post_params)){
  xdraw <- ouwie_draws_ml(df = post_params[[x]], draws = par_draws, 
                          mod = "OU1", tree = lasth,
                          states_ml)
  dump(list = "xdraw", 
       file = paste0("derived_files/", x, "_OU1_ml.txt"))
}

#OUM ML
for(x in names(post_params)){
  xdraw <- ouwie_draws_ml(df = post_params[[x]], draws = par_draws, 
                 mod = "OUM", tree = lasth,
                 states_ml)
  dump(list = "xdraw", 
       file = paste0("derived_files/", x, "_OUM_ml.txt"))
}


#with stochastic mapping

#OU1
for(x in names(post_params)){
  xdraw <- ouwie_draws(df = post_params[[x]], draws = par_draws,
                       mod = "OU1", tree = lasth,
                       sims = hab_sims)
  dump(list = "xdraw", 
       file = paste0("derived_files/", x, "_OU1_simmap.txt"))
}

for(x in names(post_params)){
  xdraw <- ouwie_draws(df = post_params[[x]], draws = par_draws, 
                       mod = "OUM", tree = lasth,
                       sims = hab_sims)
  dump(list = "xdraw", 
       file = paste0("derived_files/", x, "_OUM_simmap.txt"))
}



#ou1_simmap_all_par <- mclapply(post_params, function(x){
#  ouwie_draws( df = x, draws = par_draws,
#               mod = "OU1", tree = lasth,
#               sims = hab_sims)
#  })


#OUM
#oum_simmap_all_par <- mclapply(post_params, function(x){
#  ouwie_draws( df = x, draws = par_draws,
#               mod = "OUM", tree = lasth,
#               sims = hab_sims)
#})


#write data to files -- perserve format using dump()

#dump(c("ou1_ml_all_par", "oum_ml_all_par"), 
#     file = "derived_files/ouwie_ml_all_par.R")
#dump(c("ou1_simmap_all_par", "oum_simmap_all_par"), 
#     file = "derived_files/ouwie_simmap_all_par.R")

#source("derived_files/ouwie_ml_all_par.R")
#source("derived_files/ouwie_simmap_all_par.R")


compMat <- sapply(as.list(1:14), function(x){
  sapply(as.list(1:14), function(y){
    #val <- mean(posts$c[,y] > posts$c[,x])
    val <- mean(integraldf[,y] > integraldf[,x])
    #ifelse(  y > x & (val >= 0.95 | val <= 0.05), yes = "*", no = "-")
    if(y > x){
      if(val >= 0.95){
        ">"
      } else if(val <= 0.05){
        "<"
      } else{
        "-"
      }
    }
  })
})

length(grep(">", compMat)) + length(grep("<", compMat))
#c = 57

rownames(compMat) <- unique(emery$Species)
colnames(compMat) <- unique(emery$Species)
compMat

mean(maximadf[,3] > maximadf[,4])

rownames(compMat) <- unique(emery$Species)
colnames(compMat) <- unique(emery$Species)
compMat


#quantile(sapply(ou1_simmap_all_par$c, function(y){
#  sapply(y, function(x) x$AICc)}))

#quantile(sapply(oum_simmap_all_par$c, function(y){
#  sapply(y, function(x) x$AICc)}))

