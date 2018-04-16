source("tolerance_functions.R")
emery <- load_emery()
tolerance_spp <- emery$Species %>% unique
map <- purrr::map


#key between tree tip names and taxa names used in expriment
taxa_key <- read_delim("data/Isolate_Key.txt", delim = "\t", col_names = F) %>%
  filter(!is.na(X3)) %>%
  mutate(
    X3 = ifelse(X3 == "amb", 'ambl', X3),
    #X3 = paste0("'", X3, "'"), #these quotes are bad
    names = ifelse( is.na(X8),X6,X8)
    )


n_trees <- 40273 #total number of trees
n_samp <- 1000 #how many trees to read in
n_final <- 100 #how many trees to save
#the index for trees after thinning (avoid autocorrelation)
thin_samp <- seq(1, n_samp, length.out = n_final) %>% round(0)

#returns `n_final` trees -- rooted, pruned, ultrametricized
trees_post <- read_lines("data/C1.trees", n_max = n_samp, skip = n_trees-n_samp) %>%
  .[thin_samp] %>%
  map(function(x){
    tree_i <- read.tree(text = x)
    #tree_i <- read.tree(text = trees_post[1])
    tree_taxa_switch <- tree_i$tip.label %>%
      map_df(~ taxa_key %>% 
               filter(X3 == .x) %>% 
               select(X3, names)) %>%
      group_by(names) %>%
      sample_n(1)
    
    new_tree <- tree_i %>% 
      drop.tip(.$tip.label[-match(tree_taxa_switch$X3, .$tip.label)])
    
    #new_tree$tip.label_alt <- new_tree$tip.label %>% match(taxa_key$X3) %>% taxa_key$names[.]
    new_tree$tip.label <- new_tree$tip.label %>% match(taxa_key$X3) %>% taxa_key$names[.]
  
    drop <- new_tree$tip.label[!new_tree$tip.label %in% tolerance_spp]
    
    new_tree %>% root(c("congdonii", "pusillus")) %>%
      chronos(lambda = 1) %>%
      drop.tip(tip = drop)
    
    })

class(trees_post) <- "multiPhylo"

dump(list = "trees_post", file = "derived_files/lasth_100_post.R")

#dist.topo(trees_post) %>% mean
#dist.topo(trees_post, method = "PH") %>% mean

#par(mfrow = c(3,3))
#1:9 %>% map(~{
#  par(mar=c(0,0,0,0))
#  plot.phylo(trees_post[[.x]])
#  })

#sim_tree <- rtree(n = 14, br = trees_post[[1]]$edge.length)
#sim_tree$tip.label <- trees_post[[1]]$tip.label
#dist.topo(trees_post[[1]], sim_tree, method = "PH85")
#plot(sim_tree$edge.length, trees_post[[1]]$edge.length)

#dist.topo(trees_post[[1]], sim_tree)
#dist.topo(trees_post) %>% mean
#dist.topo(trees_post, method = "PH") %>% mean
