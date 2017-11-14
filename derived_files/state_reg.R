Genus_species <- c("debilis", "ferrisiae", "chrysantha", 
                   "fremontii", "coulteri", "microglossa",
                   "platycarpha", "conjugens", "gracilis", 
                   "minor", "glabrata", "burkei",
                   "californica", "glaberrima")

state_reg <- c("terrestrial", "vernal", "vernal",
               "vernal", "vernal", "terrestrial",
               "aqua_terr", "vernal", "aqua_terr",
               "terrestrial", "vernal", "vernal",
               "aqua_terr", "vernal")
names(state_reg) <- Genus_species

state_reg_aqua_terr2terr <- c("terrestrial", "vernal", "vernal",
                              "vernal", "vernal", "terrestrial",
                              "terrestrial", "vernal", "terrestrial",
                              "terrestrial", "vernal", "vernal",
                              "terrestrial", "vernal")
names(state_reg_aqua_terr2terr) <- Genus_species

state_reg_aqua_terr2vernal <- c("terrestrial", "vernal", "vernal",
                                "vernal", "vernal", "terrestrial",
                                "vernal", "vernal", "vernal",
                                "terrestrial", "vernal", "vernal",
                                "vernal", "vernal")
names(state_reg_aqua_terr2vernal) <- Genus_species


reg_df <- tibble(
  habit = state_reg, 
  aqua_terr2terr = state_reg_aqua_terr2terr, 
  aqua_terr2vernal = state_reg_aqua_terr2vernal, 
  Species = names(state_reg)
) 

