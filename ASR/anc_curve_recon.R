# This R script demonstrates the equivalence of ancestral curve reconstruction using 
# time warping alignment with ancestral curve reconstruction using inverse functions
#
# NOTE: generalized time warping is not currently implemented in R, but pairwise 
# alignment is implemented in the dtw package. DTW allows ancestral state
# reconstruction of any function using the fast PIC method (as in fastAnc),
# but does not align landmarks across species for distance-based multivariate
# methods such as K-mult (Adams 2014a) or D-PGLS (Adams 2014b).
# To perform alignment across species, generalized time warping may be performed
# using the current GTW MATLAB implementation: http://www.f-zhou.com/ta_code.html

require(scales)
require(GPfit)
require(dtw)
#setwd("")
source("functions.R")

# simulates a dataset with 4 species
data <- f(N = 4)
tree <- data$tree

# simulated environmental gradient (level of water deficiency for plants)
x <- data$wilting$observed$water_deficiency

# simulated proportion wilting at each level of water deficiency
Y <- t(apply(data$wilting$observed$tip_coefficients,1,function(X) logit_fx(X,x)))

# ancestral curve reconstruction of x and y axes via dynamic time warping alignment
# if gpr_fit==TRUE, Gaussian process regression is used (NOT TO BE CONFUSED WITH PGPR)
#   this option allows any shaped function to be reconstructed along the x and y axes
# if gpr_fit==FALSE, glm estimation is done using link=logit (MUCH faster, than GPR
#   but limited to datasets where y ranges from 0 to 1)
dtw_root_vals <- fast_anc_hand(x,Y,tree,root_only = TRUE,gpr_fit = FALSE)

# perform PGLS-based ancestral curve reconstruction (keep root) and evaluate function with estimated coefficients
inv_root_coefs <- pgls_curve(data$tree,data$wilting$true$tip_coefficients)[1,]
inv_root_vals <- logit_fx(inv_root_coefs,data$wilting$observed$water_deficiency)

# results from both methods are essentially identical
plot(inv_root_vals,dtw_root_vals)
summary(lm(dtw_root_vals~inv_root_vals))
