Place all files in the same folder. Be sure to install the following packages first:
	ape, phytools, geomorph, geiger, GPfit, dtw

To perform ancestral curve reconstruction of a logit function:
-source("functions.R")
-Fit coefficients for each tip using glm(y~x,family=quasibinomial("logit")) where y ranges from 0 to 1
---Or, if dealing with count data: glm(cbind(y,n-y)~x,family=binomial("logit")) where y = # successes out of n trials
-Put tip coefficients for each species in a matrix with Nspecies rows (name rows appropriately and 2 columns (b and m)
-Next, run the following code: anc_nodes <- pgls_curve(tree,tip_coefficients)
-If vals_only=TRUE, pgls_curve returns reconstructed y values which can be regressed against x; otherwise, reconstructed logit coefficients are returned

To perform ancestral curve reconstruction on a function of ANY shape (note, this function is VERY slow):
-Load ape, dtw, and GPfit packages
-source("functions.R")
-Fit functions and evaluate along a sequence (x) of x-values
-Store the resulting y-values in a matrix (Y) with Nspecies rows (name rows appropriately), and each column corresponding to an x-value evaluation
-Next, run the following code: anc_Y <- fast_anc_hand(x,Y,tree,TRUE,TRUE)
-Each row of anc_Y can be regressed against x to estimate ancestral node functions

The equivalence of inverse function and time warping alignment methods for ancestral curve reconstruction are demonstrated in anc_curve_recon.R

To reproduce simulations in paper:
-For ancestral curve simulations, also install the following packages: minqa, inline, caTools, MASS
-NOTE - must source ICA_functions_SHOWCASE.R and cxxFunctions_OU_Generic_SHOWCASE.cxx from https://github.com/fpgpr/DemoScripts/
-NOTE - the simulations take several days (each) to run
-Power simulations for Kmult and D-PGLS can be performed by sourcing phylo_signal_power.R and dpgls_power.R, respectively.
-Ancestral curve reconstruction simulations can be performed with anc_curve_recon_simulations.R