all: 	analyses_and_viz/conceptual_ReNorm.pdf \
	derived_files/lasth_100_post.R \
	bayes/samples/tolerance_v3_*.csv \
	analyses_and_viz/postpred.pdf \
	derived_files/curve_K.csv \
	figures/*pdf
	
#make concept plot
analyses_and_viz/conceptual_ReNorm.pdf: analyses_and_viz/reactionNorm_conceptPlots2.R
	Rscript analyses_and_viz/reactionNorm_conceptPlots2.R


#create 100 pruned, ultrametric posterior trees
derived_files/lasth_100_post.R: generate_post_trees.R data/C1.trees
	Rscript generate_post_trees.R


#run stan, and tidy files for downstream analysis and plotting
bayes/samples/tolerance_v3_*.csv: bayes/tolerance_v3_alt.stan bayes/tolerance_v3.R
	Rscript bayes/tolerance_v3.R

#posterior predictive checks plot
analyses_and_viz/postpred.pdf: bayes/postpred.R bayes/samples/tolerance_v3_*.csv 
	Rscript bayes/postpred.R

#phylogenetic signal for curves and parameters
derived_files/curve_K.csv: phylo_signal.R bayes/samples/tolerance_v3_*.csv bayes/stan_par1_df.csv bayes/fitted_points_mod1.csv
	Rscript phylo_signal.R

figures/*pdf: bayes/samples/tolerance_v3_*.csv analyses_and_viz/tolerance_v3_plotting.Rmd
	Rscript -e "rmarkdown::render('analyses_and_viz/tolerance_v3_plotting.Rmd')"
