all: 	figures/fig1.pdf \
	derived_files/lasth_100_post.R \
	bayes/samples/tolerance_v3_*.csv \
	figures/fig2.pdf \
	derived_files/curve_K.csv \
	figures/fig3.pdf \
	derived_files/pgls.tex \
	figures/A2.pdf

	
#make concept plot
figures/fig1.pdf: reactionNorm_conceptPlots2.R
	Rscript reactionNorm_conceptPlots2.R


#create 100 pruned, ultrametric posterior trees
derived_files/lasth_100_post.R: generate_post_trees.R data/C1.trees
	Rscript generate_post_trees.R


#run stan, and tidy files for downstream analysis and plotting
bayes/samples/tolerance_v3_*.csv: bayes/tolerance_v3_alt.stan bayes/tolerance_v3.R
	Rscript bayes/tolerance_v3.R

#posterior predictive checks plot and pvalue
figures/fig2.pdf: bayes/postpred.R bayes/samples/tolerance_v3_*.csv 
	Rscript bayes/postpred.R

#phylogenetic signal for curves and parameters
derived_files/curve_K.csv: phylo_signal.R bayes/samples/tolerance_v3_*.csv bayes/stan_par1_df.csv bayes/fitted_points_mod1.csv
	Rscript phylo_signal.R

figures/fig3.pdf: bayes/samples/tolerance_v3_*.csv derived_files/curve_K.csv tolerance_v3_plotting.Rmd tolerance_functions.R
	Rscript -e "rmarkdown::render('tolerance_v3_plotting.Rmd')"

derived_files/pgls.tex: moisture_pic.R bayes/samples/tolerance_v3_*.csv
	        Rscript moisture_pic.R

figures/A2.pdf: bayes/tolerance_v3_alt.stan tolerance_functions.R bayes/simulate_tolerance.R
	Rscript bayes/simulate_tolerance.R

