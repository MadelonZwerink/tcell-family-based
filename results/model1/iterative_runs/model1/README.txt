This folder contains all plots and statistics for the iterated run with 100 simulations, each consisting of 500 families.
Order of execution:
1) Run model1_multiple_runs_data_generation.R (done on cluster, output was several R objects containing all necessary information to get stats and figures)
2) Run model1_multiple_runs_plots_stats.R
3) Run model1_multiple_runs_panelplot_combi.R

Dependencies:
- speed_discrete_burst.R
- functions_multiple_sims.R