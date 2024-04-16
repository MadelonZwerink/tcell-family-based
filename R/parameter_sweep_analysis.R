setwd("C:/Users/madel/OneDrive/Documenten/BiBC/Major_internship/tcell-family-based/")

source("R/speed_discrete_burst.R")

parameter_sweep <- read.delim("data/processed/parameter_sweep_results.tsv", sep = "\t")
parameter_sweep <- parameter_sweep[order(parameter_sweep$total_SSR),]
parameter_sweep_exsmall <- parameter_sweep[order(parameter_sweep$total_SSR_exsmall),]

parameter_sweep_small <- read.delim("C:/Users/madel/OneDrive/Documenten/BiBC/Major_internship/Tcell_familybased/model_3_discrete_burst_PQ/output/parameter_sweep6_day78.tsv", sep = "\t")
parameter_sweep_small <- parameter_sweep_small[order(parameter_sweep_small$total_SSR),]

top_results <- 20
families <- 500
fam_counter <- 0
prim_parameters <- list()

for(r in seq(top_results)) {
  prim_parameters[[r]] <- pick_parameters(bp_rule = parameter_sweep$bp_dist[r],
                                          dp_rule = 0.3,
                                          rq_rule = parameter_sweep$rq_dist[r],
                                          nr_of_families = families,
                                          response_nr = 1,
                                          t_run_rule = parameter_sweep$t_run_dist[r],
                                          t_start_dist = parameter_sweep$t_start_dist[r],
                                          nr_burst_divs = parameter_sweep$b_div[r])
}
plots_famsize_dist <-  list()

for(r in seq(top_results)){
  plots_famsize_dist[[r]] <- plot_grid_famsize_dist(prim_parameters[[r]]) 
}

plots_famsize_dist

#-------------------------------------------------------------------------------

fam_counter <- 0
prim_parameters_exsmall <- list()

for(r in seq(top_results)) {
  prim_parameters_exsmall[[r]] <- pick_parameters(bp_rule = parameter_sweep$bp_dist[r],
                                          dp_rule = 0.3,
                                          rq_rule = parameter_sweep$rq_dist[r],
                                          nr_of_families = families,
                                          response_nr = 1,
                                          t_run_rule = parameter_sweep$t_run_dist[r],
                                          t_start_dist = parameter_sweep$t_start_dist[r],
                                          nr_burst_divs = parameter_sweep$b_div[r])
}
plots_famsize_dist_exsmall <-  list()

for(r in seq(top_results)){
  plots_famsize_dist_exsmall[[r]] <- plot_grid_famsize_dist(prim_parameters_exsmall[[r]]) 
}

plots_famsize_dist_exsmall

#-------------------------------------------------------------------------------

prim_parameters_best <- list()
families <- 1000
fam_counter <- 0
best_results_index <- c(1, 5, 8, 13)
top_results <- length(best_results_index)

for(r in seq(top_results)) {
  p <- best_results_index[r]
  prim_parameters_best[[r]] <- pick_parameters(bp_rule = parameter_sweep$bp_dist[p],
                                          dp_rule = 0.3,
                                          rq_rule = parameter_sweep$rq_dist[p],
                                          nr_of_families = families,
                                          response_nr = 1,
                                          t_run_rule = parameter_sweep$t_run_dist[p],
                                          t_start_dist = parameter_sweep$t_start_dist[p],
                                          nr_burst_divs = parameter_sweep$b_div[p])
}
plots_famsize_dist_best <-  list()

for(r in seq(top_results)){
  plots_famsize_dist_best[[r]] <- plot_grid_famsize_dist(prim_parameters_best[[r]]) 
}

plots_famsize_dist_best
