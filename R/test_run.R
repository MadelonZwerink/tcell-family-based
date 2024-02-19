parameter_sweep <- read.delim("C:/Users/madel/OneDrive/Documenten/BiBC/Major_internship/Tcell_familybased/model_3_discrete_burst_PQ/output/parameter_sweep6_day78.tsv", sep = "\t")
parameter_sweep <- parameter_sweep[order(parameter_sweep$mean_SSR_day7), ]
r = 3
families = 1000

prim_parameters <- pick_parameters(bp_rule = "runif(1, min = 0.5, max = 2.5)",
                                   dp_rule = "0.3",
                                   nr_of_families = families,
                                   response_nr = 1,
                                   t_run_rule = "rnorm(1, mean = 2.5, sd = 0.7)",
                                   t_start_dist = "rlnorm(families, meanlog = 1.35, sdlog = 0.2)",
                                   nr_burst_divs = "sample(c(2,3,4), nr_of_families, replace = TRUE)",
                                   #quality_dist = parameter_sweep$quality_dist[r],
                                   rq_rule = "rbeta(1, shape1 = 2, shape2 = 1.6)")

plot_grid_famsize_dist(prim_parameters)

sec_parameters <- pick_parameters(prev_parameters = prim_parameters,
                                        bp_rule = parameter_sweep$bp_dist[r],
                                        dp_rule = "0.3",
                                        response_nr = 2,
                                        t_run_rule = parameter_sweep$t_run_dist[r],
                                        t_start_dist = "rlnorm(nr_of_families, meanlog = 1.25, sdlog = 0.2)",
                                        nr_burst_divs = "sample(c(1,2,3,4), nr_of_families, replace = TRUE)", 
                                        #quality_dist = parameter_sweep$quality_dist[result],
                                        rq_rule = parameter_sweep$rq_dist[r])

ter_parameters <- pick_parameters(prev_parameters = sec_parameters,
                                        bp_rule = parameter_sweep$bp_dist[r],
                                        dp_rule = "0.3",
                                        response_nr = 3,
                                        t_run_rule = parameter_sweep$t_run_dist[r],
                                        t_start_dist = "rlnorm(nr_of_families, meanlog = 1.25, sdlog = 0.2)",
                                        nr_burst_divs = "sample(c(1,2,3,4), nr_of_families, replace = TRUE)", 
                                        #quality_dist = parameter_sweep$quality_dist[result],
                                        rq_rule = parameter_sweep$rq_dist[r])
