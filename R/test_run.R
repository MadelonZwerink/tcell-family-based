parameter_sweep <- read.delim("C:/Users/madel/OneDrive/Documenten/BiBC/Major_internship/Tcell_familybased/model_3_discrete_burst_PQ/output/parameter_sweep6_day78.tsv", sep = "\t")
parameter_sweep <- parameter_sweep[order(parameter_sweep$total_SSR), ]
r = 3
families = 500

prim_parameters <- pick_parameters(bp_rule = "runif(1, min = quality[i], max = 2 + quality[i])",
                                   dp_rule = "0.3",
                                   nr_of_families = families,
                                   response_nr = 1,
                                   t_run_rule = "0.5 + (3 * quality[i])",
                                   t_start_dist = "rlnorm(families, meanlog = 1.35, sdlog = 0.2)",
                                   nr_burst_divs = "sample(c(2,3,4), nr_of_families, replace = TRUE)",
                                   quality_dist = "runif(nr_of_families, min = 0, max = 1)",
                                   rq_rule = "runif(1, min = 0, max = 1)")

plot_grid_famsize_dist(prim_parameters)

sec_parameters <- pick_parameters(prev_parameters = prim_parameters,
                                        bp_rule = "runif(1, min = 1.5, max = 2.5)",
                                        dp_rule = "0.3",
                                        response_nr = 2,
                                        t_run_rule = "runif(1, min = 1, max = 2)",
                                        t_start_dist = "rlnorm(nr_of_families, meanlog = 1.25, sdlog = 0.2)",
                                        nr_burst_divs = "sample(c(2,3), nr_of_families, replace = TRUE)", 
                                        #quality_dist = parameter_sweep$quality_dist[result],
                                        rq_rule = "runif(1, min = 0, max = 1)")
                                        #ASD = T)

ter_parameters <- pick_parameters(prev_parameters = sec_parameters,
                                        bp_rule = "runif(1, min = 1.5, max = 2.5)",
                                        dp_rule = "0.3",
                                        response_nr = 3,
                                        t_run_rule = "runif(1, min = 1, max = 2)",
                                        t_start_dist = "rlnorm(nr_of_families, meanlog = 1.25, sdlog = 0.2)",
                                        nr_burst_divs = "sample(c(2,3), nr_of_families, replace = TRUE)", 
                                        #quality_dist = parameter_sweep$quality_dist[result],
                                        rq_rule = "runif(1, min = 0, max = 1)")
                                        # ASD = T)
generate_max_fam_plot(prim_parameters, sec_parameters, ter_parameters, type = "max")
