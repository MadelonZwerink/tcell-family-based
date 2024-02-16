parameter_sweep <- read.delim("C:/Users/madel/OneDrive/Documenten/BiBC/Major_internship/Tcell_familybased/model_3_discrete_burst_PQ/output/parameter_sweep6_day78.tsv", sep = "\t")
r = 1
families = 200

prim_parameters <- pick_parameters(bp_rule = parameter_sweep$bp_dist[r],
                                   dp_rule = "0.3",
                                   nr_of_families = families,
                                   response_nr = 1,
                                   t_run_rule = parameter_sweep$t_run_dist[r],
                                   t_start_dist = parameter_sweep$t_start_dist[r],
                                   nr_burst_divs = parameter_sweep$b_div[r],
                                   #quality_dist = parameter_sweep$quality_dist[r],
                                   rq_rule = parameter_sweep$rq_dist[r])

plot_famsize_dist(prim_parameters, timepoint = 8, method = "time", show_title = F)

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
