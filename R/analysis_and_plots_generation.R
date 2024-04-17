
families = 10000
set.seed(11)
prim_parameters <- pick_parameters(bp_rule = 'runif(1, min = 0.5, max = 3.5)',
                                   dp_rule = 0.4,
                                   rq_rule = 0.5,
                                   nr_of_families = families,
                                   response_nr = 1,
                                   t_run_rule = 'runif(1, min = 0.5, max = 4)',
                                   t_start_dist = 'rlnorm(nr_of_families, meanlog = 1.4, sdlog = 0.2)',
                                   nr_burst_divs = 'sample(c(1,2,3,4), nr_of_families, replace = TRUE)')

plot_famsize_dist(prim_parameters, timepoint = 8, method = "time")

plot_grid_famsize_dist(prim_parameters)

plot_Q_famsize(prim_parameters)

prim_resp <- plot_response(prim_parameters, 0, 40, 0.1)
prim_resp

sec_parameters <- pick_parameters(bp_rule = 'runif(1, min = 0.5, max = 3.5)',
                                  dp_rule = 0.4,
                                  rq_rule = 0.5,
                                  prev_parameters = prim_parameters,
                                  response_nr = 2,
                                  t_run_rule = 'runif(1, min = 0.5, max = 4)',
                                  t_start_dist = 'rlnorm(nr_of_families, meanlog = 1.4, sdlog = 0.2)',
                                  nr_burst_divs = 'sample(c(1,2,3,4), nr_of_families, replace = TRUE)')

ter_parameters <- pick_parameters(bp_rule = 'runif(1, min = 0.5, max = 3.5)',
                                  dp_rule = 0.4,
                                  rq_rule = 0.5,
                                  prev_parameters = sec_parameters,
                                  response_nr = 3,
                                  t_run_rule = 'runif(1, min = 0.5, max = 4)',
                                  t_start_dist = 'rlnorm(nr_of_families, meanlog = 1.4, sdlog = 0.2)',
                                  nr_burst_divs = 'sample(c(1,2,3,4), nr_of_families, replace = TRUE)')

generate_max_fam_plot(prim_parameters = prim_parameters, 
                      sec_parameters = sec_parameters, 
                      ter_parameters = ter_parameters,
                      type = "time",
                      time = 7,
                      relative = F,
                      show_title = T)


#-------------------------------------------------------------------------------


