source("./R/speed_discrete_burst.R")

families <- 500
fam_counter <- 0

prim_parameters <- pick_parameters(bp_rule = "runif(1, min = 0.5, max = 3)",
                                   d_p = 0.3,
                                   nr_of_families = families,
                                   response_nr = 1,
                                   t_run = "rnorm(1, mean = 3, sd = 0.2)",
                                   t_start = "rlnorm(families, meanlog = 1.3, sdlog = 0.15)",
                                   nr_burst_divs = "sample(c(1,2,3,4),families, replace = TRUE)")

famsize_dist_d5 <- plot_famsize_dist(prim_parameters, timepoint = 5, method = "time", binwidth = 1, show_title = F)
famsize_dist_d6 <- plot_famsize_dist(prim_parameters, timepoint = 6, method = "time", binwidth = 1, show_title = F)
famsize_dist_d7 <- plot_famsize_dist(prim_parameters, timepoint = 7, method = "time", binwidth = 1, show_title = F)
famsize_dist_d8 <- plot_famsize_dist(prim_parameters, timepoint = 8, method = "time", binwidth = 1, show_title = F)

plot_grid(famsize_dist_d5[[1]], famsize_dist_d6[[1]], famsize_dist_d7[[1]],  famsize_dist_d8[[1]], nrow = 2, ncol = 2)
