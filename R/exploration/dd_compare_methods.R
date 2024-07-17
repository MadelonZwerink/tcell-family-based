# This script explores the different methods to execute model 2, one of which is 
# to fixate bp, and the other takes the sqrt of the DD to calculate prolif rate
# and run time, these methods are compared to model 1 in terms of the relationship
# betwween DD and the time it takes to reach this maximum

bp_rule <- 1.75
dp_rule <- 0.5 #default
rq_rule <- 0.5 #default
t_start_dist <- 'rlnorm(nr_of_families, meanlog = 1.4, sdlog = 0.3)'
t_run_rule <- '(q*12.6/1.75) + (nr_burst_divs[i]*t_burst[i])'
nr_burst_divs <- 'sample(c(2,3,4), nr_of_families, replace = TRUE)'
quality_dist <- "rbeta(nr_of_families, shape1 = 2.29, shape2 = 3.80)"
quality_noise <- FALSE
uniform_fam <- TRUE
ASD <- FALSE
burst_time <- 0.15 #default
max_run_time <- NULL
min_t_start <- 0
resolution_resp_plot <- 1

prim_parameters <- pick_parameters(bp_rule = bp_rule,
                                   dp_rule = dp_rule,
                                   rq_rule = rq_rule,
                                   t_start_dist = t_start_dist,
                                   t_run_rule = t_run_rule,
                                   nr_of_families = nr_of_families,
                                   nr_burst_divs = nr_burst_divs,
                                   response_nr = 1,
                                   quality_dist = quality_dist,
                                   quality_noise = quality_noise,
                                   q_noise_dist = q_noise_dist,
                                   uniform_fam = uniform_fam,
                                   ASD = ASD,
                                   burst_time = burst_time,
                                   max_run_time = max_run_time,
                                   min_t_start = min_t_start)

params <- prim_parameters[prim_parameters$cell_type == "P",]
ggplot(params, aes(x = 12.6*quality, y = t_run + t_start)) +
  geom_point() +
  labs(title = "Prolif. rate fixed on 1.75",
       x = "DD", y = "Time when DD is reached")+
  scale_y_continuous(limits = c(0, 15)) +
  scale_x_continuous(limits = c(0, 12))

mean(params$t_run + params$t_start)

#-------------------------------------------------------------------------------

bp_rule <- 'sqrt(q*12.6)'
dp_rule <- 0.5 #default
rq_rule <- 0.5 #default
t_start_dist <- 'rlnorm(nr_of_families, meanlog = 1.4, sdlog = 0.3)'
t_run_rule <- 'sqrt(q*12.6) + (nr_burst_divs[i]*t_burst[i])'
nr_burst_divs <- 'sample(c(2,3,4), nr_of_families, replace = TRUE)'
quality_dist <- "rbeta(nr_of_families, shape1 = 2.29, shape2 = 3.80)"
quality_noise <- FALSE
uniform_fam <- TRUE
ASD <- FALSE
burst_time <- 0.15 #default
max_run_time <- NULL
min_t_start <- 0
resolution_resp_plot <- 1

prim_parameters <- pick_parameters(bp_rule = bp_rule,
                                   dp_rule = dp_rule,
                                   rq_rule = rq_rule,
                                   t_start_dist = t_start_dist,
                                   t_run_rule = t_run_rule,
                                   nr_of_families = nr_of_families,
                                   nr_burst_divs = nr_burst_divs,
                                   response_nr = 1,
                                   quality_dist = quality_dist,
                                   quality_noise = quality_noise,
                                   q_noise_dist = q_noise_dist,
                                   uniform_fam = uniform_fam,
                                   ASD = ASD,
                                   burst_time = burst_time,
                                   max_run_time = max_run_time,
                                   min_t_start = min_t_start)

params <- prim_parameters[prim_parameters$cell_type == "P",]
ggplot(params, aes(x = 12.6*quality, y = t_run + t_start)) +
  geom_point() +
  labs(title = "Prolif. rate and time are sqrt of DD",
       x = "DD", y = "Time when DD is reached")+
  scale_y_continuous(limits = c(0, 15)) +
  scale_x_continuous(limits = c(0, 12))

mean(params$t_run + params$t_start)

#-------------------------------------------------------------------------------

bp_rule <- 'runif(1, min = 0, max = 3.5)'
dp_rule <- 0.5 #default
rq_rule <- 0.5 #default
t_start_dist <- 'rlnorm(nr_of_families, meanlog = 1.4, sdlog = 0.3)'
t_run_rule <- 'runif(1, min = 0.6, max = 4)' 
# was 0.5 in st2_nomax, but 0.6 ensures that all families can at least to their 4 burst 
# divisions, because 0.15*4 = 0.6
nr_burst_divs <- 'sample(c(2,3,4), nr_of_families, replace = TRUE)'
quality_dist <- NULL
ASD <- FALSE
burst_time <- 0.15 #default
max_run_time <- NULL
min_t_start <- 2.5

prim_parameters <- pick_parameters(bp_rule = bp_rule,
                                   dp_rule = dp_rule,
                                   rq_rule = rq_rule,
                                   t_start_dist = t_start_dist,
                                   t_run_rule = t_run_rule,
                                   nr_of_families = nr_of_families,
                                   nr_burst_divs = nr_burst_divs,
                                   response_nr = 1,
                                   quality_dist = quality_dist,
                                   quality_noise = quality_noise,
                                   q_noise_dist = q_noise_dist,
                                   uniform_fam = uniform_fam,
                                   ASD = ASD,
                                   burst_time = burst_time,
                                   max_run_time = max_run_time,
                                   min_t_start = min_t_start)

params <- prim_parameters[prim_parameters$cell_type == "P",]
params$dd <- params$bp * (params$t_run - params$nr_burst_divs * params$t_burst)
ggplot(params, aes(x = dd, y = t_run + t_start)) +
  geom_point() +
  labs(title = "Model 1",
       x = "DD", y = "Time when DD is reached") +
  scale_y_continuous(limits = c(0, 15)) +
  scale_x_continuous(limits = c(0, 12))

mean(params$t_run + params$t_start)

