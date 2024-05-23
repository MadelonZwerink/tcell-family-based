folder <- "./results/model1/"
run_name <- "st1"
seed <- 4321
families = 500

bp_rule <- 'runif(1, min = 0, max = 3.5)'
dp_rule <- 0.4 #default
rq_rule <- 0.5 #default
t_start_dist <- 'rlnorm(nr_of_families, meanlog = 1.4, sdlog = 0.3)'
t_run_rule <- 'runif(1, min = 0.5, max = 4)'
nr_burst_divs <- 'sample(c(2,3,4), nr_of_families, replace = TRUE)'
quality_dist <- NULL
ASD <- FALSE
burst_time <- 0.15 #default
max_run_time <- 8 #default
min_t_start <- 2.5

source("./R/run_parameter_set.R")

