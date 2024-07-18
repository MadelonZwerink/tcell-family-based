folder <- "./results/model1/"
run_name <- "st6"

rm(bp_rule,
   dp_rule,
   rq_rule,
   t_start_dist,
   t_run_rule,
   nr_of_families,
   nr_burst_divs,
   quality_dist,
   quality_noise,
   q_noise_dist,
   uniform_fam,
   ASD,
   burst_time,
   max_run_time,
   min_t_start)

seed <- 4321
nr_of_families = 500

bp_rule <- 'runif(1, min = 0, max = 3.5)'
dp_rule <- 0.5 #default
rq_rule <- 0.5 #default
t_start_dist <- '1 + rlnorm(nr_of_families, meanlog = 1.2, sdlog = 0.3)'
t_run_rule <- 'runif(1, min = 0.6, max = 4)' 
# was 0.5 in st2_nomax, but 0.6 ensures that all families can at least to their 4 burst 
# divisions, because 0.15*4 = 0.6
nr_burst_divs <- 'sample(c(2,3,4), nr_of_families, replace = TRUE)'
quality_dist <- NULL
ASD <- FALSE
burst_time <- 0.15 #default
max_run_time <- NULL
min_t_start <- 0
resolution_resp_plot <- 1
uniform_fam <- FALSE

source("./R/run_parameter_set.R")

