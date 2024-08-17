#!/usr/bin/env Rscript

generate_dd_branch <- function(mean_nr, cells) {
  target_sum <- cells * exp(mean_nr)
  
  # Randomly distribute the target sum across the cells
  random_distribution <- runif(cells)
  random_distribution <- random_distribution / sum(random_distribution) * target_sum
  
  # Calculate the all_nr values from the random distribution
  all_nr <- log(random_distribution)
  
  # Ensure values are within the desired range [0, 12]
  all_nr <- pmax(pmin(all_nr, 12), 0)
  
  return(round(all_nr, digits = 2))
}

if(Sys.info()[[4]]=="LAPTOP-3RJSLMKV") {
  print("Working from R-project")
  setwd("~/BiBC/Major_internship/tcell-family-based")
} else {setwd('~/Documents')
  .libPaths("./R")}

folder <- "./results/model3/"
run_name <- "m3_DD_try_2"
seed <- 4321
families = 500
n <- 100

bp_rule <- 'sqrt(q)'
dp_rule <- 0.5 #default
rq_rule <- 1 #only P cells will form
t_start_dist <- 'rlnorm(nr_of_families, meanlog = 1.4, sdlog = 0.3)'
t_run_rule <- 'sqrt(q) + (nr_burst_divs[i]*t_burst[i])' 
# was 0.5 in st2_nomax, but 0.6 ensures that all families can at least to their 4 burst 
# divisions, because 0.15*4 = 0.6
nr_burst_divs <- 'sample(c(2,3,4), nr_of_families, replace = TRUE)'
quality_dist <- "11.31 * rbeta(nr_of_families, shape1 = 1.45, shape2 = 2.40)"
quality_noise <- TRUE
q_noise_dist <- 'generate_dd_branch(quality[i], cells[i])'
uniform_fam <- FALSE
ASD <- FALSE
burst_time <- 0.15 #default
max_run_time <- NULL
min_t_start <- 0

recruitment_mean <- 8
recruitment_sd <- 2

source("./R/model3_run_parameter_set.R")
