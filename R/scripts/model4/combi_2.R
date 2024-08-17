#!/usr/bin/env Rscript

if(Sys.info()[[4]]=="LAPTOP-3RJSLMKV") {
  print("Working from R-project")
  setwd("~/BiBC/Major_internship/tcell-family-based")
} else {setwd('~/Documents')
  .libPaths("./R")}

mean(prim_parameters$quality)
median(prim_parameters$quality)
hist(prim_parameters$quality)
local <- runif(4896, 0.5, 1.5)
mean(prim_parameters$quality * local)
median(prim_parameters$quality * local)
hist(prim_parameters$quality * local)

folder <- "./results/model4/"
run_name <- "combi_2"
seed <- 4321
nr_of_families = 500

#-------------------------------------------------------------------------------

bp_rule <- 'sqrt(q*11.13)'
dp_rule <- 0.5 #default
rq_rule <- 0.5 #default
t_start_dist <- 'rlnorm(nr_of_families, meanlog = 1.4, sdlog = 0.3)'
t_run_rule <- 'sqrt(q*11.13) + (nr_burst_divs[i]*t_burst[i])'
nr_burst_divs <- 'sample(c(2,3,4), nr_of_families, replace = TRUE)'
quality_dist <- "rbeta(nr_of_families, shape1 = 2.165, shape2 = 3.295)"
quality_dist_2 <- 'generate_fixed_q_vector(prev_parameters)'
quality_noise <- TRUE
q_noise_dist <- "round(rnorm(cells[i], mean = quality[i], sd = 0.1), digits = 2)"
uniform_fam <- FALSE
ASD <- FALSE
burst_time <- 0.15 #default
max_run_time <- NULL
min_t_start <- 0
resolution_resp_plot <- 0.1

source("./R/run_parameter_set_fixed_q.R")

#-------------------------------------------------------------------------------

