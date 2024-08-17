#!/usr/bin/env Rscript

if(Sys.info()[[4]]=="LAPTOP-3RJSLMKV") {
  print("Working from R-project")
  setwd("~/BiBC/Major_internship/tcell-family-based")
} else {setwd('~/Documents')
  .libPaths("./R")}

folder <- "./results/model3/"
run_name <- "m3_DD_branch"
seed <- 4321
families = 500
n <- 100

bp_rule <- 'sqrt(q)'
dp_rule <- 0.5 #default
rq_rule <- 1 #only P cells will form
t_start_dist <- 'rlnorm(nr_of_families, meanlog = 1.4, sdlog = 0.3)'
t_run_rule <- 'sqrt(q) + (nr_burst_divs[i]*t_burst[i])' 
nr_burst_divs <- 'sample(c(2,3,4), nr_of_families, replace = TRUE)'
quality_dist <- "12.7 * rbeta(nr_of_families, shape1 = 0.227, shape2 = 1.625)"
quality_noise <- TRUE
q_noise_dist <- "12.7 * rbeta(cells[i], shape1 = 0.227, shape2 = 1.625)"
uniform_fam <- FALSE
ASD <- FALSE
burst_time <- 0.15 #default
max_run_time <- NULL
min_t_start <- 0

recruitment_mean <- 8
recruitment_sd <- 2

source("./R/model3_run_parameter_set.R")
