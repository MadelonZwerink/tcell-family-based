#!/usr/bin/env Rscript

if(Sys.info()[[4]]=="LAPTOP-3RJSLMKV") {
  print("Working from R-project")
  setwd("~/BiBC/Major_internship/tcell-family-based")
} else {setwd('~/Documents')
  .libPaths("./R")}

folder <- "./results/model2/"
run_name <- "q15_norm"
seed <- 4321
families = 500

#-------------------------------------------------------------------------------

bp_rule <- 'q*3.5'
dp_rule <- 0.5 #default
rq_rule <- 0.5 #default
t_start_dist <- 'rlnorm(nr_of_families, meanlog = 1.4, sdlog = 0.3)'
t_run_rule <- '(q*(4 - 0.6)) + 0.6'
# was 0.5 in st2_nomax, but 0.6 ensures that all families can at least to their 4 burst 
# divisions, because 0.15*4 = 0.6
nr_burst_divs <- 'sample(c(2,3,4), nr_of_families, replace = TRUE)'
quality_dist <- "rnorm(nr_of_families, mean = 0.52, sd = 0.127)"
quality_noise <- FALSE
q_noise_dist <- NULL
uniform_fam <- FALSE
ASD <- FALSE
burst_time <- 0.15 #default
max_run_time <- NULL
min_t_start <- 2.5
resolution_resp_plot <- 1

source("./R/run_parameter_set.R")

#-------------------------------------------------------------------------------
