#!/usr/bin/env Rscript

if(Sys.info()[[4]]=="LAPTOP-3RJSLMKV") {
  print("Working from R-project")
  setwd("~/BiBC/Major_internship/tcell-family-based")
} else {setwd('~/Documents')
  .libPaths("./R")}

folder <- "./results/model2/"
run_name <- "q20"
seed <- 4321
nr_of_families = 500

#-------------------------------------------------------------------------------

bp_rule <- 'q*12.7/3'
dp_rule <- 0.5 #default
rq_rule <- 0.5 #default
t_start_dist <- 'rlnorm(nr_of_families, meanlog = 1.4, sdlog = 0.3)'
t_run_rule <- '3 + (nr_burst_divs[i]*t_burst[i])'
nr_burst_divs <- 'sample(c(2,3,4), nr_of_families, replace = TRUE)'
quality_dist <- "rbeta(nr_of_families, shape1 = 2.33, shape2 = 3.92)"
quality_noise <- FALSE
uniform_fam <- TRUE
ASD <- FALSE
burst_time <- 0.15 #default
max_run_time <- NULL
min_t_start <- 0
resolution_resp_plot <- 1

source("./R/run_parameter_set.R")

#-------------------------------------------------------------------------------