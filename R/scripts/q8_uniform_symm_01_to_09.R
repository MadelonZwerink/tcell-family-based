#!/usr/bin/env Rscript

if(Sys.info()[[4]]=="LAPTOP-3RJSLMKV") {
  print("Working from R-project")
  setwd("~/BiBC/Major_internship/tcell-family-based")
} else {setwd('~/Documents')
  .libPaths("./R")}

folder <- "./results/model2/"
run_name <- "q8_uniform_symm_01_to_09"
seed <- 4321
families = 500

#-------------------------------------------------------------------------------

min_q <- 0.1
med_q <- 0.5
max_q <- 0.9

get_alpha <- function(quality){
  alpha <- 1 / (1 - quality)
  return(alpha)
}

get_beta <- function(quality){
  beta <- 1 / quality
  return(beta)
}
#-------------------------------------------------------------------------------

bp_rule <- 'skewed_uniform(1, 0, 3.5, q)'
dp_rule <- 0.5 #default
rq_rule <- 0.5 #default
t_start_dist <- 'rlnorm(nr_of_families, meanlog = 1.4, sdlog = 0.3)'
t_run_rule <- 'skewed_uniform(1, 0.6, 4, q)'
# was 0.5 in st2_nomax, but 0.6 ensures that all families can at least to their 4 burst 
# divisions, because 0.15*4 = 0.6
nr_burst_divs <- 'sample(c(2,3,4), nr_of_families, replace = TRUE)'
quality_dist <- "runif(nr_of_families, min = min_q, max = max_q)"
quality_noise <- FALSE
uniform_fam <- TRUE
ASD <- FALSE
burst_time <- 0.15 #default
max_run_time <- NULL
min_t_start <- 2.5
resolution_resp_plot <- 1

source("./R/run_parameter_set.R")

#-------------------------------------------------------------------------------
