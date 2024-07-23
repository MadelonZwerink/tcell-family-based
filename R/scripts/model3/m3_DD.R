#!/usr/bin/env Rscript

if(Sys.info()[[4]]=="LAPTOP-3RJSLMKV") {
  print("Working from R-project")
  setwd("~/BiBC/Major_internship/tcell-family-based")
} else {setwd('~/Documents')
  .libPaths("./R")}

source("./R/m3_functions_multiple_sims.R") 

# Define parameters
nr_sims <- 100
folder <- "./results/model3/"
run_name <- "m3_DD"
seed <- 4321
nr_of_families <- 500
families <- 500
set.seed(seed)

bp_rule <- 'runif(1, min = 0, max = 3.5)'
dp_rule <- 0.5 #default
rq_rule <- 1 
t_start_dist <- 'rlnorm(nr_of_families, meanlog = 1.4, sdlog = 0.3)'
t_run_rule <- 'runif(1, min = 0.6, max = 4)' 
nr_burst_divs <- 'sample(c(2,3,4), nr_of_families, replace = TRUE)'
quality_dist <- "rbeta(nr_of_families, alpha = 1.42, beta = 2.66)"
quality_noise <- FALSE
ASD <- FALSE
burst_time <- 0.15 #default
uniform_fam <- FALSE
max_run_time <- NULL
min_t_start <- 0
resolution_resp_plot <- 0.1

recruitment_mean <- 6
recruitment_sd <- 2

# Initialize lists
df_famsizes <- vector("list", nr_sims)
total_resp <- vector("list", nr_sims)
max_fam_table <- vector("list", nr_sims)
all_parameters_list <- vector("list", nr_sims)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Run simulations
lapply(seq_len(nr_sims), generate_response_multi)
famsizes_freq_table <- generate_freq_famsize_table_multidays_multi(df_famsizes)

# Calculate statistics
print("Calculating statistics...")

famsize_stats <- get_famsize_stats_multidays_multi(df_famsizes)
max_fam_stats <- sapply(max_fam_table, get_max_fam_stats, table = F)
max_fam_stats_long <- t(max_fam_stats)

prim_response_stats <- convert_response_stats_part(total_resp, max_time = 50, min_time = 0)
sec_response_stats <- convert_response_stats_part(total_resp, max_time = 100, min_time = 50)
ter_response_stats <- convert_response_stats_part(total_resp, max_time = 150, min_time = 100)

#-------------------------------------------------------------------------------
print("Saving data...")

saveRDS(df_famsizes, file = paste0(folder, run_name, "df_famsizes.rds"))
saveRDS(total_resp, file = paste0(folder, run_name, "total_resp.rds"))
saveRDS(max_fam_table, file = paste0(folder, run_name, "max_fam_table.rds"))
saveRDS(famsizes_freq_table, file = paste0(folder, run_name, "famsizes_freq_table.rds"))

saveRDS(famsize_stats, file = paste0(folder, run_name, "famsizes_stats.rds"))
saveRDS(max_fam_stats_long, file = paste0(folder, run_name, "max_fam_stats_long.rds"))
saveRDS(prim_response_stats, file = paste0(folder, run_name, "prim_response_stats.rds"))
saveRDS(sec_response_stats, file = paste0(folder, run_name, "sec_response_stats.rds"))
saveRDS(ter_response_stats, file = paste0(folder, run_name, "ter_response_stats.rds"))



