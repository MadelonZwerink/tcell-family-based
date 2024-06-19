#!/usr/bin/env Rscript

if(Sys.info()[[4]]=="LAPTOP-3RJSLMKV") {
  print("Working from R-project")
  setwd("~/BiBC/Major_internship/tcell-family-based")
} else {setwd('~/Documents')
  .libPaths("./R")}

source("./R/speed_discrete_burst.R") 

# Define parameters
nr_sims <- 2
folder <- "./results/model1/"
run_name <- "st3_multiple_"
seed <- 4321
families <- 2

set.seed(seed)

bp_rule <- 'runif(1, min = 0, max = 3.5)'
dp_rule <- 0.5 # default
rq_rule <- 0.5 # default
t_start_dist <- 'rlnorm(nr_of_families, meanlog = 1.4, sdlog = 0.3)'
t_run_rule <- 'runif(1, min = 0.6, max = 4)' 
nr_burst_divs <- 'sample(c(2,3,4), nr_of_families, replace = TRUE)'
quality_dist <- NULL
ASD <- FALSE
burst_time <- 0.15 # default
max_run_time <- NULL
min_t_start <- 2.5

# Initialize lists
df_famsizes <- vector("list", nr_sims)
total_resp <- vector("list", nr_sims)
max_fam_table <- vector("list", nr_sims)
all_parameters_list <- vector("list", nr_sims)

#-------------------------------------------------------------------------------

# Function to generate parameters and responses
generate_parameters_multi <- function() {
  prim_parameters <- pick_parameters(bp_rule = bp_rule,
                                     dp_rule = dp_rule,
                                     rq_rule = rq_rule,
                                     t_start_dist = t_start_dist,
                                     t_run_rule = t_run_rule,
                                     nr_of_families = families,
                                     nr_burst_divs = nr_burst_divs,
                                     response_nr = 1,
                                     quality_dist = quality_dist,
                                     ASD = ASD,
                                     burst_time = burst_time,
                                     max_run_time = max_run_time,
                                     min_t_start = min_t_start)
  
  sec_parameters <- pick_parameters(response_nr = 2,
                                    prev_parameters = prim_parameters,
                                    bp_rule = bp_rule,
                                    dp_rule = dp_rule,
                                    rq_rule = rq_rule,
                                    t_start_dist = t_start_dist,
                                    t_run_rule = t_run_rule,
                                    nr_burst_divs = nr_burst_divs,
                                    quality_dist = quality_dist,
                                    ASD = ASD,
                                    burst_time = burst_time,
                                    max_run_time = max_run_time,
                                    min_t_start = min_t_start)
  
  ter_parameters <- pick_parameters(response_nr = 3,
                                    prev_parameters = sec_parameters,
                                    bp_rule = bp_rule,
                                    dp_rule = dp_rule,
                                    rq_rule = rq_rule,
                                    t_start_dist = t_start_dist,
                                    t_run_rule = t_run_rule,
                                    nr_burst_divs = nr_burst_divs,
                                    quality_dist = quality_dist,
                                    ASD = ASD,
                                    burst_time = burst_time,
                                    max_run_time = max_run_time,
                                    min_t_start = min_t_start)
  
  return(list(prim_parameters = prim_parameters,
              sec_parameters = sec_parameters,
              ter_parameters = ter_parameters))
}

# Distribution family sizes primary response 
generate_df_famsizes_multi <- function(prim_parameters) {
  df_famsizes <- generate_famsize_table_multidays(prim_parameters)
  return(df_famsizes)
}

# Total response table
generate_total_response_table_multi <- function(all_parameters_list) {
  prim_resp <- generate_total_response_table(all_parameters_list$prim_parameters, 0, 50, 0.1)
  sec_resp <- generate_total_response_table(all_parameters_list$sec_parameters, 0, 50, 0.1)
  sec_resp$time <- sec_resp$time + 50
  ter_resp <- generate_total_response_table(all_parameters_list$ter_parameters, 0, 50, 0.1)
  ter_resp$time <- ter_resp$time + 100
  
  total_resp <- rbind(prim_resp, sec_resp, ter_resp)
  
  return(total_resp)
}

# Max family table
generate_max_fam_table_multi <- function(all_parameters_list) {
  max_fam_table <- generate_max_fam_table(all_parameters_list$prim_parameters, 
                                          all_parameters_list$sec_parameters, 
                                          all_parameters_list$ter_parameters, 
                                          type = "time", 
                                          timepoint = 7)
  return(max_fam_table)
}

# Get famsize dist stats for multiple runs
get_famsize_stats_multidays_multi <- function(df_famsizes){
  # Get the statistics for each day and each simulation
  famsizes_stats <- lapply(1:4, function(day_index) {
    sapply(seq_along(df_famsizes), function(sim_index) {
      get_famsize_stats(df_famsizes[[sim_index]][[day_index]])
    })
  })  
  
  # The sapply function puts each observation in a separate column, but I would 
  # like a long format where each row is a simulation
  famsizes_stats_long <- lapply(1:4, function(day_index) {
    rownames(famsizes_stats[[day_index]]) <- c("mean", "median")
    colnames(famsizes_stats[[day_index]]) <- seq(ncol(famsizes_stats[[day_index]]))
    t(famsizes_stats[[day_index]])
  })
  
  names(famsizes_stats_long) <- c("d5", "d6", "d7", "d8")
  
  return(famsizes_stats_long)
}

generate_freq_famsize_table_multidays_multi <- function(df_famsizes){
  famsizes_freq <<- lapply(1:4, function(day_index) {
    sapply(seq_along(df_famsizes), function(sim_index) {
      freq_table <- generate_freq_famsize_table(df_famsizes[[sim_index]][[day_index]]$logfamsize)
      if(setequal(freq_table$logfamsize, c(0, seq(max(freq_table$logfamsize)))) == F){
        print(paste0("Warning: frequencies are off, simulation nr: ", sim_index,
                     " day: ", day_index + 4))}
      if(length(freq_table$freq) <= 20){
        freq_vect <- c(freq_table$freq, rep(0, 20 - length(freq_table$freq)))
      } else {
        print(paste0("Warning: unusually large family detected in simulation ", sim_index))
      }
      freq_vect
    })
  }) 
  
  # The sapply function puts each observation in a separate column, but I would 
  # like a long format where each row is a simulation
  famsizes_freq_long <- lapply(1:4, function(day_index) {
    rownames(famsizes_freq[[day_index]]) <- seq(nrow(famsizes_freq[[day_index]]))
    colnames(famsizes_freq[[day_index]]) <- seq(ncol(famsizes_freq[[day_index]]))
    as.data.frame(t(famsizes_freq[[day_index]]))
  })
  
  names(famsizes_freq_long) <- c("d5", "d6", "d7", "d8")
  
  return(famsizes_freq_long)
}

# Combine all above functions
generate_response_multi <- function(sim) {
  print(sim)
  all_parameters_list[[sim]] <<- generate_parameters_multi()
  df_famsizes[[sim]] <<- generate_df_famsizes_multi(all_parameters_list[[sim]]$prim_parameters)
  total_resp[[sim]] <<- generate_total_response_table_multi(all_parameters_list[[sim]])
  max_fam_table[[sim]] <<- generate_max_fam_table_multi(all_parameters_list[[sim]])
}

#-------------------------------------------------------------------------------

# Run simulations
lapply(seq_len(nr_sims), generate_response_multi)
famsizes_freq_table <- generate_freq_famsize_table_multidays_multi(df_famsizes)

# Calculate statistics
print("Calculating statistics...")

famsize_stats <- get_famsize_stats_multidays_multi(df_famsizes)
max_fam_stats <- sapply(max_fam_table, get_max_fam_stats, table = F)
max_fam_stats_long <- t(max_fam_stats)
total_response_stats <- sapply(total_resp, get_total_response_stats)

#-------------------------------------------------------------------------------
print("Saving data...")

saveRDS(df_famsizes, file = paste0(folder, run_name, "df_famsizes.rds"))
saveRDS(total_resp, file = paste0(folder, run_name, "total_resp.rds"))
saveRDS(max_fam_table, file = paste0(folder, run_name, "max_fam_table.rds"))
saveRDS(famsizes_freq_table, file = paste0(folder, run_name, "famsizes_freq_table.rds"))

saveRDS(famsize_stats, file = paste0(folder, run_name, "famsizes_stats.rds"))
saveRDS(max_fam_stats_long, file = paste0(folder, run_name, "max_fam_stats_long.rds"))
saveRDS(total_response_stats, file = paste0(folder, run_name, "total_response_stats.rds"))



