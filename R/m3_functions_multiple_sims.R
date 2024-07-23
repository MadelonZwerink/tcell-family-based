#!/usr/bin/env Rscript

if(Sys.info()[[4]]=="LAPTOP-3RJSLMKV") {
  print("Working from R-project")
  setwd("~/BiBC/Major_internship/tcell-family-based")
} else {setwd('~/Documents')
  .libPaths("./R")}

source("./R/speed_discrete_burst.R")
source("./R/model3_functions.R")


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# FUNCTIONS USED FOR MULTIPLE SIMS
#-------------------------------------------------------------------------------
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
  
  div_index_prim <- generate_div_index_df(prim_parameters)
  
  recruited_cells_prim_parameters <- generate_recruited_cells(div_index_prim, 
                                                              recruitment_mean,
                                                              recruitment_sd)
  recruited_cells_prim <- generate_processed_recruited_cells(recruited_cells_prim_parameters)
  
  #-------------------------------------------------------------------------------
  
  sec_parameters <- pick_parameters(response_nr = 2,
                                    prev_parameters = recruited_cells_prim,
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
  
  div_index_sec <- generate_div_index_df(sec_parameters)
  
  recruited_cells_sec_parameters <- generate_recruited_cells(div_index_sec, 
                                                             recruitment_mean,
                                                             recruitment_sd)
  recruited_cells_sec <- generate_processed_recruited_cells(recruited_cells_sec_parameters)
  
  
  #-------------------------------------------------------------------------------
  
  ter_parameters <- pick_parameters(response_nr = 3,
                                    prev_parameters = recruited_cells_sec,
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
  
  div_index_ter <- generate_div_index_df(ter_parameters)
  
  recruited_cells_ter_parameters <- generate_recruited_cells(div_index_ter, 
                                                             recruitment_mean,
                                                             recruitment_sd)
  recruited_cells_ter <- generate_processed_recruited_cells(recruited_cells_ter_parameters)
  
  
  return(list(prim_parameters = prim_parameters,
              sec_parameters = sec_parameters,
              ter_parameters = ter_parameters,
              recruited_cells_prim_parameters = recruited_cells_prim_parameters,
              recruited_cells_sec_parameters = recruited_cells_sec_parameters,
              recruited_cells_ter_parameters = recruited_cells_ter_parameters,
              recruited_cells_prim = recruited_cells_prim,
              recruited_cells_sec = recruited_cells_sec,
              recruited_cells_ter = recruited_cells_ter
              ))
}

# Distribution family sizes primary response 
generate_df_famsizes_multi <- function(prim_parameters) {
  df_famsizes <- generate_famsize_table_multidays(prim_parameters)
  return(df_famsizes)
}

# Total response table
generate_total_response_table_multi <- function(all_parameters_list) {
  prim_resp <- generate_total_response_table(all_parameters_list$prim_parameters, 0, 50, 0.1,
                                             model = "noQ", 
                                             recruited_cells = all_parameters_list$recruited_cells_prim_parameters)
  sec_resp <- generate_total_response_table(all_parameters_list$sec_parameters, 0, 50, 0.1,
                                            model = "noQ", 
                                            recruited_cells = all_parameters_list$recruited_cells_sec_parameters)
  sec_resp$time <- sec_resp$time + 50
  ter_resp <- generate_total_response_table(all_parameters_list$ter_parameters, 0, 50, 0.1,
                                            model = "noQ", 
                                            recruited_cells = all_parameters_list$recruited_cells_ter_parameters)
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
                                          timepoint = 7,
                                          model = "noQ",
                                          recruited_cells_prim = all_parameters_list$recruited_cells_prim,
                                          recruited_cells_sec = all_parameters_list$recruited_cells_sec,
                                          recruited_cells_ter = all_parameters_list$recruited_cells_ter)
  return(max_fam_table)
}

# Combine all above functions
generate_response_multi <- function(sim) {
  print(sim)
  all_parameters_list[[sim]] <<- generate_parameters_multi()
  df_famsizes[[sim]] <<- generate_df_famsizes_multi(all_parameters_list[[sim]]$prim_parameters)
  total_resp[[sim]] <<- generate_total_response_table_multi(all_parameters_list[[sim]])
  max_fam_table[[sim]] <<- generate_max_fam_table_multi(all_parameters_list[[sim]])
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
    rownames(famsizes_stats[[day_index]]) <- c("mean", "median","sd", "disparity")
    colnames(famsizes_stats[[day_index]]) <- seq(ncol(famsizes_stats[[day_index]]))
    as.data.frame(t(famsizes_stats[[day_index]]))
  })
  
  names(famsizes_stats_long) <- c("d5", "d6", "d7", "d8")
  
  return(famsizes_stats_long)
}

generate_freq_famsize_table_multidays_multi <- function(df_famsizes) {
  # Initialize a list to store frequency tables for each day
  famsizes_freq <<- lapply(1:4, function(day_index) {
    # Process each simulation for the current day
    sapply(seq_along(df_famsizes), function(sim_index) {
      # Generate frequency table for logfamsize on the current day
      freq_table <- generate_freq_famsize_table(df_famsizes[[sim_index]][[day_index]]$logfamsize)
      
      # Check for missing bins in the frequency table
      expected_bins <- c(0, seq(max(freq_table$logfamsize)))
      if (!setequal(freq_table$logfamsize, expected_bins)) {
        warning(paste0("Warning: frequencies are off, simulation nr: ", sim_index,
                       " day: ", day_index + 4))
      }
      
      # Ensure frequency vector length is 20, pad with zeros if necessary
      if (length(freq_table$freq) <= 20) {
        freq_vect <- c(freq_table$freq, rep(0, 20 - length(freq_table$freq)))
      } else {
        warning(paste0("Warning: unusually large family detected in simulation ", sim_index))
        freq_vect <- freq_table$freq 
      }
      
      return(freq_vect)
    })
  })
  
  # Transform the frequency tables to long format where each row is a simulation
  famsizes_freq_long <- lapply(1:4, function(day_index) {
    day_data <- famsizes_freq[[day_index]]
    rownames(day_data) <- seq(nrow(day_data))
    colnames(day_data) <- seq(0, ncol(day_data) - 1)
    as.data.frame(t(day_data))
  })
  
  # Name the list elements for each day
  names(famsizes_freq_long) <- c("d5", "d6", "d7", "d8")
  
  return(famsizes_freq_long)
}

# Otherwise the function returns the statistics for the whole response, we want
# statistics for each response
get_response_stats_part <- function(total_resp, max_time, min_time){
  get_total_response_stats(total_resp[total_resp$time < max_time & total_resp$time > min_time, ])
}

convert_cols_to_numeric <- function(df){
  if(is.data.frame(df) == F){
    df <- as.data.frame(df)
  }
  i <- seq(ncol(df))
  df[, i] <- apply(df[, i], 2, function(x) as.numeric(as.character(x)))
  return(df)
}

convert_response_stats_part <- function(total_resp, max_time, min_time){
  part_response_stats <- as.data.frame(t(sapply(1:length(total_resp), function(i) {
    get_response_stats_part(total_resp[[i]], max_time, min_time)
    })))
  part_response_stats <- convert_cols_to_numeric(part_response_stats)
  return(part_response_stats)
}
