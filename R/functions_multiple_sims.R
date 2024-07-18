#!/usr/bin/env Rscript

if(Sys.info()[[4]]=="LAPTOP-3RJSLMKV") {
  print("Working from R-project")
  setwd("~/BiBC/Major_internship/tcell-family-based")
} else {setwd('~/Documents')
  .libPaths("./R")}

source("./R/speed_discrete_burst.R")


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
        freq_vect <- freq_table$freq 
      }
      # Ensure frequency vector length is 20, pad with zeros if necessary
      if (length(freq_table$freq) <= 25) {
        freq_vect <- c(freq_table$freq, rep(0, 25 - length(freq_table$freq)))
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
    rownames(day_data) <- seq(0, nrow(day_data) - 1)
    as.data.frame(t(day_data))
  })
  
  # Name the list elements for each day
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

# Otherwise the function returns the statistics for the whole response, we want
# statistics for each response
get_response_stats_part <- function(total_resp, max_time){
  get_total_response_stats(total_resp[total_resp$time < max_time, ])
}

convert_cols_to_numeric <- function(df){
  if(is.data.frame(df) == F){
    df <- as.data.frame(df)
  }
  i <- seq(ncol(df))
  df[, i] <- apply(df[, i], 2, function(x) as.numeric(as.character(x)))
  return(df)
}

convert_response_stats_part <- function(total_resp, max_time){
  part_response_stats <- as.data.frame(t(sapply(total_resp, get_response_stats_part, c(max_time = max_time))))
  part_response_stats <- convert_cols_to_numeric(part_response_stats)
  return(part_response_stats)
}

#-------------------------------------------------------------------------------
# These functions are used in multiple_runs_plots_stats
#-------------------------------------------------------------------------------

# Function to get mean and sd for each day
generate_mean_freq_table <- function(famsizes_freq_table){
  sd_famsizes_freq <- famsizes_freq_table %>%
    summarize(across(everything(), sd, na.rm = TRUE)) %>%
    t()
  
  famsizes_freq_table_mean <- data.frame(logfamsize = as.numeric(colnames(famsizes_freq_table)),
                                         freq_mean = colMeans(famsizes_freq_table),
                                         freq_sd = sd_famsizes_freq,
                                         freq_confint = 1.96*sd_famsizes_freq)
  
  return(famsizes_freq_table_mean)
}

get_processed_famsize_stats <- function(famsize_stats){
  mean <- mean(famsize_stats$mean)  
  sd_mean <- sd(famsize_stats$mean)
  median <- mean(famsize_stats$median)  
  sd_median <- sd(famsize_stats$median)

  return(c(mean, sd_mean, median, sd_median))
}

#This function returns mean, sd and lower and upper bound of 95% CI for 
#dataframes in which each column is a variable and each row an observation
get_processed_stats <- function(max_fam_stats_long, scientific = F){
  means <- colMeans(max_fam_stats_long)
  sd <- apply(max_fam_stats_long, 2, sd)
  margin_of_error <- 1.96 * (sd / sqrt(nrow(max_fam_stats_long)))
  lower_bound <- means - margin_of_error
  upper_bound <- means + margin_of_error
  
  if(scientific == F){
    for(output_stats in c("means", "sd", "lower_bound", "upper_bound")){
      df <- get(output_stats)  # Retrieve the data frame by its name
      df <- as.numeric(format(round(df, digits = 4), scientific = FALSE))  # Round and format
      assign(output_stats, df)  # Assign the modified data frame back to the original name
    }
  }
  return(data.frame("variable" = colnames(max_fam_stats_long),
                    "means" = c(means), 
                    "sd" = c(sd), 
                    "lower_bound" = c(lower_bound), 
                    "upper_bound" = c(upper_bound)))
}

get_Q_famsize_stats_multi <- function(max_fam_table){
  Q_famsize_stats <- get_Q_famsize_stats(max_fam_table)
  Q_famsize_stats_row <- c(Q_famsize_stats$r, Q_famsize_stats$p_value)
  return(Q_famsize_stats_row)
}
