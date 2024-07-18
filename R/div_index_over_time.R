
# This function returns the expected number of divisions for each branch at a 
# specific timepoint
get_exp_divs <- function(parameters, time) {
  exp_divs <- parameters %>%
    rowwise() %>%
    mutate(
      time_start_prolif = t_start + t_burst * nr_burst_divs,
      time_stop_prolif = t_start + t_run - t_correct,
      exp_divs = case_when(
        time < t_start ~ div_counter - nr_burst_divs,
        time <= time_start_prolif ~ div_counter,
        time > time_start_prolif & time < time_stop_prolif ~ (bp * (time - time_start_prolif)) + div_counter,
        TRUE ~ div_counter  # Default case, not strictly necessary here
      )
    ) %>%
    pull(exp_divs)  # Extract the exp_divs column as a vector
  
  return(exp_divs)
}

get_famsize <- function(parameters, timepoint) {
  # Helper function to compute famsize for a single row
  compute_famsize <- function(params, timepoint) {
    t_start <- as.numeric(params["t_start"])
    t_burst <- as.numeric(params["t_burst"])
    nr_burst_divs <- as.numeric(params["nr_burst_divs"])
    t_run <- as.numeric(params["t_run"])
    t_correct <- as.numeric(params["t_correct"])
    bp <- as.numeric(params["bp"])
    dp <- as.numeric(params["dp"])
    
    time_start_prolif <- t_start + t_burst * nr_burst_divs
    time_stop_prolif <- t_start + t_run - t_correct
    
    if (timepoint <= t_start) { 
      famsize <- 1 / 2^nr_burst_divs
    } else if (timepoint <= time_start_prolif) { 
      div_times <- seq(t_start + t_burst, 
                       t_start + t_burst * nr_burst_divs, 
                       by = t_burst)
      index <- sum(timepoint >= div_times)
      famsize <- 2^index / 2^nr_burst_divs
    } else if (timepoint <= time_stop_prolif) { 
      famsize <- exp(bp * (timepoint - time_start_prolif))
    } else { 
      max_cells <- exp(bp * (time_stop_prolif - time_start_prolif))
      famsize <- max_cells * exp(-dp * (timepoint - time_stop_prolif))
    }
    
    return(famsize)
  }
  
  # Check if parameters is a data.frame or a single row (vector)
  if (is.data.frame(parameters)) {
    famsizes <- apply(parameters, 1, function(row) compute_famsize(as.list(row), timepoint))
  } else {
    famsizes <- compute_famsize(as.list(parameters), timepoint)
  }
  
  return(famsizes)
}

get_div_index <- function(parameters, min_time, max_time, interval){
  times <- seq(min_time, max_time, interval)
  div_index <- sapply(times, function(t){
    exp_divs <- get_exp_divs(parameters, time = t)
    divs <- exp_divs - parameters$div_counter
    divs[divs < 0] <- 0
    nr_cells <- get_famsize(parameters, t)
    mean_div_index <- sum(exp_divs*nr_cells)/sum(nr_cells)
    mean_div_index
  })
  return(data.frame(time = times, mean_div_index = div_index))
}

prim_div_indexes <- get_div_index(prim_parameters, 0, 50, 0.1)
sec_div_indexes <- get_div_index(sec_parameters, 0, 50, 0.1)
sec_div_indexes$time <- sec_div_indexes$time + 50
ter_div_indexes <- get_div_index(ter_parameters, 0, 50, 0.1)
ter_div_indexes$time <- ter_div_indexes$time + 100

total_div_indexes <- rbind(prim_div_indexes, sec_div_indexes, ter_div_indexes)

ggplot(total_div_indexes, aes(x = time, y = mean_div_index)) +
  geom_point() +
  labs(title = "Mean division index over time",
       x = "Time (days)",
       y = "Division index") +
  theme_clean() + th
