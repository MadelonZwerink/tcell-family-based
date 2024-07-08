SSR_famsize_dist <- function(prim_parameters) {
  exsmall_size <- 1
  pday5 <- plot_famsize_dist(prim_parameters, timepoint = 5, method = "time", binwidth = 1)
  pday6 <- plot_famsize_dist(prim_parameters, timepoint = 6, method = "time", binwidth = 1)
  pday7 <- plot_famsize_dist(prim_parameters, timepoint = 7, method = "time", binwidth = 1)
  pday8 <- plot_famsize_dist(prim_parameters, timepoint = 8, method = "time", binwidth = 1)
  
  model_famsizes_day5 <- as.data.frame(table(round(pday5[[3]]$logfamsizes)))
  names(model_famsizes_day5) <- c("Amount", "Freq")
  famsizes_day5 <- left_join(pday5[[4]], model_famsizes_day5, join_by(Amount), suffix = c("_Exp", "_Model"))
  famsizes_day5$Freq_Model[is.na(famsizes_day5$Freq_Model)] <- 0
  famsizes_day5_exsmall <- famsizes_day5[as.numeric(famsizes_day5$Amount) > exsmall_size, ]
  
  famsizes_day5$SR <- (famsizes_day5$Freq_Exp - famsizes_day5$Freq_Model)^2
  SSR_day5 <- sum(famsizes_day5$SR)
  
  famsizes_day5_exsmall$SR <- (famsizes_day5_exsmall$Freq_Exp - famsizes_day5_exsmall$Freq_Model)^2
  SSR_day5_exsmall <- sum(famsizes_day5_exsmall$SR)
  
  model_famsizes_day6 <- as.data.frame(table(round(pday6[[3]]$logfamsizes)))
  names(model_famsizes_day6) <- c("Amount", "Freq")
  famsizes_day6 <- left_join(pday6[[5]], model_famsizes_day6, join_by(Amount), suffix = c("_Exp", "_Model"))
  famsizes_day6$Freq_Model[is.na(famsizes_day6$Freq_Model)] <- 0
  
  famsizes_day6_exsmall <- famsizes_day6[as.numeric(famsizes_day6$Amount) > exsmall_size, ]
  
  famsizes_day6$SR <- (famsizes_day6$Freq_Exp - famsizes_day6$Freq_Model)^2
  SSR_day6 <- sum(famsizes_day6$SR)
  
  famsizes_day6_exsmall$SR <- (famsizes_day6_exsmall$Freq_Exp - famsizes_day6_exsmall$Freq_Model)^2
  SSR_day6_exsmall <- sum(famsizes_day6_exsmall$SR)
  
  model_famsizes_day7 <- as.data.frame(table(round(pday7[[3]]$logfamsizes)))
  names(model_famsizes_day7) <- c("Amount", "Freq")
  famsizes_day7 <- left_join(pday7[[6]], model_famsizes_day7, join_by(Amount), suffix = c("_Exp", "_Model"))
  famsizes_day7$Freq_Model[is.na(famsizes_day7$Freq_Model)] <- 0
  famsizes_day7_exsmall <- famsizes_day7[as.numeric(famsizes_day7$Amount) > exsmall_size, ]
  
  famsizes_day7$SR <- (famsizes_day7$Freq_Exp - famsizes_day7$Freq_Model)^2
  SSR_day7 <- sum(famsizes_day7$SR)
  
  famsizes_day7_exsmall$SR <- (famsizes_day7_exsmall$Freq_Exp - famsizes_day7_exsmall$Freq_Model)^2
  SSR_day7_exsmall <- sum(famsizes_day7_exsmall$SR)
  
  model_famsizes_day8 <- as.data.frame(table(round(pday8[[3]]$logfamsizes)))
  names(model_famsizes_day8) <- c("Amount", "Freq")
  famsizes_day8 <- left_join(pday8[[7]], model_famsizes_day8, join_by(Amount), suffix = c("_Exp", "_Model"))
  famsizes_day8$Freq_Model[is.na(famsizes_day8$Freq_Model)] <- 0
  
  famsizes_day8_exsmall <- famsizes_day8[as.numeric(famsizes_day8$Amount) > exsmall_size, ]
  
  famsizes_day8$SR <- (famsizes_day8$Freq_Exp - famsizes_day8$Freq_Model)^2
  SSR_day8 <- sum(famsizes_day8$SR)
  
  famsizes_day8_exsmall$SR <- (famsizes_day8_exsmall$Freq_Exp - famsizes_day8_exsmall$Freq_Model)^2
  SSR_day8_exsmall <- sum(famsizes_day8_exsmall$SR)
  
  return(list(
    list(SSR_day5, SSR_day6, SSR_day7, SSR_day8),
    list(SSR_day5_exsmall, SSR_day6_exsmall, SSR_day7_exsmall, SSR_day8_exsmall),
    list(famsizes_day5, famsizes_day6, famsizes_day7, famsizes_day8)
  ))
}

#-------------------------------------------------------------------------------

run_model <- function(parameters,
                      t_total = 100,
                      response_nr = 1,
                      t_previous = 0,
                      prev_results = NULL) {
  n_sol <- list()
  parameters <<- parameters
  
  for (i in 1:nrow(parameters)) {
    t_run <- parameters$t_run[i]
    t_correct <- parameters$t_correct[i]
    
    input_vectors <- generate_input_vectors(n = n,
                                            parameters = parameters,
                                            i = i)
    # make sure that cells do not proliferate after tRun
    this_run <- run(tmax = t_total - parameters$t_start[i] - 0.01,
                    tstep = 0.01,
                    state = input_vectors,
                    parms = c(dp = parameters$dp[i],
                              bp = parameters$bp[i],
                              dq = parameters$dq[i],
                              bq = parameters$bq[i]),
                    table = TRUE,
                    arrest = (t_run - t_correct),
                    after = "if(t==t_run - t_correct)parms[\"bp\"]<-0",
                    timeplot = FALSE)
    
    # add time before the start of replication to the time
    this_run$time <- this_run$time + parameters$t_start[i]
    
    # add the values before start of proliferation and set to input vectors
    for (time_point in seq(0, parameters$t_start[i] - 0.005, 0.01)) {
      this_run[nrow(this_run) + 1, ] <- c(time_point, input_vectors)
    }
    
    this_run$time <- round(this_run$time, digits = 2) + t_previous + 0.01
    this_run <- this_run[order(this_run$time), ]
    this_run %<>% mutate(fam_nr = parameters$fam_nr[i],
                         fam_nr_2 = 0,
                         fam_nr_3 = 0)
    
    if (response_nr == 2) {
      this_run %<>% mutate(fam_nr_2 = i)
    }
    if (response_nr == 3) {
      this_run %<>% mutate(fam_nr_2 = prev_results$fam_nr_2[i],
                           fam_nr_3 = i)
    }
    
    n_sol[[i]] <- this_run
    print(i)
  }
  return(n_sol)
}

#-------------------------------------------------------------------------------

generate_input_vectors <- function(n,
                                   parameters,
                                   i) {
  P <- rep(0, n)
  Q <- rep(0, n)
  if (parameters$cell_type[i] == "P") {
    P[parameters$div_counter[i]] <- 1
  } else {
    Q[parameters$div_counter[i]] <- 1
  }
  names(P) <- paste0("P", seq(1, n))
  names(Q) <- paste0("Q", seq(1, n))
  return(c(P, Q))
}

#-------------------------------------------------------------------------------

plot_h_grid_famsize_dist <- function(famsizes_table, x_axis_max = NULL, y_axis_max = NULL){
  frequencies_table <- generate_freq_famsize_table_multidays(famsizes_table)
  
  if(is.null(x_axis_max)) {
    x_axis_max <- rep(max(frequencies_table[[4]]$logfamsize), 4) + 1
  }
  if(is.null(y_axis_max)){
    y_axis_max <- rep(max(frequencies_table[[1]]$freq) + 0.05, 4)
  }
  
  plots <- lapply(5:8, function(day) {
    plot_famsize_distribution(frequencies_table[[day-4]], day, 
                              x_axis_max = x_axis_max[day-4], y_axis_max = y_axis_max[day-4])
  })
  plots_no_ylabel <- lapply(2:4, function(day) {
    plots[[day]] + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
  })
  plots_one_ylabel <- list(plots[[1]], plots_no_ylabel[[1]], plots_no_ylabel[[2]], plots_no_ylabel[[3]])
  
  gt <- arrangeGrob(grobs = plots_one_ylabel, ncol = 4, widths = c(4.6,4,4,4))
  plot <- as_ggplot(gt) +                                # transform to a ggplot
    draw_plot_label(label = c("A - day 5", "B - day 6", "C - day 7", "D - day 8"), size = 13,
                    x = c(0.015, 0.25, 0.49, 0.73), y = c(1, 1, 1, 1))
  
  return(plot)
}

#-------------------------------------------------------------------------------

plot_sq_grid_famsize_dist <- function(famsizes_table, x_axis_max = NULL, y_axis_max = NULL){
  frequencies_table <- lapply(1:4, function(day) {
    generate_freq_famsize_table(famsizes_table[[day]]$logfamsize)
  })
  if(is.null(x_axis_max)){
    x_axis_max <- c(max(frequencies_table[[1]]$logfamsize, frequencies_table[[3]]$logfamsize) + 1,
                    max(frequencies_table[[2]]$logfamsize, frequencies_table[[4]]$logfamsize) + 1,
                    max(frequencies_table[[1]]$logfamsize, frequencies_table[[3]]$logfamsize) + 1,
                    max(frequencies_table[[2]]$logfamsize, frequencies_table[[4]]$logfamsize) + 1)
  }
  if(is.null(y_axis_max)){
    y_axis_max <- c(max(frequencies_table[[1]]$freq, frequencies_table[[2]]$freq) + 0.05, 
                    max(frequencies_table[[1]]$freq, frequencies_table[[2]]$freq) + 0.05,
                    max(frequencies_table[[3]]$freq, frequencies_table[[4]]$freq) + 0.05,
                    max(frequencies_table[[3]]$freq, frequencies_table[[4]]$freq) + 0.05)
  }
  
  plots <- lapply(5:8, function(day) {
    plot_famsize_distribution(frequencies_table[[day-4]], day, 
                              x_axis_max = x_axis_max[day-4], y_axis_max = y_axis_max[day-4])
  })
  
  plots_two_ylabel <- list(plots[[1]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
                                              plot.margin = margin(t = 0.5, r = 0, b = 0, l = 0, unit = "cm")), 
                           plots[[2]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
                                              axis.text.y = element_blank(), axis.title.y = element_blank(),
                                              plot.margin = margin(t = 0.5, r = 0, b = 0, l = 0, unit = "cm")),
                           plots[[3]] + theme(plot.margin = margin(t = 0.2, r = 0, b = 0.5, l = 0, unit = "cm")), 
                           plots[[4]] + theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
                                              plot.margin = margin(t = 0.2, r = 0, b = 0.5, l = 0, unit = "cm")))
  
  gt <- arrangeGrob(grobs = plots_two_ylabel, ncol = 2, widths = c(4.8, 4), heights = c(2, 2.5))
  plot <- as_ggplot(gt) +                                # transform to a ggplot
    draw_plot_label(label = c("A - day 5", "B - day 6", "C - day 7", "D - day 8"), size = 12,
                    x = c(0.04, 0.515, 0.04, 0.515), y = c(0.98, 0.98, 0.55, 0.55))
  
  return(plot)
}

#-------------------------------------------------------------------------------

get_famsize_old <- function(parameters, timepoint){
  time_start_prolif <- parameters$t_start + parameters$t_burst*parameters$nr_burst_divs
  time_stop_prolif <- parameters$t_start + parameters$t_run - parameters$t_correct
  
  if (timepoint <= parameters$t_start) { # before family starts dividing
    famsize <- 1 / 2^parameters$nr_burst_divs 
    # divide by the number of cells formed at the end of the burst divisions
    # because otherwise each branch will contain all burst cells
  } else if (timepoint <= time_start_prolif) { # while family is in burst division
    burst_divs_v <- seq(parameters$nr_burst_divs)
    div_times <- (burst_divs_v * parameters$t_burst) + parameters$t_start # creates a vector with the times where a cell divides
    differences <- div_times - timepoint
    index <- which(min(differences[differences >= 0]) == differences)
    famsize <- (2^burst_divs_v[index]) / 2^parameters$nr_burst_divs 
    # divide by the number of cells formed at the end of the burst divisions
    # because otherwise each branch will contain all burst cells
  } else if (timepoint <= time_stop_prolif && 
             timepoint > time_start_prolif) { # while family is dividing
    famsize <- exp((parameters$bp) * (timepoint - time_start_prolif))
  } else { # after family stopped dividing
    max_cells <- exp((parameters$bp) * (time_stop_prolif - time_start_prolif))
    famsize <- max_cells * exp(-1 * parameters$dp * (timepoint - time_stop_prolif))
  }
  
  return(famsize)
}

#-------------------------------------------------------------------------------

solve_function_old <- function(parameters, min_time, max_time, interval) {
  num_timepoints <- length(seq(min_time, max_time, interval))
  P_total <- data.frame(matrix(0, nrow = nrow(parameters), ncol = num_timepoints + 3))
  colnames(P_total) <- c("fam_nr", "fam_nr_2", "fam_nr_3", seq(min_time, max_time, by = interval))
  P_total$fam_nr <- as.factor(parameters$fam_nr)
  P_total$fam_nr_2 <- as.factor(parameters$fam_nr_2)
  P_total$fam_nr_3 <- as.factor(parameters$fam_nr_3)
  
  for (cell in 1:nrow(parameters)) {
    for (timepoint in seq(min_time, max_time, by = interval)) {
      col_index <- which(colnames(P_total) == timepoint)
      P_total[cell, col_index] <- get_famsize(parameters[cell,], timepoint)
    }
    
  }
  return(P_total)
}

#-------------------------------------------------------------------------------
