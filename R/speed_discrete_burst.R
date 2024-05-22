speed_version <- "31-01-2024"

library(tidyverse)  # CRAN v2.0.0
library(magrittr)   # CRAN v2.0.3
library(ggplot2)    # CRAN v3.4.4
library(purrr)      # CRAN v1.0.2
library(data.table) # CRAN v1.15.0
library(cowplot)    # CRAN v1.1.3
library(tibble)     # CRAN v3.2.1
library(ggthemes)
library(gridExtra)

source("./R/data_load.R")

#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------

pick_parameters <- function(
    bp_rule = "runif(1, min = 0.5, max = 3)",
    dp_rule = 0.4,
    rq_rule = 0.5,
    t_start_dist,
    t_run_rule,
    nr_of_families = NULL,
    nr_burst_divs = 3,
    correction_t_burst = NULL,
    response_nr = 1,
    prev_parameters = NULL,
    quality_dist = NULL,
    ASD = FALSE,
    burst_time = 0.15,
    max_run_time = 8,
    min_t_start = NULL) {
  
  parameters <- data.frame(cell_type = factor(),
                           div_counter = integer(),
                           t_start = numeric(),
                           quality = numeric(),
                           fam_nr = integer(),
                           fam_nr_2 = integer(),
                           fam_nr_3 = integer(),
                           nr_burst_divs = integer(),
                           t_burst = numeric(),
                           t_run = numeric(),
                           t_correct = numeric(),
                           bp = numeric(),
                           dp = numeric())
  
  if (is.null(prev_parameters) == FALSE) {
    prev_q <- prev_parameters[prev_parameters$cell_type == "Q", ]
    nr_of_families <- nrow(prev_q)
    prev_divs <- prev_q$div_counter
  } else {
    prev_divs <- rep(0, nr_of_families)
  }
  
  t_start_expr <- substitute(eval(parse(text = t_start_dist)))
  t_start_val <- eval(t_start_expr)
  t_start <- round(t_start_val, digits = 2)
  
  if (!is.null(quality_dist)) {
    quality <- round(eval(parse(text = paste(quality_dist))), digits = 2)
    # beta distribution: rbeta(nr_of_families, shape1 = 2, shape2 = 5)
    # uniform distribution: runif(nr_of_families, 0, 1)
  } else {
    quality <- rep(NA, nr_of_families)
  }
  
  nr_burst_divs <- eval(parse(text = nr_burst_divs))
  cells <- 2^nr_burst_divs
  div_counter <- nr_burst_divs + prev_divs
  t_burst <- rep(burst_time, nr_of_families)

  
  for (i in seq(nr_of_families)) {
    # Make sure the families have the correct family number
    if (response_nr == 1) {
      fam_nr <- i
      fam_nr_2 <- 0
      fam_nr_3 <- 0
    }
    if (response_nr == 2) {
      fam_nr <- prev_q$fam_nr[i]
      fam_nr_2 <- i
      fam_nr_3 <- 0
    }
    if (response_nr == 3) {
      fam_nr <- prev_q$fam_nr[i]
      fam_nr_2 <- prev_q$fam_nr_2[i]
      fam_nr_3 <- i
    }
    if(ASD == F) {
      cell_type <- rbinom(cells[i], 1, eval(parse(text = rq_rule)))
     # if(sum(cell_type == 0)){cell_type[1] = 1}
      cell_type <- replace(cell_type, which(cell_type == 1), "P")
      cell_type <- replace(cell_type, which(cell_type == 0), "Q")
    } else if (ASD == T){
      cell_type <- c("Q", rep("P", cells[i] - 1))
    }
    
    for(cell in seq(cells[i])){
      # First define the expression and then evaluate the expression
      # If this is not included, the expression is only evaluated once,
      # leading to the same t_run for all (sub)families
      t_run_expr <- substitute(eval(parse(text = t_run_rule)))
      t_run_val <- eval(t_run_expr)
      t_run <- round(t_run_val, digits = 2)
      
      # this is to make sure all proliferation stops at max_run_time
      if (t_start[i] + t_run >= max_run_time) {
        t_correct <- round((t_start[i] + t_run - max_run_time), digits = 2)
      } else {
        t_correct <- 0
      }
      if (t_start[i] < 0) {
        t_start[i] <- 0.01
      }
      if (!is.null(min_t_start)){
        if (t_start[i] < min_t_start) {
          t_start[i] <- runif(1, min = min_t_start, max = (min_t_start + 1.5))
      }}
      
      if (cell_type[cell] == "P"){
        bp <- round(eval(parse(text = bp_rule)), digits = 2)
        dp <- round(eval(parse(text = dp_rule)), digits = 2)
      } else if (cell_type[cell] == "Q"){
        bp <- 0
        dp <- 0
      }
      
      parameters %<>% add_row(
        cell_type = cell_type[cell],
        div_counter = div_counter[i],
        t_start = t_start[i],
        quality = quality[i],
        fam_nr = fam_nr,
        fam_nr_2 = fam_nr_2,
        fam_nr_3 = fam_nr_3,
        nr_burst_divs = nr_burst_divs[i],
        t_burst = t_burst[i],
        t_run = t_run,
        t_correct = t_correct,
        bp = bp,
        dp = dp
      )
    }
  }
  return(parameters)
}


#-------------------------------------------------------------------------------

get_famsize <- function(parameters, timepoint){
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

get_max_famsize <- function(parameters){
  time_start_prolif <- parameters$t_start + parameters$t_burst*parameters$nr_burst_divs
  time_stop_prolif <- parameters$t_start + parameters$t_run - parameters$t_correct
  max_cells <- exp((parameters$bp) * (time_stop_prolif - time_start_prolif))
  return(max_cells)
}

#-------------------------------------------------------------------------------

solve_function <- function(parameters, min_time, max_time, interval) {
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

generate_total_response_table <- function(parameters, min_time, max_time, interval){
  famsizes_table <- solve_function(parameters, min_time, max_time, interval)
  timepoints <- seq(min_time, max_time, interval)
  num_timepoints <- length(timepoints)
  total_resp <- apply(famsizes_table[,4 : (num_timepoints + 3)], 2, FUN = sum)
  response <- data.frame(time = timepoints, cells = total_resp, log_cells = log10(total_resp))
  
  return(response)
}

#-------------------------------------------------------------------------------

plot_response <- function(response){
  plot <- ggplot(data = response, aes(x = time, y = log_cells)) +
    geom_point() +
    labs(
      x = "Time (days)",
      y = "Number of cells (log)"
    ) + theme_clean()
  return(plot)
}

#-------------------------------------------------------------------------------

generate_max_fam_table <- function(prim_parameters, 
                                  sec_parameters, 
                                  ter_parameters, 
                                  type = c("time", "max"), 
                                  timepoint = NULL, 
                                  remove_zeroes = FALSE) {
  
  nr_of_families <- max(prim_parameters$fam_nr)
  max_prim <- rep(0, nr_of_families)
  max_sec <- rep(0, nr_of_families)
  max_ter <- rep(0, nr_of_families)
  
  if (type == "time"){
    for (i in 1:nrow(prim_parameters)){
      max_prim[prim_parameters$fam_nr[i]] <- max_prim[prim_parameters$fam_nr[i]] + get_famsize(prim_parameters[i,], timepoint)
    }
    for (i in 1:nrow(sec_parameters)){
      max_sec[sec_parameters$fam_nr[i]] <- max_sec[sec_parameters$fam_nr[i]] + get_famsize(sec_parameters[i,], timepoint)
    }
    for (i in 1:nrow(ter_parameters)){
      max_ter[ter_parameters$fam_nr[i]] <- max_ter[ter_parameters$fam_nr[i]] + get_famsize(ter_parameters[i,], timepoint)
    }
  }
  if (type == "max"){
    for (i in 1:nrow(prim_parameters)){
      max_prim[prim_parameters$fam_nr[i]] <- max_prim[prim_parameters$fam_nr[i]] + get_max_famsize(prim_parameters[i,])
    }
    for (i in 1:nrow(sec_parameters)){
      max_sec[sec_parameters$fam_nr[i]] <- max_sec[sec_parameters$fam_nr[i]] + get_max_famsize(sec_parameters[i,])
    }
    for (i in 1:nrow(ter_parameters)){
      max_ter[ter_parameters$fam_nr[i]] <- max_ter[ter_parameters$fam_nr[i]] + get_max_famsize(ter_parameters[i,])
    }
  }
  
  nr_burst_divs <- prim_parameters$nr_burst_divs[match(1:nr_of_families, prim_parameters$fam_nr)]
  
  max_cells <- data.frame(
    fam_nr = seq(1, nr_of_families),
    cells_prim = max_prim,
    cells_sec = max_sec,
    cells_ter = max_ter,
    nr_burst_divs = factor(nr_burst_divs)
  )
  if (remove_zeroes == TRUE) {
    max_cells <- max_cells[max_cells$cells_ter > 0, ]
  }
  
  return(max_cells)
}

#-------------------------------------------------------------------------------

generate_max_fam_stats <- function(max_cells){
  # Calculate Spearman correlation coefficients
  cor_prim_sec <- cor.test(max_cells$cells_prim, max_cells$cells_sec, method = "spearman")
  cor_sec_ter <- cor.test(max_cells$cells_sec, max_cells$cells_ter, method = "spearman")
  
  # Calculate mean and median fam size
  mean_prim <- mean(max_cells$cells_prim)
  mean_sec <- mean(max_cells$cells_sec)
  mean_ter <- mean(max_cells$cells_ter)
  median_prim <- median(max_cells$cells_prim)
  median_sec <- median(max_cells$cells_sec)
  median_ter <- median(max_cells$cells_ter)
  stats <- data.frame(mean = c(mean_prim, mean_sec, mean_ter),
                      median = c(median_prim, median_sec, median_ter),
                      correlation = c(cor_prim_sec$estimate, 0, cor_sec_ter$estimate),
                      p_value = c(cor_prim_sec$p.value, 0, cor_sec_ter$p.value))
  return(stats)
}

#-------------------------------------------------------------------------------

plot_prim_sec_max_fam <- function(max_cells, 
                         show_title = TRUE){ 
  
  plot <- ggplot(data = max_cells, 
                          aes(x = log10(cells_prim), y = log10(cells_sec))) +
    geom_point(aes(col = nr_burst_divs)) +
    labs(title = "Primary vs. Secondary Response",
         x = "Family size primary (log-scale)",
         y = "Family size secondary (log-scale)") +
    theme_clean() +
    theme(legend.position = "none")
  

  
  if (show_title == FALSE) {
    plot <- plot + labs(title = NULL, subtitle = NULL)
  }
  
  return(plot)
} 

#-------------------------------------------------------------------------------

plot_sec_ter_max_fam <- function(max_cells, 
                                show_title = TRUE){ 

  plot <- ggplot(data = max_cells, 
                         aes(x = log10(cells_sec), y = log10(cells_ter))) +
    geom_point(aes(col = nr_burst_divs)) +
    labs(
      title = "Secondary vs. Tertiary Response",
      x = "Family size secondary (log-scale)",
      y = "Family size tertiary (log-scale)") + 
    theme_clean() +
    theme(legend.position = "none")
  
  
  if (show_title == FALSE) {
    plot <- plot + labs(title = NULL, subtitle = NULL)
  }
  
  return(plot)
} 

#-------------------------------------------------------------------------------

generate_famsize_table <- function(parameters, timepoint = NULL) {
  fam_nr <- c()
  famsize <- c()
  
  for (i in 1:nrow(parameters)) {
    fam_nr[i] <- parameters$fam_nr[i]
    famsize[i] <- get_famsize(parameters[i,], timepoint)
  }

  df_famsizes <- data.frame(fam_nr = fam_nr, famsize = famsize)
  df_famsizes <- aggregate(famsize ~ fam_nr, data = df_famsizes, sum) 
  df_famsizes <- df_famsizes[order(df_famsizes$famsize, decreasing = TRUE), ]
  df_famsizes$logfamsize <- log2(df_famsizes$famsize)
  
  return(df_famsizes)
}

#-------------------------------------------------------------------------------

generate_famsize_stats <- function(df_famsizes){
  stat_mean <- round(mean(df_famsizes$famsize[df_famsizes$famsize > 1]), digits = 0)
  stat_median <- round(median(df_famsizes$famsize[df_famsizes$famsize > 1]), digits = 0)
  
  return(c(stat_mean, stat_median))
}

#-------------------------------------------------------------------------------

generate_freq_famsize_table <- function(logfamsizes){
  freq_famsizes <- data.frame(table(round(logfamsizes)))
  colnames(freq_famsizes) <- c("logfamsize", "freq")
  freq_famsizes$freq <- freq_famsizes$freq/sum(freq_famsizes$freq)
  while (nrow(freq_famsizes) < 16) {
    freq_famsizes <- add_row(freq_famsizes, logfamsize = as.factor(nrow(freq_famsizes)+1), freq = 0)}
  
  return(freq_famsizes)
}

#-------------------------------------------------------------------------------

plot_famsize_distribution <- function(freq_famsizes, timepoint, show_title = F){

  if(timepoint %in% c(5,6,7,8)){
    expected_famsizes <- get(paste0("day", timepoint))}
  else{expected_famsizes <- data.frame(Number_2log = 0, Frequency = 0)}
  # To make sure the function doesnt break when the timepoint is not 5,6,7 or 8
  
  plot <- ggplot(data = freq_famsizes, aes(x = logfamsize)) + 
    geom_bar(aes(y = freq), stat = "identity", fill = "skyblue") +
    geom_point(data = expected_famsizes, aes(x = Number_2log, y = Frequency), color = "darkorange", size = 2) + 
    labs(title = paste("Day", timepoint),
         x = "Family size (2log)",
         y = "Frequency") +
    theme(plot.margin = margin(t = 1, r = 1, 0, 0, unit = "cm")) +
    theme_clean() 
    
  if (show_title == F) {
    plot <- plot + labs(title = NULL, subtitle = NULL)
  }
    
  return(plot)
}

#-------------------------------------------------------------------------------

plot_cum_famsize <- function(df_famsizes, timepoint, show_title = F){
  df_famsizes$nr <- seq(nrow(df_famsizes))
  
  zeros_row <- c(0, 0, 0, 0) 
  # We would like to have this row as the first value, so that the cumplot 
  # starts at 0,0
  df_famsizes <- rbind(zeros_row, df_famsizes)

  plot <- ggplot(data = df_famsizes,
                     aes(x = nr / max(nr) * 100,
                         y = cumsum(famsize) / sum(famsize) * 100)) +
    geom_point() +
    labs(title = paste("Cumulative family sizes day", timepoint),
         x = "% of cumulative size",
         y = "% of accumulated progenies") +
    scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    geom_vline(xintercept = 5, colour = "red") +
    theme_clean()
  
  if (show_title == F) {
    plot <- plot + labs(title = NULL, subtitle = NULL)
  }
  return(plot)
}

#-------------------------------------------------------------------------------

generate_famsize_dist_multidays <- function(parameters){
  
  famsizes_table <- lapply(5:8, function(day) {
    generate_famsize_table(parameters, timepoint = day)
  })
  
  return(famsizes_table)
} 

#-------------------------------------------------------------------------------

generate_famsize_stats_multidays <- function(famsizes_table){
  famsizes_stats <- lapply(5:8, function(day) {
    c(day, generate_famsize_stats(famsizes_table[[day-4]]))
  })
  famsizes_stats <- data.frame(do.call(rbind, famsizes_stats))
  colnames(famsizes_stats) <- c("timepoint", "mean", "median")
  
  return(famsizes_stats)
}

#-------------------------------------------------------------------------------

plot_grid_famsize_dist <- function(famsizes_table){
  frequencies_table <- lapply(1:4, function(day) {
    generate_freq_famsize_table(famsizes_table[[day]]$logfamsize)
  })
  plots <- lapply(5:8, function(day) {
    plot_famsize_distribution(frequencies_table[[day-4]], day, show_title = T)
  })
  plot <- grid.arrange(grobs = plots, ncol = 2)
  
  return(plot)
}

# _______________________________________________________________________________
# This function plots the size of the family as a function of the number of Q cells that are formed

generate_Q_famsize_table <- function(parameters) {
  max_famsizes <- get_max_famsize(parameters)
  Q_famsize_table <- data.frame(fam_nr = parameters$fam_nr, 
                                cell_type = parameters$cell_type, 
                                nr_burst_divs = parameters$nr_burst_divs, 
                                max_famsize = max_famsizes)
  Q_famsize_table %<>%
    group_by(fam_nr) %>%
    summarize(
      famsize = sum(max_famsize),
      Q_cells = sum(cell_type == "Q"),
      nr_burst_divs = as.factor(first(nr_burst_divs))  # Since nr_burst_divs is the same for each fam_nr
    )
  
  Q_famsize_table$log_famsize <- log10(Q_famsize_table$famsize)
  
  return(Q_famsize_table)
}

#-------------------------------------------------------------------------------

get_Q_famsize_stats <- function(Q_famsize_table){
  correlation <- cor.test(Q_famsize$famsize, Q_famsize$Q_cells, method = "spearman")
  
  return(correlation)
}

#-------------------------------------------------------------------------------

plot_Q_famsize <- function(Q_famsize_table, show_title = T, show_legend = T){
  plot <- ggplot(data = Q_famsize) +
    geom_point(aes(x = log_famsize, y = Q_cells, col = burst_divs)) +
    labs(
      x = "10log(max. family size)",
      y = "Nr. of Q cells",
      color = "Number of burst divisions"
    ) + theme_clean() +
    theme(legend.position = "top",
          legend.justification = 0,
          legend.background = element_rect(fill="white",
                                           size=0.5, linetype="solid", 
                                           colour ="black")) 
  
  if (show_title == F) {
    plot <- plot + labs(title = NULL, subtitle = NULL)
  }
  if (show_legend == F) {
    plot <- plot + theme(legend.position = "none")
  }
  
  return(plot)
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Functions below are not used anymore

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

wide_to_long <- function(solve_table) {
  solve_table %<>% pivot_longer(cols = 4:ncol(solve_table), 
                                names_to = "time", 
                                values_to = "cells")
  solve_table$time <- as.numeric(solve_table$time)
  return(solve_table)
}


#-------------------------------------------------------------------------------

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
