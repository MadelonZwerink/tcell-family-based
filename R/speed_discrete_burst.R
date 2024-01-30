speed_version <- "29-08-2023"

path_famsizes <- "C:/Users/madel/OneDrive/Documenten/BiBC/Major_internship/Literature/Data_distribution_famsizes/data/"

library(tidyverse)
library(magrittr)
library(ggplot2)
library(purrr)
library(data.table)
library(cowplot)

source(paste0(path_famsizes, "data load.R"))
#-------------------------------------------------------------------------------

pick_parameters <- function(
    bp_rule = "runif(1, min = 0.5, max = 3)",
    d_p,
    sd_dp = 0,
    r_q = NULL, 
    b_q = 0,
    sd_bq = 0,
    d_q = 0,
    sd_dq = 0,
    t_start_dist, 
    t_run_rule,    
    nr_of_families = NULL, 
    nr_burst_divs = 3,
    correction_t_burst = NULL,
    response_nr = 1,
    prev_parameters = NULL,
    primary_parameters = NULL,
    secondary_parameters = NULL,
    quality_dist = NULL,
    Q1 = NULL,
    Q2 = NULL,
    Q3 = NULL,
    Q4 = NULL,
    ASD = F) {
  
  # for (name in c("mean_t_start", "sd_t_start", "mean_t_run", "sd_t_run", "mean_b", "sd_b", "nr_of_families")){
  #  overview_params[name, response_nr] <<- eval(parse(text = name))
  # }
  max_run_time <- 8
  
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
                           dp = numeric(),
                           bq = numeric(),
                           dq = numeric())
  
  
  if(is.null(prev_parameters) == FALSE){
    prev_Q <<- prev_parameters[prev_parameters$cell_type == "Q",]
    nr_of_families <- nrow(prev_Q)
    families <- nr_of_families
    prev_divs <- prev_Q$div_counter
  }
  else {prev_divs <- rep(0, nr_of_families)}
  
  if(!is.null(quality_dist)){
    quality <- round(eval(parse(text = paste(quality_dist))), digits = 2)
      # beta distribution: rbeta(nr_of_families, shape1 = 2, shape2 = 5)
      # uniform distribution: runif(nr_of_families, 0, 1)
  } else {quality <- rep(NA, nr_of_families)}
  
  ifelse(!is.null(Q1) & !is.null(quality_dist), Q_1 <- eval(parse(text = Q1)), Q_1 <- rep(1, nr_of_families))            #rq: max fraction (min is 0) 

  nr_burst_divs <- eval(parse(text = nr_burst_divs))
  t_burst <- rep(0.2, nr_of_families)    
  t_start_expr <- substitute(eval(parse(text = t_start_dist)))
  t_start_val <- eval(t_start_expr)
  t_start <- round(t_start_val, digits =2)
    
  for (i in 1:nr_of_families){    
    #Make sure the families have the correct family number
    if(response_nr == 1){fam_nr <- i
                        fam_nr_2 <- 0
                        fam_nr_3 <- 0}
    if(response_nr == 2){fam_nr <- prev_Q$fam_nr[i]
                        fam_nr_2 <- i
                        fam_nr_3 <- 0}
    if(response_nr == 3){fam_nr <- prev_Q$fam_nr[i]
                        fam_nr_2 <- prev_Q$fam_nr_2[i]
                        fam_nr_3 <- i}
    
    div_counter = nr_burst_divs[i] + prev_divs[i]
    
    if(!is.null(quality_dist) & !is.null(Q4)){t_start[i] <- eval(parse(text = Q4))}
    
    frac_rq <- ifelse(is.null(r_q), runif(1, min = 0, max = Q_1[i]), eval(parse(text = r_q)))
    q_cells <- ifelse(ASD == F, min(max(ceiling(frac_rq * 2^nr_burst_divs[i]), 0), 2^nr_burst_divs[i] ), 1)
    p_cells <- (2^nr_burst_divs[i]) - q_cells

    if(p_cells != 0){
      for (cell in 1 : p_cells){
        #rnorm(1, mean = ifelse(quality_parameter == T, Q_3[i], mean_t_run), sd = sd_t_run)
        
        #First define the expression and then evaluate the expression
        #If tbis is not included, the expression is only evaluated once, leading to the same t_run for all (sub)families

        t_run_expr <- substitute(eval(parse(text = t_run_rule)))
        t_run_val <- eval(t_run_expr)
        t_run <- round(t_run_val, digits = 2)
        bp <- round(eval(parse(text = bp_rule)), digits = 2)
        
        #this is to make sure all proliferation stops at t=7
        if (t_start[i] + t_run >= max_run_time){
          t_correct <- round((-max_run_time + t_start[i] + t_run), digits = 2)
        } else { t_correct <- 0 }
        if (t_start[i] < 0){
          t_start[i] = 0.01}
        if (t_start[i] > 6.5){
          t_start[i] <- 6.5}
        if (t_start[i] < 2){
          t_start[i] <- 2.5}

        if(!is.null(quality_dist)){
          if(!is.null(Q2)){
            bp <- round(eval(parse(text = Q2)), digits = 2)
          }
          if(!is.null(Q3)){
            t_run <- round(eval(parse(text = Q3)), digits = 2)
          }
        }
        #bp <- round(rnorm(1, mean = b_p, sd = sd_bp), digits = 2)
        dp <- round(rnorm(1, mean = d_p, sd = sd_dp), digits = 2)
        bq <- round(rnorm(1, mean = eval(b_q), sd = sd_bq), digits = 2)
        dq <- round(rnorm(1, mean = eval(d_q), sd = sd_dq), digits = 2)

        
        parameters %<>% add_row(cell_type = "P",
                                div_counter = div_counter,
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
                                dp = dp,
                                bq = bq,
                                dq = dq)
      }}
    
    if(q_cells != 0){
      for (cell in 1 : q_cells){
        t_run <- 0
        #this is to make sure all proliferation stops at t=7
        if (t_start[i] + t_run  >= max_run_time){
          t_correct <- round((-max_run_time + t_start[i] + t_run), digits = 2)} 
        else { t_correct <- 0 }
        if (t_start[i] < 0){
          t_start[i] = 0.01}
        if (t_start[i] > 6.5){
          t_start[i] <- 6.5}
        
        bp <- 0
        dp <- 0
        bq <- round(rnorm(1, mean = eval(b_q), sd = sd_bq), digits = 2)
        dq <- round(rnorm(1, mean = eval(d_q), sd = sd_dq), digits = 2)
        
        parameters %<>% add_row(cell_type = "Q",
                                div_counter = div_counter,
                                t_start = t_start[i], 
                                quality = quality[i],
                                fam_nr = fam_nr,
                                fam_nr_2 = fam_nr_2,
                                fam_nr_3 = fam_nr_3,
                                nr_burst_divs = nr_burst_divs[i],
                                t_burst = t_burst[i],
                                t_run = 0, 
                                t_correct = t_correct,
                                bp = bp,
                                dp = dp,
                                bq = bq,
                                dq = dq)
    }}
    
  }
  return(parameters)
}

#-------------------------------------------------------------------------------

#this is to make sure all proliferation stops at t=10
#if (t_start[i] + t_run >= 8){
 # t_correct <- round((-8 + t_start[i] + t_run), digits = 2)
#} else { t_correct <- 0 }
#if (t_start[i] < 0){
 # t_start[i] = 0.01
#}
#-------------------------------------------------------------------------------

run_model <- function(parameters, 
                      t_total = 100, 
                      response_nr = 1, 
                      t_previous = 0, 
                      prev_results = NULL) {
  n_sol = list()
  parameters <<- parameters
  
  for (i in 1 : nrow(parameters)){ 
    t_run <<- parameters$t_run[i]
    t_correct <<- parameters$t_correct[i]
    
    input_vectors <- generate_input_vectors(n = n, 
                                            parameters = parameters, 
                                            i = i)
    #make sure that cells do not proliferate after tRun
    this_run <- run(tmax = t_total - parameters$t_start[i] - 0.01, 
                    tstep = 0.01, 
                    state = input_vectors, 
                    parms = c(dp = parameters$dp[i],
                              bp = parameters$bp[i],
                              dq = parameters$dq[i],
                              bq = parameters$bq[i]),
                    table = T, 
                    arrest = (t_run - t_correct),
                    after = "if(t==t_run - t_correct)parms[\"bp\"]<-0", 
                    timeplot = F)
    
    #add time before the start of replication to the time
    this_run$time <- this_run$time + parameters$t_start[i] 
    
    #add the values before start of proliferation and set to input vectors
    for (time_point in seq(0, parameters$t_start[i] - 0.005, 0.01)) {
      this_run[nrow(this_run) + 1, ] <- c(time_point, input_vectors)
    }
    
    this_run$time <- round(this_run$time, digits=2) + t_previous + 0.01
    this_run <- this_run[order(this_run$time), ]
    this_run %<>% mutate(fam_nr = parameters$fam_nr[i],
                         fam_nr_2 = 0,
                         fam_nr_3 = 0)
    
    if (response_nr == 2){
      this_run %<>% mutate(fam_nr_2 = i)
    }
    if (response_nr == 3){
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
                                   i){
  P <- rep(0, n)
  Q <- rep(0, n)
  if(parameters$cell_type[i] == "P"){
    P[parameters$div_counter[i]] <- 1
  } else {
    Q[parameters$div_counter[i]] <- 1
  }
  names(P) <- paste0("P", seq(1, n))
  names(Q) <- paste0("Q", seq(1, n))
  return(c(P, Q))
}

#-------------------------------------------------------------------------------

solve_function <- function(parameters, min_time, max_time, interval) {
  num_timepoints <- length(seq(min_time, max_time, interval)) + 3
  P_total <- data.frame(matrix(0, nrow = nrow(parameters), ncol = num_timepoints))
  colnames(P_total) <- c("fam_nr", "fam_nr_2", "fam_nr_3", seq(min_time, max_time, by = interval))
  
  for (cell in 1:nrow(parameters)) {
    P_total$fam_nr[cell] <- parameters$fam_nr[cell]
    P_total$fam_nr_2[cell] <- parameters$fam_nr_2[cell]
    P_total$fam_nr_3[cell] <- parameters$fam_nr_3[cell]
    t_max <- round(min_time + parameters$t_run[cell] - parameters$t_correct[cell], digits = 2)
    max_cells <- 0
    burst_time <- parameters$t_burst[cell] * parameters$nr_burst_divs[cell]
    
    for (timepoint in seq(min_time, max_time, by = interval)) {
      
      if(timepoint <= parameters$t_start[cell]){ #before family starts dividing
        P_total[cell, col_index] <- 0 
      } else if (timepoint <= burst_time){ #while family is in burst division
        burst_divs_v <- burst_divs_vector[1:parameters$nr_burst_divs[cell]]
        div_times <- burst_divs_v * parameters$t_burst[cell] #creates a vector with the times where a cell divides
        differences <- div_times - timepoint
        index <- which(min(differences[differences >= 0]) == differences)
        P_total[cell, col_index] <- 2^burst_divs_v[index]
      } else if(time[i] <= parameters$t_run[i] - parameters$t_correct[i] & time[i] > burst_time[i]){ #while family is dividing
        famsizes[parameters$fam_nr[i]] <- famsizes[parameters$fam_nr[i]] + exp((parameters$bp[i]) * (time[i] - burst_time[i]))
      } else { #after family stopped dividing
        max_cells <- exp((parameters$bp[i]) * (parameters$t_run[i] - parameters$t_correct[i]))
        after_contraction <- max_cells * exp(- 1 * parameters$dp[i] * (time[i] - parameters$t_run[i] - burst_time[i] + parameters$t_correct[i]))
        famsizes[parameters$fam_nr[i]] <- famsizes[parameters$fam_nr[i]] + after_contraction
      }
      
      
      timepoint <- round(timepoint, digits = 2)
      cor_time <- timepoint + parameters$t_start[cell]
      col_index <- which(colnames(P_total) == cor_time)
      if (abs(timepoint - t_max) < (interval/2)){
        max_cells <- exp(parameters$bp[cell] * timepoint)
      } else if (timepoint < t_max){
        P_total[cell, col_index] <- exp(parameters$bp[cell] * timepoint)
      } else if (timepoint > t_max){
        P_total[cell, col_index] <- max_cells * exp(- 1 * parameters$dp[cell] * (timepoint - parameters$t_run[cell]))
      }
    }
    
    if (max_cells == 0){
      print(paste("Warning: cell", cell, "from family", parameters$fam_nr[cell], ".", 
                  parameters$fam_nr_2[cell], ".", parameters$fam_nr_3[cell], "skipped statement to set max_cells"))
    }
  }
  P_total$fam_nr <- as.factor(P_total$fam_nr)
  return(P_total)
}

#-------------------------------------------------------------------------------

wide_to_long <- function(solve_table){
  solve_table %<>% pivot_longer(cols = 4 : ncol(solve_table), names_to = "time", values_to = "cells")
  solve_table$time <- as.numeric(solve_table$time)
  return(solve_table)
}

#-------------------------------------------------------------------------------

generate_max_fam_plot <- function(prim_parameters, sec_parameters, ter_parameters, type = c("time", "max"), time = NULL, remove_zeroes = FALSE, relative = FALSE, show_title = TRUE, show_legend = TRUE){
  library(cowplot)
  
  nr_of_families <- max(prim_parameters$fam_nr)
  max_prim <- rep(0, nr_of_families)
  max_sec <- rep(0, nr_of_families)
  max_ter <- rep(0, nr_of_families)
  
  if(type == "time"){
    time_prim <- min(time - prim_parameters$t_start, prim_parameters$t_run - prim_parameters$t_correct)
    time_sec <- min(time - sec_parameters$t_start, sec_parameters$t_run - sec_parameters$t_correct)
    time_ter <- min(time - ter_parameters$t_start, ter_parameters$t_run - ter_parameters$t_correct)
  }
  if(type == "max"){
    time_prim <- prim_parameters$t_run - prim_parameters$t_correct
    time_sec <- sec_parameters$t_run - sec_parameters$t_correct
    time_ter <- ter_parameters$t_run - ter_parameters$t_correct
  }
  for(i in 1:nrow(prim_parameters)){
    if(prim_parameters$t_run[i] - prim_parameters$t_correct[i] >= time_prim[i]){
      max_prim[prim_parameters$fam_nr[i]] <- max_prim[prim_parameters$fam_nr[i]] + exp((prim_parameters$bp[i]) * time_prim[i])
    } else {
      max_cells <- exp((prim_parameters$bp[i]) * prim_parameters$t_run[i])
      after_contraction <- max_cells * exp(- 1 * prim_parameters$dp[i] * (time_prim[i] - prim_parameters$t_run[i]))
      max_prim[prim_parameters$fam_nr[i]] <- max_prim[prim_parameters$fam_nr[i]] + after_contraction
    }
  }
  for(i in 1:nrow(sec_parameters)){
    if(sec_parameters$t_run[i] - sec_parameters$t_correct[i] >= time_sec[i]){
      max_sec[sec_parameters$fam_nr[i]] <- max_sec[sec_parameters$fam_nr[i]] + exp((sec_parameters$bp[i]) * time_sec[i])
    } else {
      max_cells <- exp((sec_parameters$bp[i]) * sec_parameters$t_run[i])
      after_contraction <- max_cells * exp(- 1 * sec_parameters$dp[i] * (time_sec[i] - sec_parameters$t_run[i]))
      max_sec[sec_parameters$fam_nr[i]] <- max_sec[sec_parameters$fam_nr[i]] + after_contraction
    }
  }
  for(i in 1:nrow(ter_parameters)){
    if(ter_parameters$t_run[i] - ter_parameters$t_correct[i] >= time_ter[i]){
      max_ter[ter_parameters$fam_nr[i]] <- max_ter[ter_parameters$fam_nr[i]] + exp((ter_parameters$bp[i]) * time_ter[i])
    } else {
      max_cells <- exp((ter_parameters$bp[i]) * ter_parameters$t_run[i])
      after_contraction <- max_cells * exp(- 1 * ter_parameters$dp[i] * (time_ter[i] - ter_parameters$t_run[i]))
      max_ter[ter_parameters$fam_nr[i]] <- max_ter[ter_parameters$fam_nr[i]] + after_contraction
    }
  }
  nr_burst_divs <- c()
  for(i in 1:nr_of_families){
    nr_burst_divs[i] <- prim_parameters$nr_burst_divs[first(which(prim_parameters$fam_nr == i))]
  }
  max_cells <- data.frame(fam_nr = seq(1, nr_of_families),
                          cells_prim = max_prim, 
                          cells_sec = max_sec,
                          cells_ter = max_ter,
                          nr_burst_divs = factor(nr_burst_divs),
                          rel_cells_prim = max_prim/sum(max_prim),
                          rel_cells_sec = max_sec/sum(max_sec),
                          rel_cells_ter = max_ter/sum(max_ter))
  if(remove_zeroes == TRUE){
    max_cells <- max_cells[max_cells$cells_ter > 0,]
  }
  
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
                      median = c(median_prim, median_sec, median_ter))
  
  plot_prim_sec <- ggplot(data = max_cells, aes(x=log10(cells_prim), y=log10(cells_sec))) + 
    geom_point(aes(col = nr_burst_divs)) +
    labs(title = "Primary vs. Secondary Response",
         subtitle = paste("Nr. of families:", nr_of_families, 
                          "\nMethod:", type, 
                          if(type=="time"){paste( "\nTime:", time)}),
         x = "Primary Max Value",
         y = "Secondary Max Value") +
    annotate("text", x = log10(min(max_cells$cells_prim)), 
             y = log10(min(max_cells$cells_sec)) + 0.4,
             hjust = 0,
             label = paste("r =", round(cor_prim_sec$estimate, digits = 2), 
                           "\nP =", format(cor_prim_sec$p.value, scientific = TRUE, digits = 2),
                           "\nTotal primary:", format(sum(max_prim), scientific = TRUE, digits = 2),
                           "\nTotal secondary:", format(sum(max_sec), scientific = TRUE, digits = 2)))
  legend <- get_legend(
    # create some space to the left of the legend
    plot_prim_sec + theme(legend.box.margin = margin(0, 0, 0, 12))
  )
  
  plot_prim_sec <- plot_prim_sec +
    theme(legend.position = "none")
  
  plot_sec_ter <- ggplot(data = max_cells, aes(log10(cells_sec), y= log10(cells_ter))) + 
    geom_point(aes(col = nr_burst_divs)) +
    labs(title = "Secondary vs. Tertiary Response",
         subtitle = paste("Nr. of families:", nr_of_families, 
                          "\nMethod:", type, 
                          if(type=="time"){paste( "\nTime:", time)}),
         x = "Secondary Max Value",
         y = "Tertiary Max Value") +
    annotate("text", x = log10(min(max_cells$cells_sec)), 
             y = log10(min(max_cells$cells_ter)) + 0.4,
             hjust = 0,
             label = paste("r =", round(cor_sec_ter$estimate, digits = 2), 
                           "\nP =", format(cor_sec_ter$p.value, scientific = TRUE, digits = 2),
                           "\nTotal primary:", format(sum(max_sec), scientific = TRUE, digits = 2),
                           "\nTotal secondary:", format(sum(max_ter), scientific = TRUE, digits = 2))) +
    theme(legend.position = "none")
  
  plot_rel_prim_sec <- ggplot(data = max_cells, aes(x=log10(rel_cells_prim), y=log10(rel_cells_sec))) + 
    geom_point(aes(col = nr_burst_divs)) +
    labs(title = "Primary vs. Secondary Response",
         subtitle = paste("Nr. of families:", nr_of_families, 
                          "\nMethod:", type, 
                          if(type=="time"){paste( "\nTime:", time)}),
         x = "Primary Max Value",
         y = "Secondary Max Value") +
    annotate("text", x = log10(min(max_cells$rel_cells_prim)), 
             y = log10(max(max_cells$rel_cells_sec)) - 0.3,
             hjust = 0,
             label = paste("r =", round(cor_prim_sec$estimate, digits = 2), 
                           "\nP =", format(cor_prim_sec$p.value, scientific = TRUE, digits = 2),
                           "\nTotal primary:", format(sum(max_prim), scientific = TRUE, digits = 2),
                           "\nTotal secondary:", format(sum(max_sec), scientific = TRUE, digits = 2))) +
    theme(legend.position = "none")
  
  plot_rel_sec_ter <- ggplot(data = max_cells, aes(x=log10(rel_cells_sec), y=log10(rel_cells_ter))) + 
    geom_point(aes(col = nr_burst_divs)) +
    labs(title = "Secondary vs. Tertiary Response",
         subtitle = paste("Nr. of families:", nr_of_families, 
                          "\nMethod:", type, 
                          if(type=="time"){paste( "\nTime:", time)}),
         x = "Secondary Max Value",
         y = "Tertiary Max Value") +
    annotate("text", x = log10(min(max_cells$rel_cells_sec)), 
             y = log10(max(max_cells$rel_cells_ter)) - 0.3,
             hjust = 0,
             label = paste("r =", round(cor_sec_ter$estimate, digits = 2), 
                           "\nP =", format(cor_sec_ter$p.value, scientific = TRUE, digits = 2),
                           "\nTotal secondary:", format(sum(max_sec), scientific = TRUE, digits = 2),
                           "\nTotal tertiary:", format(sum(max_ter), scientific = TRUE, digits = 2))) +
    theme(legend.position = "none")
  
  if(show_title == FALSE){
    plot_prim_sec <- plot_prim_sec + labs(title = NULL, subtitle = NULL)
    plot_sec_ter <- plot_sec_ter + labs(title = NULL, subtitle = NULL)
    plot_rel_prim_sec <- plot_rel_prim_sec + labs(title = NULL, subtitle = NULL)
    plot_rel_sec_ter <- plot_rel_sec_ter + labs(title = NULL, subtitle = NULL)
  }
  
  if(relative == T){
    plots <- plot_grid(plot_rel_prim_sec, plot_rel_sec_ter)
  } else {plots <- plot_grid(plot_prim_sec, plot_sec_ter)}
  
  if(show_legend == T){
    plots <- plot_grid(plots, legend, rel_widths = c(4, 0.5)) 
  }
  
  return(list(max_cells, plots, stats, cor_prim_sec, cor_sec_ter, legend))
}

#-------------------------------------------------------------------------------

plot_famsize_dist <- function(parameters, timepoint = NULL, method = c("time", "max"), binwidth = 0.2, show_title = T){
  famsizes <- rep(0, max(parameters$fam_nr))
  burst_time <- parameters$t_burst * parameters$nr_burst_divs
  burst_divs_vector <- c(1,2,3,4,5,6,7)
  
  if(method == "time"){
    time <- round(timepoint - parameters$t_start, digits = 2)
  }
  if(method == "max"){
    time <- parameters$t_run - parameters$t_correct
  }
  for(i in 1:nrow(parameters)){
    if(time[i] <= 0){ #before family starts dividing
      famsizes[parameters$fam_nr[i]] <- 1 
    } else if (time[i] <= burst_time[i]){ #while family is in burst division
      burst_divs_v <- burst_divs_vector[1:parameters$nr_burst_divs[i]]
      div_times <- burst_divs_v * parameters$t_burst[i] #creates a vector with the times where a cell divides
      differences <- div_times - time[i]
      index <- which(min(differences[differences >= 0]) == differences)
      famsizes[parameters$fam_nr[i]] <- 2^burst_divs_v[index]
    } else if(time[i] <= (parameters$t_run[i] - parameters$t_correct[i]) & time[i] > burst_time[i]){ #while family is dividing
      famsizes[parameters$fam_nr[i]] <- famsizes[parameters$fam_nr[i]] + exp((parameters$bp[i]) * (time[i] - burst_time[i]))
    } else { #after family stopped dividing
      max_cells <- exp((parameters$bp[i]) * (parameters$t_run[i] - parameters$t_correct[i]))
      after_contraction <- max_cells * exp(- 1 * parameters$dp[i] * (time[i] - parameters$t_run[i] + parameters$t_correct[i]))
      famsizes[parameters$fam_nr[i]] <- famsizes[parameters$fam_nr[i]] + after_contraction
    }
  }
  mean_val <- mean(famsizes)
  ci_lower <- mean_val - 1.96 * sd(famsizes) / sqrt(length(famsizes))
  ci_upper <- mean_val + 1.96 * sd(famsizes) / sqrt(length(famsizes))
  famsizes <- sort(famsizes, decreasing = T)
  df_famsizes <- data.frame(nr = seq(0, length(famsizes)), famsizes = c(0, famsizes), logfamsizes = c(0, log2(famsizes)))
  expected_famsizes_day5 <- data.frame(Amount = day5$Number_2log, Freq = day5$Freq * (nrow(df_famsizes) - 1))
  expected_famsizes_day6 <- data.frame(Amount = day6$Number_2log, Freq = day6$Freq * (nrow(df_famsizes) - 1))
  expected_famsizes_day7 <- data.frame(Amount = day7$Number_2log, Freq = day7$Freq * (nrow(df_famsizes) - 1))
  expected_famsizes_day8 <- data.frame(Amount = day8$Number_2log, Freq = day8$Freq * (nrow(df_famsizes) - 1))

  plot_distribution <- ggplot(data = df_famsizes, aes(logfamsizes)) +
    geom_histogram(binwidth = binwidth)  +
    geom_point(data = expected_famsizes_day5, aes(x = Amount, y = Freq), color = "green", size = 2) +  
    geom_point(data = expected_famsizes_day6, aes(x = Amount, y = Freq), color = "yellow", size = 2) + 
    geom_point(data = expected_famsizes_day7, aes(x = Amount, y = Freq), color = "orange", size = 2) +  
    geom_point(data = expected_famsizes_day8, aes(x = Amount, y = Freq), color = "red", size = 2) + 
    annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1,
             label = paste("Mean:", round(mean_val, digits = 1),
                           "\n95% CI:", round(ci_lower, digits = 1), "-", round(ci_upper, digits = 1),
                           "\nMedian:", round(median(df_famsizes$famsizes), digits = 1))) +
    labs(title = "Distribution of family sizes",
         subtitle = paste("Nr. of families:", nrow(df_famsizes) - 1, "\nTime: day", timepoint, "\nDay 5: blue, Day 6: red"),
         x = "Family size",
         y = "Frequency") +
    theme(plot.margin = margin(t = 1, r = 1, unit = "cm"))
  
  cum_plot <- ggplot(data = df_famsizes, aes(x = nr/max(nr)*100, y = cumsum(famsizes)/sum(famsizes) * 100)) + 
    geom_point() +
    labs(title = "Cumulative family sizes",
         subtitle = paste("Nr. of families:", nrow(df_famsizes) - 1, "\nTime: day", timepoint),
         x = "% of cumulative size",
         y = "% of accumulated progenies") +
    scale_x_continuous(limits = c(0,100), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0,100), expand = c(0, 0)) +
    geom_vline(xintercept = 5, colour = "red")
  
  if(show_title == F){
    plot_distribution <- plot_distribution+ labs(title = NULL, subtitle = NULL)
    cum_plot <- cum_plot + labs(title = NULL, subtitle = NULL)
  }
  
  expected_famsizes_day5$Amount <- as.factor(expected_famsizes_day5$Amount)
  expected_famsizes_day6$Amount <- as.factor(expected_famsizes_day6$Amount)
  expected_famsizes_day7$Amount <- as.factor(expected_famsizes_day7$Amount)
  expected_famsizes_day8$Amount <- as.factor(expected_famsizes_day8$Amount)
    
  return(list(plot_distribution, cum_plot, df_famsizes[2:nrow(df_famsizes),], expected_famsizes_day5, expected_famsizes_day6, expected_famsizes_day7, expected_famsizes_day8))
}

#_______________________________________________________________________________
#This function plots the size of the family as a function of the number of Q cells that are formed

plot_Q_famsize <- function(parameters, show_legend = T, show_title = T){
  famsize <- c()
  Q_cells <- c()
  burst_divs <- c()
  index_P_cells <- which(parameters$cell_type == "P")
  for(fam_nr in 1:max(parameters$fam_nr)){
    #For each family, first determine how many Q cells are formed
    Q_cells[fam_nr] <- sum(parameters$cell_type[parameters$fam_nr == fam_nr] == "Q")
    famsize[fam_nr] <- Q_cells[fam_nr]
    burst_divs[fam_nr] <- first(parameters$nr_burst_divs[parameters$fam_nr == fam_nr])
    
    #Then determine how many P cells are formed and calculate how large the progeny of each of these P cells will be
    fam_cells <- which(parameters$fam_nr == fam_nr)
    P_cells <- intersect(index_P_cells, fam_cells)
    
    for(P_cell in P_cells){
      famsize[fam_nr] <- famsize[fam_nr] + exp((parameters$bp[P_cell]) * parameters$t_run[P_cell])
    }
  }
  Q_famsize <- data.frame(Q_cells = Q_cells, burst_divs = as.factor(burst_divs), famsize = famsize, log_famsize = log10(famsize))   
  
  correlation <- cor.test(Q_famsize$Q_cells, Q_famsize$log_famsize)
  
  plot <- ggplot(data = Q_famsize) + geom_point(aes(x = Q_cells, y = log_famsize, col = burst_divs)) + 
    annotate("text", x = max(Q_cells) - 1, y = max(Q_famsize$log_famsize) - 0.5,
             label = paste("r =", round(correlation$estimate, digits = 2), 
                "\nP =", format(correlation$p.value, scientific = TRUE))) +
    labs(title = "Maximum family size as a function of nr. of Q cells",
         subtitle = paste("Nr. of families:", nrow(Q_famsize)),
         x = "10log(Max. family size)",
         y = "Nr. of Q cells")
  
  if(show_title == F){
    plot <- plot + labs(title = NULL, subtitle = NULL)
  }
  if(show_legend == F){
    plot <- plot + theme(legend.position = "none")
  }
  
  return(list(plot, correlation))
}

#-------------------------------------------------------------------------------

SSR_famsize_dist <- function(prim_parameters){
  exsmall_size <- 1
  pday5 <- plot_famsize_dist(prim_parameters, timepoint = 5, method = "time", binwidth = 1)
  pday6 <- plot_famsize_dist(prim_parameters, timepoint = 6, method = "time", binwidth = 1)
  pday7 <- plot_famsize_dist(prim_parameters, timepoint = 7, method = "time", binwidth = 1)
  pday8 <- plot_famsize_dist(prim_parameters, timepoint = 8, method = "time", binwidth = 1)
  
  model_famsizes_day5 <- as.data.frame(table(round(pday5[[3]]$logfamsizes)))
  names(model_famsizes_day5) <- c("Amount", "Freq")
  famsizes_day5 <- left_join(pday5[[4]], model_famsizes_day5, join_by(Amount), suffix = c("_Exp", "_Model"))
  famsizes_day5$Freq_Model[is.na(famsizes_day5$Freq_Model)] <- 0
  famsizes_day5_exsmall <- famsizes_day5[as.numeric(famsizes_day5$Amount) > exsmall_size,]
  
  famsizes_day5$SR <- (famsizes_day5$Freq_Exp - famsizes_day5$Freq_Model)^2
  SSR_day5 <- sum(famsizes_day5$SR)
  
  famsizes_day5_exsmall$SR <- (famsizes_day5_exsmall$Freq_Exp - famsizes_day5_exsmall$Freq_Model)^2
  SSR_day5_exsmall <- sum(famsizes_day5_exsmall$SR)
  
  model_famsizes_day6 <- as.data.frame(table(round(pday6[[3]]$logfamsizes)))
  names(model_famsizes_day6) <- c("Amount", "Freq")
  famsizes_day6 <- left_join(pday6[[5]], model_famsizes_day6, join_by(Amount), suffix = c("_Exp", "_Model"))
  famsizes_day6$Freq_Model[is.na(famsizes_day6$Freq_Model)] <- 0
  
  famsizes_day6_exsmall <- famsizes_day6[as.numeric(famsizes_day6$Amount) > exsmall_size,]
  
  famsizes_day6$SR <- (famsizes_day6$Freq_Exp - famsizes_day6$Freq_Model)^2
  SSR_day6 <- sum(famsizes_day6$SR)
  
  famsizes_day6_exsmall$SR <- (famsizes_day6_exsmall$Freq_Exp - famsizes_day6_exsmall$Freq_Model)^2
  SSR_day6_exsmall <- sum(famsizes_day6_exsmall$SR)
  
  model_famsizes_day7 <- as.data.frame(table(round(pday7[[3]]$logfamsizes)))
  names(model_famsizes_day7) <- c("Amount", "Freq")
  famsizes_day7 <- left_join(pday7[[6]], model_famsizes_day7, join_by(Amount), suffix = c("_Exp", "_Model"))
  famsizes_day7$Freq_Model[is.na(famsizes_day7$Freq_Model)] <- 0
  famsizes_day7_exsmall <- famsizes_day7[as.numeric(famsizes_day7$Amount) > exsmall_size,]
  
  famsizes_day7$SR <- (famsizes_day7$Freq_Exp - famsizes_day7$Freq_Model)^2
  SSR_day7 <- sum(famsizes_day7$SR)
  
  famsizes_day7_exsmall$SR <- (famsizes_day7_exsmall$Freq_Exp - famsizes_day7_exsmall$Freq_Model)^2
  SSR_day7_exsmall <- sum(famsizes_day7_exsmall$SR)
  
  model_famsizes_day8 <- as.data.frame(table(round(pday8[[3]]$logfamsizes)))
  names(model_famsizes_day8) <- c("Amount", "Freq")
  famsizes_day8 <- left_join(pday8[[7]], model_famsizes_day8, join_by(Amount), suffix = c("_Exp", "_Model"))
  famsizes_day8$Freq_Model[is.na(famsizes_day8$Freq_Model)] <- 0
  
  famsizes_day8_exsmall <- famsizes_day8[as.numeric(famsizes_day8$Amount) > exsmall_size,]
  
  famsizes_day8$SR <- (famsizes_day8$Freq_Exp - famsizes_day8$Freq_Model)^2
  SSR_day8 <- sum(famsizes_day8$SR)
  
  famsizes_day8_exsmall$SR <- (famsizes_day8_exsmall$Freq_Exp - famsizes_day8_exsmall$Freq_Model)^2
  SSR_day8_exsmall <- sum(famsizes_day8_exsmall$SR)
  
  return(list(list(SSR_day5, SSR_day6, SSR_day7, SSR_day8), 
              list(SSR_day5_exsmall, SSR_day6_exsmall, SSR_day7_exsmall, SSR_day8_exsmall), 
              list(famsizes_day5, famsizes_day6, famsizes_day7, famsizes_day8)))
}

