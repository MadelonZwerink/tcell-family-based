#!/usr/bin/env Rscript
speed_version <- "31-01-2024"

if(Sys.info()[[4]]=="LAPTOP-3RJSLMKV") {
  print("Working from R-project")
  setwd("~/BiBC/Major_internship/tcell-family-based")
} else {setwd('~/Documents')
        .libPaths("./R")}

packages <- c("tidyverse", "magrittr", "ggplot2", "purrr", "data.table", 
              "cowplot", "tibble", "ggthemes", "gridExtra", "patchwork",
              "ggpubr", "scales", "ggtext", "dplyr", "ggsignif")

# Function to install missing packages and log any errors
install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    tryCatch({
      install.packages(p, dependencies = TRUE, lib = "~/Documents/R")
    }, error = function(e) {
      message(sprintf("Error installing package '%s': %s", p, e$message))
    })
  }
  library(p, character.only = TRUE)
}

# Install and load each package
lapply(packages, install_if_missing)

# Verify the installation
installed_packages <- rownames(installed.packages(lib.loc = "./R"))
missing_packages <- setdiff(packages, installed_packages)
if (length(missing_packages) > 0) {
  message("The following packages could not be installed: ", paste(missing_packages, collapse = ", "))
} else {
  message("All packages installed successfully.")
}
#library(tidyverse)  # CRAN v2.0.0
#library(magrittr)   # CRAN v2.0.3
#library(ggplot2)    # CRAN v3.4.4
#library(purrr)      # CRAN v1.0.2
#library(data.table) # CRAN v1.15.0
#library(cowplot)    # CRAN v1.1.3
#library(tibble)     # CRAN v3.2.1
#library(ggthemes)
#library(gridExtra)
#library(patchwork)
#library(ggpubr)
#library(scales)
#library(ggtext)

source("./R/data_load.R") 
# The package plyr is loaded in the data_load file, therefore the package 
# has to be detached after loading the data
# Detach plyr if it is loaded otherwise group_by and summarize dont work 
# and collapse everything into one row
if ("package:plyr" %in% search()) {
  detach("package:plyr", unload=TRUE)
}
library(dplyr)

#-------------------------------------------------------------------------------

th <- theme(legend.position = "bottom",
            legend.justification = "center",
            legend.background = element_blank(),
            legend.box = "horizontal",
            legend.title = element_text(size=10),
            plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"),
            panel.border = element_blank(),
            plot.background = element_blank(),
            axis.text.x = element_text(face="bold"),
            axis.text.y = element_text(face="bold"),
            axis.title.x = element_text(face="bold"),
            axis.title.y = element_text(face="bold"),
            axis.line.x = element_line(size = 0.3),
            axis.line.y = element_line(size = 0.3))

burst_colors <- c("#fde725", "#35b779", "#31688e", "#440154")
names(burst_colors) <- c(1, 2, 3, 4)

#-------------------------------------------------------------------------------

pick_parameters <- function(
    bp_rule = "runif(1, min = 0.5, max = 3)",
    dp_rule = 0.5,
    rq_rule = 0.5,
    t_start_dist,
    t_run_rule,
    nr_of_families = NULL,
    nr_burst_divs = 3,
    response_nr = 1,
    prev_parameters = NULL,
    quality_dist = NULL,
    quality_noise = FALSE,
    q_noise_dist = NULL,
    uniform_fam = FALSE,
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
      cell_type <- replace(cell_type, which(cell_type == 1), "P")
      cell_type <- replace(cell_type, which(cell_type == 0), "Q")
    } else if (ASD == T){
      cell_type <- c("Q", rep("P", cells[i] - 1))
    }      
    if (quality_noise == T){
      if (is.null(q_noise_dist)){
        quality_fam <- round(rnorm(cells[i], mean = quality[i], sd = 0.05), digits = 2)
      } else if (!is.null(q_noise_dist)){
        quality_fam <- round(eval(parse(text = q_noise_dist)), digits = 2)
      }
      if (uniform_fam == T){
        warning("Can't combine noise on the quality with a uniform family")
      }
    } else {
      quality_fam <- rep(quality[i], cells[i])
    }
    if (uniform_fam == T){
      q <- quality[i]
      #This makes all proliferating branches in families uniform
      t_run_expr <- substitute(eval(parse(text = t_run_rule)))
      t_run_val <- eval(t_run_expr)
      t_run <- round(t_run_val, digits = 2)
      bp_fam <- round(eval(parse(text = bp_rule)), digits = 2)
      dp_fam <- round(eval(parse(text = dp_rule)), digits = 2)
    }
    
    for(cell in seq(cells[i])){
      # First define the expression and then evaluate the expression
      # If this is not included, the expression is only evaluated once,
      # leading to the same t_run for all (sub)families
      q <- quality_fam[cell]
      
      if (uniform_fam == F){
        t_run_expr <- substitute(eval(parse(text = t_run_rule)))
        t_run_val <- eval(t_run_expr)
        t_run <- round(t_run_val, digits = 2)
      }
      # this is to make sure all proliferation stops at max_run_time
      if (is.null(max_run_time)) {
        t_correct <- 0
      } else if (t_start[i] + t_run >= max_run_time) {
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
      
      if (cell_type[cell] == "P" & uniform_fam == F){
        bp <- round(eval(parse(text = bp_rule)), digits = 2)
        dp <- round(eval(parse(text = dp_rule)), digits = 2)
      } else if (cell_type[cell] == "Q"){
        bp <- 0
        dp <- 0
      } else if (uniform_fam == T){
        #Still need to assign bp in this loop, otherwise bp and dp will be 0 
        #after there has been a Q cell in the family
        bp <- bp_fam
        dp <- dp_fam
      }
      
      parameters %<>% add_row(
        cell_type = cell_type[cell],
        div_counter = div_counter[i],
        t_start = t_start[i],
        quality = quality_fam[cell],
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

# Optimized get_famsize function
get_famsize <- function(parameters, timepoint) {
  time_start_prolif <- parameters$t_start + parameters$t_burst * parameters$nr_burst_divs
  time_stop_prolif <- parameters$t_start + parameters$t_run - parameters$t_correct
  
  if (timepoint <= parameters$t_start) { 
    famsize <- 1 / 2^parameters$nr_burst_divs 
  } else if (timepoint <= time_start_prolif) { 
    div_times <- seq(parameters$t_start + parameters$t_burst, 
                     parameters$t_start + parameters$t_burst * parameters$nr_burst_divs, 
                     by = parameters$t_burst)
    index <- sum(timepoint >= div_times)
    famsize <- 2^index / 2^parameters$nr_burst_divs 
  } else if (timepoint <= time_stop_prolif) { 
    famsize <- exp(parameters$bp * (timepoint - time_start_prolif))
  } else { 
    max_cells <- exp(parameters$bp * (time_stop_prolif - time_start_prolif))
    famsize <- max_cells * exp(-parameters$dp * (timepoint - time_stop_prolif))
  }
  
  return(famsize)
}

#-------------------------------------------------------------------------------

# Optimized solve_function using apply
solve_function <- function(parameters, min_time, max_time, interval,
                           model = "classic", recruited_cells = NULL) {
  timepoints <- seq(min_time, max_time, by = interval)
  num_timepoints <- length(timepoints)
  
  if(!is.null(recruited_cells)){
    recruited_cells <- rowSums(recruited_cells[, -1])
  }
  
  # Create a data frame with pre-filled columns
  P_total <- data.frame(matrix(0, nrow = nrow(parameters), ncol = num_timepoints + 3))
  colnames(P_total) <- c("fam_nr", "fam_nr_2", "fam_nr_3", timepoints)
  P_total$fam_nr <- as.factor(parameters$fam_nr)
  P_total$fam_nr_2 <- as.factor(parameters$fam_nr_2)
  P_total$fam_nr_3 <- as.factor(parameters$fam_nr_3)
  
  # Precompute times for burst division phases
  div_times_list <- lapply(1:nrow(parameters), function(i) {
    seq(parameters$t_start[i] + parameters$t_burst[i], 
        parameters$t_start[i] + parameters$t_burst[i] * parameters$nr_burst_divs[i], 
        by = parameters$t_burst[i])
  })
  
  # Function to compute famsizes for a single cell across all timepoints
  compute_famsizes <- function(cell_index) {
    params <- parameters[cell_index, ]
    time_start_prolif <- params$t_start + params$t_burst * params$nr_burst_divs
    time_stop_prolif <- params$t_start + params$t_run - params$t_correct
    
    sapply(timepoints, function(timepoint) {
      if (timepoint <= params$t_start) {
        1 / 2^params$nr_burst_divs
      } else if (timepoint <= time_start_prolif) {
        index <- sum(timepoint >= div_times_list[[cell_index]])
        2^index / 2^params$nr_burst_divs
      } else if (timepoint <= time_stop_prolif) {
        exp(params$bp * (timepoint - time_start_prolif))
      } else {
        max_cells <- exp(params$bp * (time_stop_prolif - time_start_prolif))
        if(model == "classic"){
          max_cells * exp(-params$dp * (timepoint - time_stop_prolif))
        } else if(model == "noQ"){
          recruited_cells[cell_index] + 
            ((max_cells - recruited_cells[cell_index]) * exp(-params$dp * (timepoint - time_stop_prolif)))
        }
      }
    })
  }
  
  # Apply the compute_famsizes function to each cell
  famsize_matrix <- t(sapply(1:nrow(parameters), compute_famsizes))
  
  # Populate P_total with the computed famsizes
  P_total[, 4:(num_timepoints + 3)] <- famsize_matrix
  
  return(P_total)
}

#-------------------------------------------------------------------------------

get_max_famsize <- function(parameters){
  time_start_prolif <- parameters$t_start + parameters$t_burst*parameters$nr_burst_divs
  time_stop_prolif <- parameters$t_start + parameters$t_run - parameters$t_correct
  max_cells <- exp((parameters$bp) * (time_stop_prolif - time_start_prolif))
  return(max_cells)
}

#-------------------------------------------------------------------------------

generate_total_response_table <- function(parameters, min_time, max_time, interval, 
                                          model = "classic", recruited_cells = NULL){
  famsizes_table <- solve_function(parameters, min_time, max_time, interval, 
                                   model = model, recruited_cells = recruited_cells)
  timepoints <- seq(min_time, max_time, interval)
  num_timepoints <- length(timepoints)
  total_resp <- apply(famsizes_table[,4 : (num_timepoints + 3)], 2, FUN = sum)
  response <- data.frame(time = timepoints, cells = total_resp, log_cells = log10(total_resp))
  
  return(response)
}

#-------------------------------------------------------------------------------

get_total_response_stats <- function(response){
  max_cells <- max(response$cells)
  time_max <- response$time[response$cells == max_cells]
  # Assume that the duration of the table is long enough so that towards the end
  # there's only Q cells left
  Q_cells <- round(response$cells[response$time == max(response$time)])
  
  stats <- data.frame(max_cells = as.numeric(max_cells), 
                      time_max = as.numeric(time_max), 
                      Q_cells = as.numeric(Q_cells))
  
  return(stats)
}

#-------------------------------------------------------------------------------

plot_response <- function(response){
  plot <- ggplot(data = response, aes(x = time, y = log_cells)) +
    geom_point(size = 1) +
    labs(
      x = "Time (days)",
      y = "Nr. of cells (log-scale)"
    ) + theme_clean() + 
    th
  
  return(plot)
}

#-------------------------------------------------------------------------------

fill_missing_families <- function(Q_cells, nr_of_families = 500){
  missing_fams <- setdiff(seq(nr_of_families), Q_cells$fam_nr)
  remaining_values <- data.frame(fam_nr = as.numeric(missing_fams), 
                                 Q_cells = rep(0, length(missing_fams)))
  Q_cells <- rbind(Q_cells, remaining_values) 
  Q_cells <- Q_cells[order(Q_cells$fam_nr),]
  return(Q_cells)
}

#-------------------------------------------------------------------------------

generate_max_fam_table <- function(prim_parameters, 
                                   sec_parameters, 
                                   ter_parameters, 
                                   type = c("time", "max"), 
                                   timepoint = NULL, 
                                   remove_zeroes = FALSE,
                                   model = "classic",
                                   recruited_cells_prim = NULL,
                                   recruited_cells_sec = NULL,
                                   recruited_cells_ter = NULL) {
  
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
  
  if(model == "classic"){
    Q_cells <- prim_parameters %>% group_by(fam_nr) %>%
      summarize(Q_cells = sum(cell_type == "Q"))
    Q_cells_sec <- sec_parameters %>% group_by(fam_nr) %>%
      summarize(Q_cells = sum(cell_type == "Q"))
    Q_cells_ter <- ter_parameters %>% group_by(fam_nr) %>%
      summarize(Q_cells = sum(cell_type == "Q"))
  } else if(model == "noQ"){
    print("No Q model")
    if (is.null(recruited_cells_prim) | is.null(recruited_cells_sec) | is.null(recruited_cells_ter)){
      warning("Please provide the dataframes with recruited cells when working with the noQ model")}
    Q_cells <- get_recruited_cells_per_fam(recruited_cells_prim)
    Q_cells_sec <- get_recruited_cells_per_fam(recruited_cells_sec)
    Q_cells_ter <- get_recruited_cells_per_fam(recruited_cells_ter)
  }
  
  Q_cells <- fill_missing_families(Q_cells)
  Q_cells_sec <- fill_missing_families(Q_cells_sec)
  Q_cells_ter <- fill_missing_families(Q_cells_ter)
  
  max_cells <- data.frame(
    fam_nr = seq(1, nr_of_families),
    cells_prim = max_prim,
    cells_sec = max_sec,
    cells_ter = max_ter,
    nr_burst_divs = factor(nr_burst_divs),
    Q_cells = Q_cells$Q_cells,
    Q_cells_sec = Q_cells_sec$Q_cells,
    Q_cells_ter = Q_cells_ter$Q_cells
  )
  if (remove_zeroes == TRUE) {
    max_cells <- max_cells[max_cells$cells_ter > 0, ]
  }
  
  return(max_cells)
}

#-------------------------------------------------------------------------------

get_max_fam_stats <- function(max_cells, table = T){
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
  if(table == T){
    stats <- data.frame(mean = c(mean_prim, mean_sec, mean_ter),
                    median = c(median_prim, median_sec, median_ter),
                    correlation = c(cor_prim_sec$estimate, 0, cor_sec_ter$estimate),
                    p_value = c(cor_prim_sec$p.value, 0, cor_sec_ter$p.value))}
  else if(table == F){
    stats <- data.frame(mean_prim, mean_sec, mean_ter, 
                        median_prim, median_sec, median_ter,
                        cor_prim_sec = cor_prim_sec$estimate, p_value_prim_sec = cor_prim_sec$p.value,
                        cor_sec_ter = cor_sec_ter$estimate, p_value_sec_ter = cor_sec_ter$p.value)}

  return(stats)
}

#-------------------------------------------------------------------------------

plot_prim_sec_max_fam <- function(max_cells, 
                                  show_title = FALSE,
                                  label_burst_divs = c("col", "shape", "both")){ 
  
  # Create the geom layer based on the condition
  geom_layer <- 
    if(label_burst_divs == "col") {
      list(
        geom_point(aes(col = as.factor(nr_burst_divs), size = Q_cells), alpha = 0.75),
        scale_size(range = c(0.5, 2.5), breaks = seq(0, max(max_cells$Q_cells))),
        scale_color_manual(values = burst_colors),
        labs(color = "Nr. of burst divisions (prim. resp.)")) 
  } else if(label_burst_divs == "shape") {
      list(
        geom_point(aes(col = Q_cells, shape = as.factor(nr_burst_divs)), 
                   size = 1, alpha = 0.75),
        scale_color_viridis_c(option = "inferno", direction = -1,
                              limits = c(0, 15),
                              labels = c(0, 5, 10, 15),
                              breaks = c(0, 5, 10, 15),
                              alpha = 0.75), 
        guides(colour = guide_colourbar(show.limits = TRUE, 
                                 title.position = "top",
                                 barwidth = 10,
                                 barheight = 0.5)),
        labs(color = "Nr. of Q cells (prim. resp.)"))
  } else if(label_burst_divs == "both") {
      list(
        geom_point(aes(shape = as.factor(nr_burst_divs), col = as.factor(nr_burst_divs)), 
                   size = 1, alpha = 0.75),
        labs(color = "Nr. of burst divisions (prim. resp.)"),
        scale_color_manual(values = burst_colors)) 
  } else {
    stop("Invalid label") 
  }
  
  plot <- ggplot(data = max_cells, 
                          aes(x = log10(cells_prim), y = log10(cells_sec))) +
    geom_layer +
    labs(title = "Primary vs. Secondary Response",
         x = "Family size primary (log-scale)",
         y = "Family size secondary (log-scale)",
         shape = "Nr. of burst divisions (prim. resp.)") +
    theme_clean() +
    th + 
    theme(legend.box.background = element_rect(color="lightgrey")) +
    guides(size = "none", shape = guide_legend(title.position = "top"))  +
    scale_shape_manual(values = c(20, 17, 15))
  
  if (show_title == FALSE) {
    plot <- plot + labs(title = NULL, subtitle = NULL)
  }
  
  return(plot)
} 


#-------------------------------------------------------------------------------

plot_sec_ter_max_fam <- function(max_cells, 
                                 show_title = FALSE,
                                 label_burst_divs = c("col", "shape", "both")){ 
  
  # Create the geom layer based on the condition
  geom_layer <- 
    if(label_burst_divs == "col") {
      list(
        geom_point(aes(col = as.factor(nr_burst_divs), size = Q_cells), alpha = 0.75),
        scale_color_manual(values = burst_colors),
      scale_size(range = c(0.5, 2.5), breaks = seq(0, max(max_cells$Q_cells)))) 
    } else if(label_burst_divs == "shape") {
      list(
        geom_point(aes(col = Q_cells, shape = as.factor(nr_burst_divs)), size = 1, alpha = 0.75),
        scale_color_viridis_c(option = "inferno", direction = -1,
                              limits = c(0, 15),
                              labels = c(0, 5, 10, 15),
                              breaks = c(0, 5, 10, 15),
                              alpha = 0.75))
    } else if(label_burst_divs == "both") {
      list(
        geom_point(aes(shape = as.factor(nr_burst_divs), 
                       col = as.factor(nr_burst_divs)), 
                   size = 1, alpha = 0.75),
        scale_color_manual(values = burst_colors)) 
    } else {
      stop("Invalid label") 
    }
  
  plot <- ggplot(data = max_cells, 
                         aes(x = log10(cells_sec), y = log10(cells_ter))) +
    geom_layer +
    labs(
      title = "Secondary vs. Tertiary Response",
      x = "Family size secondary (log-scale)",
      y = "Family size tertiary (log-scale)") +
    theme_clean() +
    th +  
    guides(size = "none", color = "none", 
           shape = "none") +
    scale_shape_manual(values = c(20, 17, 15))

  
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
  df_famsizes$fraction <- df_famsizes$famsize / sum(df_famsizes$famsize)
  
  return(df_famsizes)
}

#-------------------------------------------------------------------------------

get_famsize_stats <- function(df_famsizes){
  stat_mean <- round(mean(df_famsizes$famsize[df_famsizes$famsize > 1]), digits = 0)
  stat_sd <- round(sd(df_famsizes$famsize[df_famsizes$famsize > 1]), digits = 3)
  stat_median <- round(median(df_famsizes$famsize[df_famsizes$famsize > 1]), digits = 0)
  top_5perc_families <- nrow(df_famsizes) * 0.05
  stat_disparity <- round(sum(df_famsizes$fraction[1:top_5perc_families]), digits = 3)
  
  return(c(stat_mean, stat_median, stat_sd, stat_disparity))
}

#-------------------------------------------------------------------------------

generate_freq_famsize_table <- function(logfamsizes){
  freq_famsizes <- data.frame(table(round(logfamsizes)))
  colnames(freq_famsizes) <- c("logfamsize", "freq")
  freq_famsizes$logfamsize <- as.numeric(as.character(freq_famsizes$logfamsize))
  freq_famsizes$freq <- freq_famsizes$freq/sum(freq_famsizes$freq)
  if(setequal(freq_famsizes$logfamsize, seq(max(freq_famsizes$logfamsize))) == F){
    missing_bins <- setdiff(c(0, seq(max(freq_famsizes$logfamsize))), 
                            freq_famsizes$logfamsize)
    remaining_values <- data.frame(logfamsize = as.numeric(missing_bins), 
                                   freq = rep(0, length(missing_bins)))
    freq_famsizes <- rbind(freq_famsizes, remaining_values) 
    freq_famsizes <- freq_famsizes[order(freq_famsizes$logfamsize),]
  }
  return(freq_famsizes)
}

#-------------------------------------------------------------------------------

plot_famsize_distribution <- function(freq_famsizes, timepoint, show_title = FALSE, 
                                      x_axis_max = NULL, y_axis_max = NULL, 
                                      multi_plot = FALSE){
  if (multi_plot == TRUE) {
    freq_famsizes$freq <- freq_famsizes$means
    freq_famsizes$logfamsize <- as.numeric(c(0, seq(nrow(freq_famsizes) - 1)))
  } else {
    freq_famsizes$lower_bound <- NA
    freq_famsizes$upper_bound <- NA
    freq_famsizes$sd <- NA
  }
  
  if (!is.data.frame(freq_famsizes)) {
    stop("freq_famsizes must be a data frame")
  }
  
  if (timepoint %in% c(5, 6, 7, 8)) {
    expected_famsizes <- get(paste0("day", timepoint))
  } else {
    expected_famsizes <- data.frame(Number_2log = 0, Frequency = 0)
  }
  
  expected_famsizes$Number_2log <- as.numeric(expected_famsizes$Number_2log)
  expected_famsizes$Frequency <- as.numeric(expected_famsizes$Frequency)
  
  freq_famsizes$logfamsize <- as.numeric(freq_famsizes$logfamsize)
  freq_famsizes$freq <- as.numeric(freq_famsizes$freq)
  freq_famsizes$lower_bound <- as.numeric(freq_famsizes$lower_bound)
  freq_famsizes$upper_bound <- as.numeric(freq_famsizes$upper_bound)
  
  if (is.null(x_axis_max)) {
    x_axis_max <- max(freq_famsizes$logfamsize, expected_famsizes$Number_2log, na.rm = TRUE) + 1
  }
  
  if (is.null(y_axis_max)) {
    y_axis_max <- max(freq_famsizes$freq, expected_famsizes$Frequency, na.rm = TRUE) + 0.05
  }
  
  if (max(freq_famsizes$logfamsize) < x_axis_max){
    nr_rows <- x_axis_max - max(freq_famsizes$logfamsize)
    bins <- seq(x_axis_max)
    addition <- data.frame(logfamsize = tail(bins, nr_rows), 
                           freq = rep(0, nr_rows),
                           lower_bound = rep(NA, nr_rows),
                           upper_bound = rep(NA, nr_rows),
                           sd = rep(NA, nr_rows))
    freq_famsizes <- rbind(freq_famsizes, addition)
  }
  
  plot <- ggplot(data = freq_famsizes, aes(x = logfamsize)) + 
    geom_bar(aes(y = freq), stat = "identity", fill = "#9f2a63", col = "black") +
    geom_errorbar(aes(ymin = freq - sd, ymax = freq + sd), width = 0.2) +
    geom_point(data = expected_famsizes, aes(x = Number_2log, y = Frequency), 
               color = "#9c9797", size = 1.5, shape = 19) + 
    labs(title = paste("Day", timepoint),
         x = "Family size (2log)",
         y = "Frequency") +
    theme_clean() +
    th +
    theme(plot.margin = margin(t = 0.7, r = 0, b = 0.5, l = 0, unit = "cm"),
          axis.text.x = element_text(face = NULL, size = 7.5)) +
    scale_x_continuous(breaks = 0:x_axis_max, guide = guide_axis(angle = 60)) +
    scale_y_continuous(labels = label_number(accuracy = 0.01),
                       expand = c(0, 0),
                       limits = c(0, y_axis_max)) 
  
  if (!show_title) {
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
    theme_clean() +
    th
  
  if (show_title == F) {
    plot <- plot + labs(title = NULL, subtitle = NULL)
  }
  return(plot)
}

#-------------------------------------------------------------------------------

generate_famsize_table_multidays <- function(parameters){
  famsizes_table <- lapply(5:8, function(day) {
    generate_famsize_table(parameters, timepoint = day)
  })
  
  return(famsizes_table)
} 

#-------------------------------------------------------------------------------

get_famsize_stats_multidays <- function(famsizes_table){
  famsizes_stats <- lapply(5:8, function(day) {
    c(day, get_famsize_stats(famsizes_table[[day-4]]))
  })
  famsizes_stats <- data.frame(do.call(rbind, famsizes_stats))
  colnames(famsizes_stats) <- c("timepoint", "mean", "median", "sd", "disparity")
  
  return(famsizes_stats)
}

#-------------------------------------------------------------------------------

generate_freq_famsize_table_multidays <- function(famsizes_table){
  frequencies_table <- lapply(1:4, function(day) {
    generate_freq_famsize_table(famsizes_table[[day]]$logfamsize)
  })
  
  return(frequencies_table)
}

#-------------------------------------------------------------------------------

plot_grid_famsize_dist <- function(famsizes_table = NULL, 
                                     frequencies_table = NULL, 
                                     x_axis_max = NULL, 
                                     y_axis_max = NULL, 
                                     multi_plot = FALSE){
  if (multi_plot == TRUE) {
    for (i in seq(4)) {
      frequencies_table[[i]]$freq <- frequencies_table[[i]]$means
      frequencies_table[[i]]$logfamsize <- seq(0, nrow(frequencies_table[[i]]) - 1)
    }
  }
  
  if (is.null(frequencies_table) && !is.null(famsizes_table)) {
    frequencies_table <- generate_freq_famsize_table_multidays(famsizes_table)
  }
  
  if (is.null(frequencies_table)) {
    stop("frequencies_table is NULL. Check the generate_freq_famsize_table_multidays function.")
  }
  
  if (!all(sapply(frequencies_table, is.data.frame))) {
    stop("frequencies_table must be a list of data frames")
  }
  
  # Determine x_axis_max for all plots
  if (is.null(x_axis_max)) {
    x_axis_max <- max(sapply(frequencies_table, function(df) max(as.numeric(df$logfamsize), na.rm = TRUE))) + 1
    x_axis_max <- rep(x_axis_max, 4)  # Replicate for each plot
  }
  
  # Determine y_axis_max for all plots
  if (is.null(y_axis_max)) {
    max_freqs <- sapply(frequencies_table, function(df) max(as.numeric(df$freq), na.rm = TRUE))
    y_axis_max <- c(max_freqs[1], max_freqs[2], max_freqs[3], max_freqs[4])  + 0.05 # Adjust for the last plot
  }
  
  plots <- lapply(1:4, function(day) {
    plot_famsize_distribution(frequencies_table[[day]], 
                              day + 4, 
                              x_axis_max = x_axis_max[day], 
                              y_axis_max = y_axis_max[day],
                              multi_plot = multi_plot)
  })
  
  plots_no_xlabel <- lapply(1:3, function(day) {
    plots[[day]] + theme(axis.title.x = element_blank(),
                         plot.margin = margin(t = 0.4, r = 0.5, b = 0, l = 0.5, unit = "cm"))
  })
  
  plots_one_xlabel <- list(plots_no_xlabel[[1]], 
                           plots_no_xlabel[[2]], 
                           plots_no_xlabel[[3]], 
                           plots[[4]] + theme(plot.margin = margin(t = 0.4, r = 0.5, b = 0.2, l = 0.5, unit = "cm")))
  
  plots[[4]] <- plots[[4]] + theme(plot.margin = margin(t = 0.4, r = 0.5, b = 0.2, l = 0.5, unit = "cm"))
  
  gt <- arrangeGrob(grobs = list(plots_no_xlabel[[1]], 
                                 plots_no_xlabel[[2]], 
                                 plots_no_xlabel[[3]], 
                                 plots[[4]]), nrow = 4, heights = c(4,4,4,4.6))
  
  plot <- as_ggplot(gt) + 
    draw_plot_label(label = c("A - d5", 
                              "B - d6", 
                              "C - d7", 
                              "D - d8"), size = 13,
                    y = c(1.002, 0.76, 0.53, 0.28), x = rep(-0.03, 4))
  
  return(plot)
}
# _______________________________________________________________________________
# This function plots the size of the family as a function of the number of Q cells that are formed

generate_Q_famsize_table <- function(parameters) {
  max_famsizes <- get_max_famsize(parameters)
  Q_famsize <- as_tibble(data.frame(fam_nr = as.factor(parameters$fam_nr), 
                                cell_type = parameters$cell_type, 
                                nr_burst_divs = parameters$nr_burst_divs, 
                                max_famsize = max_famsizes))
  
  Q_famsize_table <- Q_famsize %>%
    group_by(fam_nr) %>%
    summarise(
      famsize = sum(max_famsize),
      Q_cells = sum(cell_type == "Q"),
      nr_burst_divs = first(nr_burst_divs)  # Since nr_burst_divs is the same for each fam_nr
    )
  
  Q_famsize_table$log_famsize <- log10(Q_famsize_table$famsize)
  
  return(Q_famsize_table)
}

#-------------------------------------------------------------------------------

get_recruited_cells_per_fam <- function(recruited_cells){
  recruited_cells_per_fam <- recruited_cells %>%
    group_by(fam_nr) %>%
    summarise(Q_cells = n())
  return(recruited_cells_per_fam)
}

#-------------------------------------------------------------------------------

generate_recruited_famsize_table <- function(parameters, recruited_cells){
  recruited_cells_per_fam <- get_recruited_cells_per_fam(recruited_cells)
  
  max_famsizes <- get_max_famsize(parameters)
  famsize_df <- as_tibble(data.frame(fam_nr = as.factor(parameters$fam_nr), 
                                     nr_burst_divs = parameters$nr_burst_divs, 
                                     max_famsize = max_famsizes))
  famsize_df <- famsize_df %>%
    group_by(fam_nr) %>%
    summarise(famsize = sum(max_famsize),
              nr_burst_divs = first(nr_burst_divs))
  
  famsize_df <- merge(famsize_df, recruited_cells_per_fam, by = "fam_nr", all.x = T)
  famsize_df[is.na(famsize_df)] <- 0
  
  famsize_df$log_famsize <- log10(famsize_df$famsize)
  
  return(famsize_df)
}

#-------------------------------------------------------------------------------

get_Q_famsize_stats <- function(max_fam_table){
  prim_cor <- cor.test(max_fam_table$cells_prim, max_fam_table$Q_cells, method = "spearman")
  sec_cor <- cor.test(max_fam_table$cells_sec, max_fam_table$Q_cells, method = "spearman")
  ter_cor <- cor.test(max_fam_table$cells_ter, max_fam_table$Q_cells, method = "spearman")
  
  correlations <- data.frame(r = c(prim_cor$estimate, sec_cor$estimate, ter_cor$estimate),
                             p_value = c(prim_cor$p.value, sec_cor$p.value, ter_cor$p.value))
  
  return(correlations)
}

#-------------------------------------------------------------------------------

get_Q_prim_famsize_stats <- function(Q_famsize_table){
  prim_cor <- cor.test(Q_famsize_table$famsize, Q_famsize_table$Q_cells, method = "spearman")
  
  return(prim_cor)
}

#-------------------------------------------------------------------------------

plot_Q_famsize <- function(Q_famsize_table, 
                           show_legend = F, 
                           label_burst_divs = c("shape", "col", "both"),
                           linear_model = F){
  # Create the geom layer based on the condition
  geom_layer <- if(label_burst_divs == "col") {
    list(
      geom_point(aes(col = as.factor(nr_burst_divs)), size = 0.6, alpha = 0.8),
      scale_color_manual(values = burst_colors),
      guides(col = guide_legend(title.position = "top", ncol = 1)),
      theme(legend.position = "right",
            #legend.justification = c("left", "top"),
            legend.box.background = element_rect(color="lightgrey")))
  } else if(label_burst_divs == "shape") {
    geom_point(aes(shape = as.factor(nr_burst_divs)), size = 0.6, alpha = 0.8)
  } else if(label_burst_divs == "both") {
    list(
      geom_point(aes(shape = as.factor(nr_burst_divs), 
                     col = as.factor(nr_burst_divs)), 
                     size = 0.8, alpha = 0.8),
      scale_color_manual(values = burst_colors),
      guides(col = guide_legend(title.position = "top", ncol = 1)),
      theme(legend.position = "right",
            #legend.justification = c("left", "top"),
            legend.box.background = element_rect(color="lightgrey")))
  } else {
    stop("Invalid label") 
  }
  lm_layer <- if(linear_model == T){
    list(
      geom_smooth(method = lm, aes(col = as.factor(nr_burst_divs))),
      geom_smooth(method = lm, col = "black"))
  }
  
  plot <- ggplot(data = Q_famsize_table, (aes(x = log_famsize, y = Q_cells))) +
    theme_clean() + th +
    geom_layer +
    labs(
      x = "Max. fam. size prim. response (log-scale)",
      y = "Nr. of Q cells",
      shape = "Nr. of burst divisions",
      col = "# b. div."
    ) + 
    scale_shape_manual(values = c(20, 17, 15)) +
    lm_layer
  
  if (show_legend == F) {
    plot <- plot + theme(legend.position = "none")
  } 
  
  return(plot)
}

#-------------------------------------------------------------------------------

generate_Q_mean_famsize_table <- function(max_fam_table){
  mean_max_fam_table <- max_fam_table %>% 
    group_by(Q_cells) %>% 
    summarise(mean_prim = mean(cells_prim), 
              mean_sec = mean(cells_sec), 
              mean_ter = mean(cells_ter),
              sd_prim = sd(cells_prim),
              sd_sec = sd(cells_sec),
              sd_ter = sd(cells_ter),
              n = n())
  
  return(mean_max_fam_table)
}

#-------------------------------------------------------------------------------

plot_Q_famsize_boxplots <- function(max_fam_table){
  mean_max_fam_table <- generate_Q_mean_famsize_table(max_fam_table)
  myxlab <- paste(as.factor(mean_max_fam_table$Q_cells), "\nN=", mean_max_fam_table$n, sep="")
  
  plot_Q_prim <- ggplot(max_fam_table, aes(x = as.factor(Q_cells), fill = Q_cells)) + 
    geom_boxplot(aes(y = log10(cells_prim))) +
    scale_fill_viridis_c(option = "inferno", direction = -1,
                         limits = c(0, 15),
                         labels = c(0, 5, 10, 15),
                         breaks = c(0, 5, 10, 15)) +
    guides(fill = guide_colourbar(show.limits = TRUE, 
                                  title.position = "top",
                                  barwidth = 10,
                                  barheight = 0.5)) +
    labs(
      x = "Q cells in prim. resp.",
      y = "log10(family size)",
      fill = "Nr. of Q cells (prim. resp.)") + 
    theme_clean() + th +
    theme(legend.box.background = element_rect(color="lightgrey"),
          legend.position = "none",
          axis.text.x = element_blank(), 
          axis.title.x = element_blank()) +
    guides(size = "none", shape = guide_legend(title.position = "top")) 
  
  plot_Q_sec <- ggplot(max_fam_table, aes(x = as.factor(Q_cells), fill = Q_cells)) + 
    geom_boxplot(aes(y = log10(cells_sec))) +
    scale_fill_viridis_c(option = "inferno", direction = -1,
                         limits = c(0, 15),
                         labels = c(0, 5, 10, 15),
                         breaks = c(0, 5, 10, 15)) +
    guides(fill = guide_colourbar(show.limits = TRUE, 
                                  title.position = "top",
                                  barwidth = 10,
                                  barheight = 0.5)) +
    labs(
      x = "Q cells in prim. resp.",
      y = "log10(family size)",
      fill = "Nr. of Q cells (prim. resp.)") + 
    theme_clean() + th +
    theme(legend.box.background = element_rect(color="lightgrey"),
          legend.position = "none",
          axis.text.x = element_blank(), 
          axis.title.x = element_blank()) +
    guides(size = "none", shape = guide_legend(title.position = "top")) 
  
  plot_Q_ter <- ggplot(max_fam_table, aes(x = as.factor(Q_cells), fill = Q_cells)) + 
    geom_boxplot(aes(y = log10(cells_ter))) +
    scale_fill_viridis_c(option = "inferno", direction = -1,
                         limits = c(0, 15),
                         labels = c(0, 5, 10, 15),
                         breaks = c(0, 5, 10, 15)) +
    guides(fill = guide_colourbar(show.limits = TRUE, 
                                  title.position = "top",
                                  barwidth = 10,
                                  barheight = 0.5)) +
    labs(
      x = "Q cells in prim. resp.",
      y = "log10(family size)",
      fill = "Nr. of Q cells (prim. resp.)") + 
    theme_clean() + th +
    theme(legend.box.background = element_rect(color="lightgrey")) +
    guides(size = "none", shape = guide_legend(title.position = "top"))  +
    scale_x_discrete(labels=myxlab)
  
  grob_plots <- arrangeGrob(grobs = list(plot_Q_prim, plot_Q_sec, plot_Q_ter), 
                            ncol = 1, heights = c(0.285, 0.285, 0.43))
  
  plots <- as_ggplot(grob_plots) +
    draw_plot_label(label = c("A", "B", "C"),
                    x = c(0,0,0), 
                    y = c(1, 0.715, 0.43))
  
  return(plots)
}

#-------------------------------------------------------------------------------

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#-------------------------------------------------------------------------------

get_largest_fam <- function(max_fam_table){
  largest_fam <- max(max_fam_table$cells_prim)
  median_fam <- median(max_fam_table$cells_prim)
  return(c("largest_fam" = largest_fam, "median_fam" = median_fam))
}

#-------------------------------------------------------------------------------

get_famsize_fraction <- function(max_fam_table){
  famsizes <- max_fam_table$cells_prim
  famsizes <- famsizes[order(famsizes, decreasing = T)]
  total_famsize <- sum(famsizes)
  famsizes_perc <- (famsizes/total_famsize) * 100 
  return(famsizes_perc)
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Functions below are for model 2 with quality
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Function to skew a uniform distribution towards higher values based on quality
skewed_uniform <- function(n, min_val, max_val, quality) {
  # Beta distribution parameters
  alpha <- get_alpha(quality)
  beta <- get_beta(quality)
  
  # Sample from beta distribution and transform to uniform distribution
  beta_sample <- rbeta(n, alpha, beta)
  return(min_val + (max_val - min_val) * beta_sample)
}

#-------------------------------------------------------------------------------

generate_quality_effect_panelplot <- function(min_q, med_q, max_q){
  steps <- seq(0, 1, 0.02)
  
  low_q <- data.frame("value" = steps, 
                      "chance" = dbeta(steps, 
                                       shape1 = get_alpha(min_q), 
                                       shape2 = get_beta(min_q)))
  plot_quality_low <- ggplot(low_q, aes(x = value, y = chance)) + 
    geom_point() + theme_clean() + 
    theme(axis.text.y = element_blank(), axis.text.x = element_blank()) +
    labs(title = paste0("Minimum quality (", min_q, ")"))
  
  mid_q <- data.frame("value" = steps, 
                      "chance" = dbeta(steps, 
                                       shape1 = get_alpha(med_q), 
                                       shape2 = get_beta(med_q)))
  plot_quality_mid <- ggplot(mid_q, aes(x = value, y = chance)) + 
    geom_point() + theme_clean() + 
    theme(axis.text.y = element_blank(), axis.text.x = element_blank()) +
    labs(title = paste0("Medium quality (", med_q, ")"))
  
  high_q <- data.frame("value" = steps, 
                       "chance" = dbeta(steps, 
                                        shape1 = get_alpha(max_q), 
                                        shape2 = get_beta(max_q)))
  plot_quality_high <- ggplot(high_q, aes(x = value, y = chance)) + 
    geom_point() + theme_clean() + 
    theme(axis.text.y = element_blank(), axis.text.x = element_blank()) +
    labs(title = paste0("Maximum quality (", max_q, ")"))
  
  gt_quality <- arrangeGrob(grobs = list(plot_quality_low, 
                                         plot_quality_mid,
                                         plot_quality_high), 
                            ncol =3)
  
  panel_plots_quality <- as_ggplot(gt_quality)
  
  return(panel_plots_quality)
}


