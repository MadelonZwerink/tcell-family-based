#!/usr/bin/env Rscript

generate_recruited_cells <- function(div_cells, 
                                     recruitment_mean = 5,
                                     recruitment_sd = 2){
  div_cells_rounded <- div_cells
  
  for (r in 1:nrow(div_cells_rounded)){
    for (c in 1:n){
      numb <- runif(1)
      dec <- div_cells_rounded[r, c+1] %% 1
      if(numb <= dec){
        div_cells_rounded[r, c+1] <- 1 + as.integer(div_cells_rounded[r, c+1])
      }
      if(numb >= dec){
        div_cells_rounded[r, c+1] <- 0 + as.integer(div_cells_rounded[r, c+1])
      }
    }
  }
  
  #-------------------------------------------------------------------------------
  # Now, we can select each cell with a certain probability depending on 
  # division index
  recruitment_prob <- pnorm(0:n, mean = recruitment_mean, 
                            sd = recruitment_sd, lower.tail = FALSE)
  
  #-------------------------------------------------------------------------------
  # Check the number of Q cells and the division index of Q cells for different 
  # recruitment probability distributions
  
  # Select cells depending on probability
  selected_Q_cells <- data.frame(matrix(0, nrow = nrow(div_cells_rounded), 
                                        ncol = ncol(div_cells_rounded)))
  colnames(selected_Q_cells) <- colnames(div_cells_rounded)
  selected_Q_cells$fam_nr <- div_cells_rounded$fam_nr
  
  for (r in 1:nrow(selected_Q_cells)){
    for (c in 1:n){
      numb <- runif(div_cells_rounded[r, c+1])
      Q_cells <- sum(as.numeric(numb <= recruitment_prob[c]))
      selected_Q_cells[r, c+1] <- Q_cells
    }
  }
  
  return(selected_Q_cells)
}

#-------------------------------------------------------------------------------

generate_processed_recruited_cells <- function(recruited_cells){
  
  recruited_cells <- wide_to_long(recruited_cells, "Q_cells")
  
  recruited_cells <- recruited_cells[recruited_cells$Q_cells != 0,]
  
  processed_recruited_cells <- recruited_cells %>%
    uncount(Q_cells)
  processed_recruited_cells$cell_type <- "Q"
  colnames(processed_recruited_cells) <- c("fam_nr", "div_counter", "cell_type")
  
  return(processed_recruited_cells)
}

#-------------------------------------------------------------------------------

get_exp_divs <- function(parameters){
  exp_divs <- parameters$bp * (parameters$t_run - (parameters$t_burst * parameters$nr_burst_divs))
  return(exp_divs)
}

generate_div_index_df <- function(parameters, n = 100){
  exp_divs <- get_exp_divs(parameters)
  
  div_index <- data.frame(matrix(nrow = nrow(parameters), ncol = n))
  colnames(div_index) <- seq(0,n-1)
  
  for (fam in seq(nrow(parameters))) {
    dist_divs <- dpois(x = seq(0, n), lambda = exp_divs[fam])
    dist_divs_bursts <- c(rep(0, parameters$div_counter[fam]), dist_divs)
    div_index[fam,] <- dist_divs_bursts[0:n]
  }
  
  div_index <- cbind(fam_nr = parameters$fam_nr, div_index)
  
  return(div_index)
}

generate_div_nr_cells_df <- function(parameters, n = 100){
  exp_divs <- get_exp_divs(parameters)
  
  div_cells <- data.frame(matrix(nrow = nrow(parameters), ncol = n))
  colnames(div_cells) <- seq(0,n-1)
  
  for (fam in seq(nrow(parameters))) {
    dist_divs <- dpois(x = seq(0, n), lambda = exp_divs[fam])
    nr_cells <- dist_divs * 2^seq(0, n)
    
    dist_cells_bursts <- c(rep(0, parameters$div_counter[fam]), nr_cells)
    div_cells[fam,] <- dist_cells_bursts[0:n]
  }
  
  div_cells <- cbind(fam_nr = parameters$fam_nr, div_cells)
  
  return(div_cells)
}

wide_to_long <- function(div_df, names_values) {
  div_df %<>% pivot_longer(cols = 2:ncol(div_df), 
                           names_to = "div_index", 
                           values_to = eval(names_values))
  div_df$div_index <- as.numeric(div_df$div_index)
  return(div_df)
}
