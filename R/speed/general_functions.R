
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

fill_missing_families <- function(Q_cells){
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

generate_famsize_table_multidays <- function(parameters){
  famsizes_table <- lapply(5:8, function(day) {
    generate_famsize_table(parameters, timepoint = day)
  })
  
  return(famsizes_table)
} 

#-------------------------------------------------------------------------------

generate_freq_famsize_table_multidays <- function(famsizes_table){
  frequencies_table <- lapply(1:4, function(day) {
    generate_freq_famsize_table(famsizes_table[[day]]$logfamsize)
  })
  
  return(frequencies_table)
}

#-------------------------------------------------------------------------------
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

#Returns the fraction that each family contributed to the total primary response
#ordered from large to small
get_famsize_fraction <- function(max_fam_table){
  famsizes <- max_fam_table$cells_prim
  famsizes <- famsizes[order(famsizes, decreasing = T)]
  total_famsize <- sum(famsizes)
  famsizes_perc <- (famsizes/total_famsize) * 100 
  return(famsizes_perc)
}
