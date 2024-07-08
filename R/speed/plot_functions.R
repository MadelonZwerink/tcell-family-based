

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
