families <- 500
famsizes_day <- list()

for (i in seq(100)) {
  prim_parameters <- pick_parameters(bp_rule = 'runif(1, min = 0, max = 3.5)',
                                     dp_rule = 0.4,
                                     rq_rule = 0.5,
                                     nr_of_families = families,
                                     response_nr = 1,
                                     t_run_rule = 'runif(1, min = 0.5, max = 4)',
                                     t_start_dist = 'rlnorm(nr_of_families, meanlog = 1.4, sdlog = 0.3)',
                                     nr_burst_divs = 'sample(c(2,3,4), nr_of_families, replace = TRUE)',
                                     min_t_start = 2.5)
  
  famsizes_day[[i]] <- list()
  plots <- plot_grid_famsize_dist(prim_parameters)
  
  for (day in 5:8) {
    famsizes_day[[i]][[paste0("day", day)]] <- plots[[3]][day - 4]  # Adjusted index
  }
  
  print(i)
}

mean_d <- sapply(5:8, function(day) mean(2^unlist(sapply(famsizes_day, function(x) x[[paste0("day", day)]]))[which(unlist(sapply(famsizes_day, function(x) x[[paste0("day", day)]])) > 1)]))
median_d <- sapply(5:8, function(day) median(2^unlist(sapply(famsizes_day, function(x) x[[paste0("day", day)]]))[which(unlist(sapply(famsizes_day, function(x) x[[paste0("day", day)]])) > 1)]))
freq_famsizes_d <- lapply(5:8, function(day) {
  freq <- data.frame(table(round(unlist(sapply(famsizes_day, function(x) x[[paste0("day", day)]])))))
  colnames(freq) <- c("logfamsize", "freq")
  freq$freq <- freq$freq/sum(freq$freq)
  return(freq)
})

plot_distribution_d <- lapply(5:8, function(day) {
  ggplot(data = freq_famsizes_d[[day-4]], aes(x = logfamsize)) + 
    geom_col(aes(y = freq), color = "darkblue", fill = "lightblue") +
    geom_point(data = get(paste0("day", day)), aes(x = Number_2log, y = Frequency), color = "darkorange", size = 2) +  # Adjusted data frame
    labs(title = paste("Family sizes day", day),
         x = "Family size (2log)",
         y = "Frequency") +
    theme(plot.margin = margin(t = 1, r = 1, 0, 0, unit = "cm")) +
    theme_clean() 
})

plot_grid(plotlist = plot_distribution_d, align = "hv", rel_widths = c(0.5, 0.5), rel_heights = c(0.5, 0.5))
mean_d
# mean_d 63.0352  513.1990 2342.4626 3645.9246
median_d
# median_d 20.29831  49.64314 126.09589 214.80143


