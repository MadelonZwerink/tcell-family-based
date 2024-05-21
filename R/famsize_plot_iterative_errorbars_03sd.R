library(ggplot2)
library(gridExtra)

set.seed(12345)
families <- 500
famsizes_day <- list()
mean_values <- list()
median_values <- list()

for (i in seq(1000)) {
  prim_parameters <- pick_parameters(bp_rule = 'runif(1, min = 0, max = 3.5)',
                                     dp_rule = 0.4,
                                     rq_rule = 0.5,
                                     nr_of_families = families,
                                     response_nr = 1,
                                     t_run_rule = 'runif(1, min = 0.5, max = 4)',
                                     t_start_dist = 'rlnorm(nr_of_families, meanlog = 1.4, sdlog = 0.3)',
                                     nr_burst_divs = 'sample(c(2,3,4), nr_of_families, replace = TRUE)')
  
  famsizes_day[[i]] <- list()
  plots <- plot_grid_famsize_dist(prim_parameters, make_plot = F)
  
  for (day in 5:8) {
    famsizes_day[[i]][[paste0("day", day)]] <- plots[[3]][day - 4]  # Adjusted index
  }
  print(i)
}

calculate_frequencies <- function(inner_list) {
  freq <- table(factor(round(unlist(inner_list)), levels = seq(0, 20)))
  freq <- freq / sum(freq)
  return(as.data.frame(freq))
}

freq_list <- lapply(famsizes_day, function(element) {
  lapply(element, calculate_frequencies)
})

# Assuming freq_list is the list of lists containing frequencies for each run and each day

# Convert frequencies to numeric
freq_list_numeric <- lapply(freq_list, function(run) {
  lapply(run, function(day) {
    day[, "Freq"] <- as.numeric(as.character(day[, "Freq"]))
    return(day)
  })
})

# Calculate sum and standard deviation for each day across all runs
sum_freq <- lapply(1:4, function(day) {
  sum_day <- colSums(do.call(rbind, lapply(freq_list_numeric, function(run) run[[day]][, "Freq"])))
  values_day <- colMeans(do.call(rbind, lapply(freq_list_numeric, function(run) as.numeric(row.names(run[[day]]))))) # Take mean to keep the values unchanged
  return(data.frame(Sum = sum_day, Value = values_day))
})

sd_freq <- lapply(1:4, function(day) {
  sd_day <- apply(do.call(rbind, lapply(freq_list_numeric, function(run) run[[day]][, "Freq"])), 2, sd)
  return(sd_day)
})

# Calculate average frequency per day across all runs
average_freq <- lapply(1:4, function(day) {
  average_day <- colMeans(do.call(rbind, lapply(freq_list_numeric, function(run) run[[day]][, "Freq"])))
  return(average_day)
})

day_data_list <- list()

for (day in 1:4) {
  day_data_list[[day]] <- data.frame(Sum = sum_freq[[day]], SD = sd_freq[[day]], Average = average_freq[[day]])
  colnames(day_data_list[[day]]) <- c("Sum", "Value", "SD", "Average")
}

# Function to plot histogram with error bars
plot_histogram <- function(day_data_list, timepoint) {
  expected_famsizes <<- get(paste0("day", timepoint))
  
  ggplot(day_data_list, aes(x = Value, y = Average, ymin = Average - 2*SD, ymax = Average + 2*SD)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    geom_errorbar(width = 0.2) +
    geom_point(data = expected_famsizes, aes(x = Number_2log, y = Frequency), color = "darkorange", size = 2, inherit.aes = F) + 
    labs(x = "Family size (2log)", y = "Fraction", title = paste("Day", timepoint)) +
    theme_clean()
}

# Plot histograms for each day
histogram_plots <- lapply(1:4, function(day) {
  plot_histogram(day_data_list[[day]], timepoint = day + 4)
})

# Arrange and display the histograms
plot_famsizedist <- grid.arrange(grobs = histogram_plots, ncol = 2)

ggsave("./results/model1/plot_famsize_dist_95confint_1000average_03sd.png", plot = plot_famsizedist, 
       width = 2500, height = 1800, units = "px")
ggsave("./results/model1/plot_famsize_dist_95confint_1000average_03sd.pdf", plot = plot_famsizedist, 
       width = 2500, height = 1800, units = "px")

plot_famsizedist_vert <- grid.arrange(grobs = histogram_plots, ncol = 1)

ggsave("./results/model1/vert_plot_famsize_dist_95confint_1000average_03sd.png", plot = plot_famsizedist_vert, 
       width = 2500, height = 1800, units = "px")
ggsave("./results/model1/vert_plot_famsize_dist_95confint_1000average_03sd.pdf", plot = plot_famsizedist_vert, 
       width = 2500, height = 1800, units = "px")

# Function to calculate mean, median, and confidence interval
calculate_stats <- function(data) {
  means <- data.frame(day5 = numeric(), day6 = numeric(), day7 = numeric(), day8 = numeric())
  medians <- data.frame(day5 = numeric(), day6 = numeric(), day7 = numeric(), day8 = numeric())
  for (run in data) {
    day_values <- c()
    means_run <- c()
    medians_run <- c()
    for (day in run) {
      day_values <- unlist(day)
      day_values <- 2^day_values
      day_values <- day_values[day_values > 1]
      means_run <- c(means_run, mean(day_values))
      medians_run <- c(medians_run, median(day_values))
    }
    means %<>% add_row(day5 = means_run[1], day6 = means_run[2], day7 = means_run[3], day8 = means_run[4])
    medians %<>% add_row(day5 = medians_run[1], day6 = medians_run[2], day7 = medians_run[3], day8 = medians_run[4])
  }
  return(list(means, medians))
}

# Calculate stats
stats <- calculate_stats(famsizes_day)

means <- stats[[1]]

means_confint <- sapply(means, function(x) t.test(x)$conf.int)
means_values <- colMeans(means)

stats_means <- rbind(means_values, means_confint)
rownames(stats_means) <- c("Mean", "Lowerbound", "Upperbound")

medians <- stats[[2]]

medians_confint <- sapply(medians, function(x) t.test(x)$conf.int)
medians_values <- colMeans(medians)

stats_medians <- rbind(medians_values, medians_confint)
rownames(stats_medians) <- c("Median", "Lowerbound", "Upperbound")

sink(file = "./results/model1/stats_iterative_errorbars_1000average_03sd.txt")
stats_means

stats_medians
sink(file = NULL)

sink(file = "./results/model1/rawdata_iterative_errorbars_1000average_03sd.txt")

famsizes_day

sink(file = NULL)
