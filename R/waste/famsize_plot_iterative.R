families <- 500
famsizes_day5 <- list()
famsizes_day6 <- list()
famsizes_day7 <- list()
famsizes_day8 <- list()

for(i in seq(100)){
  prim_parameters <- pick_parameters(bp_rule = 'runif(1, min = 0.5, max = 3.5)',
                                   dp_rule = 0.4,
                                   rq_rule = 0.5,
                                   nr_of_families = families,
                                   response_nr = 1,
                                   t_run_rule = 'runif(1, min = 0.5, max = 4)',
                                   t_start_dist = 'rlnorm(nr_of_families, meanlog = 1.4, sdlog = 0.2)',
                                   nr_burst_divs = 'sample(c(1,2,3,4), nr_of_families, replace = TRUE)')
  plots <- plot_grid_famsize_dist(prim_parameters)
  famsizes_day5[i] <- plots[[3]][1]
  famsizes_day6[i] <- plots[[3]][2]
  famsizes_day7[i] <- plots[[3]][3]
  famsizes_day8[i] <- plots[[3]][4]
  print(i)
}
mean_d5 <- mean(2^unlist(famsizes_day5)[which(unlist(famsizes_day5) > 1)])
median_d5 <- median(2^unlist(famsizes_day5)[which(unlist(famsizes_day5) > 1)])
freq_famsizes_d5 <- data.frame(table(round(unlist(famsizes_day5))))
colnames(freq_famsizes_d5) <- c("logfamsize", "freq")
freq_famsizes_d5$freq <- freq_famsizes_d5$freq/sum(freq_famsizes_d5$freq)

mean_d6 <- mean(2^unlist(famsizes_day6)[which(unlist(famsizes_day6) > 1)])
median_d6 <- median(2^unlist(famsizes_day6)[which(unlist(famsizes_day6) > 1)])
freq_famsizes_d6 <- data.frame(table(round(unlist(famsizes_day6))))
colnames(freq_famsizes_d6) <- c("logfamsize", "freq")
freq_famsizes_d6$freq <- freq_famsizes_d6$freq/sum(freq_famsizes_d6$freq)

mean_d7 <- mean(2^unlist(famsizes_day7)[which(unlist(famsizes_day7) > 1)])
median_d7 <- median(2^unlist(famsizes_day7)[which(unlist(famsizes_day7) > 1)])
freq_famsizes_d7 <- data.frame(table(round(unlist(famsizes_day7))))
colnames(freq_famsizes_d7) <- c("logfamsize", "freq")
freq_famsizes_d7$freq <- freq_famsizes_d7$freq/sum(freq_famsizes_d7$freq)

mean_d8 <- mean(2^unlist(famsizes_day8)[which(unlist(famsizes_day8) > 1)])
median_d8 <- median(2^unlist(famsizes_day8)[which(unlist(famsizes_day8) > 1)])
freq_famsizes_d8 <- data.frame(table(round(unlist(famsizes_day8))))
colnames(freq_famsizes_d8) <- c("logfamsize", "freq")
freq_famsizes_d8$freq <- freq_famsizes_d8$freq/sum(freq_famsizes_d8$freq)

plot_distribution_d5 <- ggplot(data = freq_famsizes_d5, aes(x = logfamsize)) + 
  geom_col(aes(y = freq), color = "darkblue", fill = "lightblue") +
  geom_point(data = day5, aes(x = Number_2log, y = Frequency), color = "darkorange", size = 2) +
  labs(title = "Family sizes day 5",
       x = "Family size (2log)",
       y = "Frequency") +
  theme(plot.margin = margin(t = 1, r = 1, 0, 0, unit = "cm")) +
  theme_clean() 
plot_distribution_d5

plot_distribution_d6 <- ggplot(data = freq_famsizes_d6, aes(x = logfamsize)) + 
  geom_col(aes(y = freq), color = "darkblue", fill = "lightblue") +
  geom_point(data = day6, aes(x = Number_2log, y = Frequency), color = "darkorange", size = 2) +
  labs(title = "Family sizes day 6",
       x = "Family size (2log)",
       y = "Frequency") +
  theme(plot.margin = margin(t = 1, r = 1, 0, 0, unit = "cm")) +
  theme_clean() 
plot_distribution_d6

plot_distribution_d7 <- ggplot(data = freq_famsizes_d7, aes(x = logfamsize)) + 
  geom_col(aes(y = freq), color = "darkblue", fill = "lightblue") +
  geom_point(data = day7, aes(x = Number_2log, y = Frequency), color = "darkorange", size = 2) +
  labs(title = "Family sizes day 7",
       x = "Family size (2log)",
       y = "Frequency") +
  theme(plot.margin = margin(t = 1, r = 1, 0, 0, unit = "cm")) +
  theme_clean() 
plot_distribution_d7

plot_distribution_d8 <- ggplot(data = freq_famsizes_d8, aes(x = logfamsize)) + 
  geom_col(aes(y = freq), color = "darkblue", fill = "lightblue") +
  geom_point(data = day8, aes(x = Number_2log, y = Frequency), color = "darkorange", size = 2) +
  labs(title = "Family sizes day 8",
       x = "Family size (2log)",
       y = "Frequency") +
  theme(plot.margin = margin(t = 1, r = 1, 0, 0, unit = "cm")) +
  theme_clean() 
plot_distribution_d8

plot_grid(plotlist = list(plot_distribution_d5,
                       plot_distribution_d6,
                       plot_distribution_d7,
                       plot_distribution_d8), 
          align = "hv", rel_widths = c(0.5, 0.5), rel_heights = c(0.5, 0.5))
