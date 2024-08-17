#This script makes a distribution plot for the start times

# Parameters for the lognormal distribution
logmean <- 1.4
logsd <- 0.3

# Generate data for the lognormal distribution
x <- seq(0.01, 10, length.out = 10000)
pdf <- dlnorm(x, meanlog = logmean, sdlog = logsd)

# Create a data frame for plotting
data <- data.frame(x = x, pdf = pdf)

# Plot the lognormal distribution
ggplot(data, aes(x = x, y = pdf)) +
  geom_line(linewidth = 1) +
  labs(title = "Start time distribution", x = "Days", y = "Probability Density") +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  theme_clean() + th
