
#-------------------------------------------------------------------------------

# Set the number of simulations
n_simulations <- 1000

# Sample quality variable
quality <- runif(n_simulations, min = 0.1, max = 1)

# Function to skew a uniform distribution towards higher values based on quality
skewed_uniform <- function(n, min_val, max_val, quality) {
  # Beta distribution parameters
  alpha <- 1 + quality
  beta <- 1 / quality # Larger quality gives more skew towards higher values
  
  # Sample from beta distribution and transform to uniform distribution
  beta_sample <- rbeta(n, alpha, beta)
  return(min_val + (max_val - min_val) * beta_sample)
}

# Sample proliferation rate and run time based on quality
proliferation_rate <- skewed_uniform(n_simulations, 0, 3.5, quality)
run_time <- skewed_uniform(n_simulations, 0.6, 4, quality)

# Combine into a data frame
simulation_data <- data.frame(proliferation_rate, run_time, quality)

# Display first few rows of the data
head(simulation_data)

ggplot(simulation_data, aes(x = run_time, y = proliferation_rate, size = quality)) + 
  geom_point()

#-------------------------------------------------------------------------------
#First make sure the skewed uniform function from speed_discrete_burst.R is loaded

get_alpha <- function(quality){
  alpha <- 1 / (1 - quality)
  return(alpha)
}

get_beta <- function(quality){
  beta <- 1 / quality
  return(beta)
}

p1 <- generate_quality_effect_panelplot(0.1, 0.2, 0.3)
p2 <- generate_quality_effect_panelplot(0.4, 0.5, 0.6)
p3 <- generate_quality_effect_panelplot(0.7, 0.8, 0.9)
cowplot::plot_grid(plotlist = list(p1, p2, p3), nrow = 3)
