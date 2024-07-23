
source("./R/speed_discrete_burst.R") 
library(fitdistrplus)

seed <- 4321
families = 10000
set.seed(seed)

bp_rule <- 'runif(1, min = 0, max = 3.5)'
dp_rule <- 0.5 #default
rq_rule <- 0.5 #default
t_start_dist <- 'rlnorm(nr_of_families, meanlog = 1.4, sdlog = 0.3)'
t_run_rule <- 'runif(1, min = 0.6, max = 4)' 
# was 0.5 in st2_nomax, but 0.6 ensures that all families can at least to their 4 burst 
# divisions, because 0.15*4 = 0.6
nr_burst_divs <- 'sample(c(2,3,4), nr_of_families, replace = TRUE)'
quality_dist <- NULL
quality_noise <- FALSE
q_noise_dist <- FALSE
uniform_fam <- FALSE
ASD <- FALSE
burst_time <- 0.15 #default
max_run_time <- NULL
min_t_start <- 0

prim_parameters <- pick_parameters(bp_rule = bp_rule,
                                   dp_rule = dp_rule,
                                   rq_rule = rq_rule,
                                   t_start_dist = t_start_dist,
                                   t_run_rule = t_run_rule,
                                   nr_of_families = families,
                                   nr_burst_divs = nr_burst_divs,
                                   response_nr = 1,
                                   quality_dist = quality_dist,
                                   quality_noise = quality_noise,
                                   q_noise_dist = q_noise_dist,
                                   uniform_fam = uniform_fam,
                                   ASD = ASD,
                                   burst_time = burst_time,
                                   max_run_time = max_run_time,
                                   min_t_start = min_t_start)

# Make dataframe that contains all division destinies for proliferating cells
dd_m1 <- prim_parameters$bp * (prim_parameters$t_run - prim_parameters$nr_burst_divs * prim_parameters$t_burst)
dd_m1 <- data.frame(fam_nr = prim_parameters$fam_nr, 
                    dd = dd_m1, 
                    cell_type = prim_parameters$cell_type)
dd_m1_P <- dd_m1[which(dd_m1$cell_type == "P"),]
dd_m1_P$cells <- exp(dd_m1_P$dd)

dd_per_fam <- dd_m1_P %>%
  group_by(fam_nr) %>%
  summarise(fam_cells = mean(cells))

dd_per_fam$fam_dd <- log(dd_per_fam$fam_cells)
max_dd <- max(dd_per_fam$fam_dd) #12.564
# Normalize the dd per family so that it is a value between 0 and 1
dd_per_fam$dd_norm <- dd_per_fam$fam_dd / max_dd

ggplot(dd_per_fam, aes(x = fam_dd)) + 
  geom_histogram(binwidth = 1) +
  labs(title = "Average DD per family (Q branches excluded)",
       x = "DD per family",
       y = "Nr. of families") +
  theme_clean() + th +
  scale_x_continuous(breaks = c(0,4,8,12), limits = c(0, 13))

ggsave("./results/model2/methods_DD_dist_e_m1.jpg", 
       width = 2000, height = 1500, units = "px")

dist_stats <- descdist(dd_per_fam$dd_norm, discrete = F)
mean_val <- dist_stats$mean
var_val <- dist_stats$sd^2

# Calculate alpha and beta using method of moments
term <- (mean_val * (1 - mean_val) / var_val) - 1
alpha <- mean_val * term
beta <- (1 - mean_val) * term

# Output the estimated parameters
alpha # 2.3332 (n = 5000) # 2.29 (n = 10000)
beta  # 3.9193 (n = 5000) # 3.80 (n = 10000)

# This will be the distribution that we use to get division indexes
dd_dist <- data.frame(dd = seq(0, max_dd, length.out = 1000),
                      density = dbeta(seq(0,1,length.out = 1000), alpha, beta))

ggplot(dd_dist, aes(x = dd, y = density)) +
  geom_point() +
  theme_clean() +
  th +
  labs(title = "Fitted DD distribution",
       x = "DD per family",
       y = "PDF") +
  scale_x_continuous(breaks = c(0,4,8,12), limits = c(0, 13))

ggsave("./results/model2/methods_fitted_DD_dist_e.jpg", 
       width = 2000, height = 1500, units = "px")

#-------------------------------------------------------------------------------

#Take the mean division rate of all proliferation branches for this model and 
#calculate the run time accordingly

mean(prim_parameters$bp[prim_parameters$cell_type == "P"])

