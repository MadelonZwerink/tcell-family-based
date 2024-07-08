
library(fitdistrplus)
n = 100000

#Generate uniform quality and calculate corresponding proliferation rate and 
#run time 
range_q <- seq(0, 1, length.out = n)
range_p <- range_q * 3.5
range_r <- range_q * 3.4
range_dd <- range_p * range_r

range_df <- data.frame(p = range_p, r = range_r, q = range_q, dd = range_dd)

get_closest_index <- function(value){
  which.min(abs(value - range_df$dd))
}

#-------------------------------------------------------------------------------
#Calculates the necessary distribution of qualities based on the DD per branch
#-------------------------------------------------------------------------------

#Generate proliferation rate and run time just like in model 1
p <- runif(n, min = 0, max = 3.5)
r <- runif(n, min = 0.6, max = 4)
b <- sample(c(2,3,4), prob = c(1/7, 8/28, 16/28), n, replace = TRUE)
dd <- p * (r - b*0.15)

#Determine which division destiny is closest to the toy model dd
#To basically reverse engineer the distribution of quality that would be 
#necessary to get the same distribution
indices <- vapply(dd, get_closest_index, numeric(1))
df_corrected <- range_df[indices,]

ggplot(df_corrected, aes(x = dd)) + geom_histogram()
ggplot(df_corrected, aes(x = q)) + geom_histogram()

hist(df_corrected$q, freq = F)
lines(density(df_corrected$q))

descdist(df_corrected$q, discrete = F)

# Given summary statistics
mean_val <- 0.4598126
var_val <- 0.2280213^2

# Calculate alpha and beta using method of moments
term <- (mean_val * (1 - mean_val) / var_val) - 1
alpha <- mean_val * term
beta <- (1 - mean_val) * term

# Output the estimated parameters
alpha
beta

#-------------------------------------------------------------------------------
#Calculates the necessary distribution of qualities based on the DD per family
#-------------------------------------------------------------------------------

# Make dataframe that contains all division destinies for proliferating cells
dd_m1 <- prim_parameters$bp * (prim_parameters$t_run - prim_parameters$nr_burst_divs * prim_parameters$t_burst)
dd_m1 <- data.frame(fam_nr = prim_parameters$fam_nr, 
                    dd = dd_m1, 
                    cell_type = prim_parameters$cell_type)
dd_m1_P <- dd_m1[which(dd_m1$cell_type == "P"),]

dd_per_fam <- dd_m1_P %>%
  group_by(fam_nr) %>%
  summarise(fam_dd = mean(dd))

indices <- vapply(dd_per_fam$fam_dd, get_closest_index, numeric(1))
df_corrected <- range_df[indices,]

ggplot(df_corrected, aes(x = dd)) + 
  geom_histogram() +
  labs(title = "Average DD per family (Q branches excluded)")
ggplot(df_corrected, aes(x = q)) + 
  geom_histogram(binwidth = 0.05) + 
  xlim(0, 1) + 
  labs(title = "Quality distribution based on average DD per family (model 1)")

hist(df_corrected$q, freq = F)
lines(density(df_corrected$q))  

descdist(df_corrected$q, discrete = F)

fitdist(df_corrected$q, distr = "norm")
plot(dnorm(seq(0,1,0.01), mean = 0.518, sd = 0.127))

#-------------------------------------------------------------------------------
#Calculates the necessary distribution of qualities based on the weighted DD per family
#-------------------------------------------------------------------------------

dd_m1_P$cells <- 2^dd_m1_P$dd

weighted_dd_per_fam <- dd_m1_P %>%
  group_by(fam_nr) %>%
  summarise(fam_cells = mean(cells))

weighted_dd_per_fam$fam_dd <- log2(weighted_dd_per_fam$fam_cells)


indices <- vapply(weighted_dd_per_fam$fam_dd, get_closest_index, numeric(1))
df_corrected <- range_df[indices,]

ggplot(df_corrected, aes(x = dd)) + 
  geom_histogram() +
  labs(title = "Average DD per family (Q branches excluded)")
ggplot(df_corrected, aes(x = q)) + 
  geom_histogram(binwidth = 0.05) + 
  xlim(0, 1) + 
  labs(title = "Quality distribution based on average DD per family (model 1)")

hist(df_corrected$q, freq = F)
lines(density(df_corrected$q))  

descdist(df_corrected$q, discrete = F)

# Given summary statistics
mean_val <- 0.62027
var_val <- 0.1627262^2

# Calculate alpha and beta using method of moments
term <- (mean_val * (1 - mean_val) / var_val) - 1
alpha <- mean_val * term
beta <- (1 - mean_val) * term

# Output the estimated parameters
alpha
beta
