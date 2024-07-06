dd_m1 <- prim_parameters$bp * (prim_parameters$t_run - prim_parameters$nr_burst_divs * prim_parameters$t_burst)
P_cells <- which(prim_parameters$cell_type == "P")
hist(dd_m1[P_cells], freq = F)
lines(density(dd_m1[P_cells]))
plot_m1 <- ggplot(as.data.frame(dd_m1[P_cells])) + geom_histogram(aes(x = dd_m1[P_cells]), binwidth = 0.5)


# Test 

p <- runif(5000, min = 0, max = 3.5)
r <- runif(5000, min = 0.6, max = 4)
b <- sample(c(2,3,4), prob = c(1/7, 8/28, 16/28), 5000, replace = TRUE)
dd <- p * (r - b*0.15)
hist(dd, freq = F)
lines(density(dd))
plot_toy <- ggplot(as.data.frame(dd)) + geom_histogram(aes(x = dd), binwidth = 0.2)


# Uniform quality

q <- runif(5000, min = 0, max = 1)

p <- q * 3.5
b <- sample(c(3.7, 3.55, 3.4), prob = c(1/7, 8/28, 16/28), 5000, replace = TRUE)
r <- q * b
dd2 <- p * r
hist(dd2, freq = F)
lines(density(dd2))
plot_q <- ggplot(as.data.frame(dd2)) + geom_histogram(aes(x = dd2), binwidth = 0.2)

# Ranges

range_p <- seq(0, 3.5, length.out = 5000)
range_r <- seq(0, 3.4, length.out = 5000)
range_q <- seq(0, 1, length.out = 5000)
range_dd <- range_p * range_r

range_df <- data.frame(p = range_p, r = range_r, q = range_q, dd = range_dd)

get_closest_index <- function(value){
  which.min(abs(value - range_df$dd))
}
indices <- vapply(dd, get_closest_index, numeric(1))
df_corrected <- range_df[indices,]
ggplot(df_corrected, aes(x = dd)) + geom_histogram()
ggplot(df_corrected, aes(x = q)) + geom_histogram()

plot_range <- ggplot(as.data.frame(range_dd)) + geom_histogram(aes(x = range_dd), binwidth = 0.2)

dd_mix <- range_r * p
plot_mix <- ggplot(as.data.frame(dd_mix)) + geom_histogram(aes(x = dd_mix), binwidth = 0.2)

df_mix <- data.frame(dd = dd_mix, q = range_r/3.4)
ggplot(df_mix) + geom_point(aes(x = dd, y = q))

plot_grid(plotlist = list(plot_m1, plot_toy, plot_q, plot_range, plot_mix), nrow = 2)

value <- dd[1]

# Beta distributed quality (alpha = 1.6, beta = 1.85)

q <- rbeta(5000, shape1 = 1.6, shape2 = 1.85)

p <- q * 3.5
b <- sample(c(3.7, 3.55, 3.4), prob = c(1/7, 8/28, 16/28), 5000, replace = TRUE)
r <- q * b
dd2 <- p * r
hist(dd2, freq = F)
lines(density(dd2))
plot_q <- ggplot(as.data.frame(dd2)) + geom_histogram(aes(x = dd2), binwidth = 0.2)
