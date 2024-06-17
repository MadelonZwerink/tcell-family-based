seed <- 4321
families = 500

bp_rule <- 'runif(1, min = 0, max = 3.5)'
dp_rule <- 0.5 #default
rq_rule <- 0.5 #default
t_start_dist <- 'rlnorm(nr_of_families, meanlog = 1.4, sdlog = 0.3)'
t_run_rule <- 'runif(1, min = 0.6, max = 4)'
nr_burst_divs <- 'sample(c(2,3,4), nr_of_families, replace = TRUE)'
quality_dist <- NULL
ASD <- FALSE
burst_time <- 0.15 #default
max_run_time <- NULL
min_t_start <- 2.5

prim_parameters <- pick_parameters(bp_rule = bp_rule,
                                   dp_rule = dp_rule,
                                   rq_rule = rq_rule,
                                   t_start_dist = t_start_dist,
                                   t_run_rule = t_run_rule,
                                   nr_of_families = families,
                                   nr_burst_divs = nr_burst_divs,
                                   response_nr = 1,
                                   quality_dist = quality_dist,
                                   ASD = ASD,
                                   burst_time = burst_time,
                                   max_run_time = max_run_time,
                                   min_t_start = min_t_start)


n <- 100
exp_divs <- prim_parameters$bp * (prim_parameters$t_run - (prim_parameters$t_burst * prim_parameters$nr_burst_divs))
div_index <- data.frame(matrix(nrow = nrow(prim_parameters), ncol = n))
div_cells <- data.frame(matrix(nrow = nrow(prim_parameters), ncol = n))
colnames(div_index) <- seq(0,n-1)
colnames(div_cells) <- seq(0,n-1)

for (fam in seq(nrow(prim_parameters))) {
  dist_divs <- dpois(x = seq(0, n), lambda = exp_divs[fam])
  nr_cells <- dist_divs * 2^seq(0, n)
  
  dist_cells_bursts <- c(rep(0, prim_parameters$nr_burst_divs[fam]), nr_cells)
  div_cells[fam,] <- dist_cells_bursts[0:n]
  
  dist_divs_bursts <- c(rep(0, prim_parameters$nr_burst_divs[fam]), dist_divs)
  div_index[fam,] <- dist_divs_bursts[0:n]
}

sum_cells <- data.frame(div_index = seq(0, n-1), 
                        nr_cells = colSums(div_cells), 
                        fraction = colSums(div_index, na.rm = TRUE)/nrow(prim_parameters))

ggplot(sum_cells, aes(x = div_index)) + 
  geom_bar(aes(y = log10(nr_cells)), stat = "identity") 

ggplot(sum_cells, aes(x = div_index)) + 
  geom_bar(aes(y = fraction), stat = "identity")

sum(sum_cells$nr_cells)
sum(sum_cells$nr_cells[0:10])

# This is to check whether the calculated number of cells matches the number of 
# cells that the function get_max_famsize returns. It is the same.

calculated_sum <- get_max_famsize(prim_parameters)
sum(calculated_sum)
