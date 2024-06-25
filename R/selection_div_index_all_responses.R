# Open this, to get the division index for the whole population

seed <- 4321
families = 500
n = 100

bp_rule <- 'runif(1, min = 0, max = 3.5)'
dp_rule <- 0.5 #default
rq_rule <- 1 #means no Q cells are formed
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

#-------------------------------------------------------------------------------

generate_div_index_df <- function(parameters, n = 100){
  exp_divs <- parameters$bp * (parameters$t_run - (parameters$t_burst * parameters$nr_burst_divs))
  div_index <- data.frame(matrix(nrow = nrow(parameters), ncol = n))
  div_cells <- data.frame(matrix(nrow = nrow(parameters), ncol = n))
  colnames(div_index) <- seq(0,n-1)
  colnames(div_cells) <- seq(0,n-1)
  
  for (fam in seq(nrow(parameters))) {
    dist_divs <- dpois(x = seq(0, n), lambda = exp_divs[fam])
    nr_cells <- dist_divs * 2^seq(0, n)
    
    dist_cells_bursts <- c(rep(0, parameters$div_counter[fam]), nr_cells)
    div_cells[fam,] <- dist_cells_bursts[0:n]
    
    dist_divs_bursts <- c(rep(0, parameters$div_counter[fam]), dist_divs)
    div_index[fam,] <- dist_divs_bursts[0:n]
  }
  
  div_cells <- cbind(fam_nr = parameters$fam_nr, div_cells)
  div_index <- cbind(fam_nr = parameters$fam_nr, div_index)
  
  return(list(div_cells, div_index))
}

sum_cells <- data.frame(div_index = seq(0, n-1), 
                        nr_cells = colSums(div_cells, na.rm = TRUE), 
                        fraction = colSums(div_index, na.rm = TRUE)/nrow(parameters))

#-------------------------------------------------------------------------------

sum_cells_prim_2 <- generate_div_index_df(prim_parameters)
sum_cells_sec <- generate_div_index_df(sec_parameters)
sum_cells_ter <- generate_div_index_df(ter_parameters)

#-------------------------------------------------------------------------------
# Primary response
#-------------------------------------------------------------------------------

ggplot(sum_cells_prim_2, aes(x = div_index)) + 
  geom_bar(aes(y = log10(nr_cells)), stat = "identity") +
  labs(title = "Total number of cells for each division index",
       x = "Division index",
       y = "log10(number of cells)") +
  theme_clean() + th + xlim(0, 50)

ggplot(sum_index_prim, aes(x = div_index)) + 
  geom_bar(aes(y = fraction), stat = "identity") +
  labs(title = "Fraction of cells at each division index",
       subtitle = "Corrected for increase in cell numbers with each division",
       x = "Division index",
       y = "Fraction") +
  theme_clean() + th + xlim(0, 50)

sum(sum_cells_prim$nr_cells)
sum(sum_cells_prim$nr_cells[0:5])

# This is to check whether the calculated number of cells matches the number of 
# cells that the function get_max_famsize returns. It is the same.

calculated_sum <- get_max_famsize(prim_parameters)
sum(calculated_sum)

#-------------------------------------------------------------------------------
# Secondary response
#-------------------------------------------------------------------------------

ggplot(sum_cells_sec, aes(x = div_index)) + 
  geom_bar(aes(y = log10(nr_cells)), stat = "identity") +
  labs(title = "Total number of cells for each division index",
       x = "Division index",
       y = "log10(number of cells)") +
  theme_clean() + th + xlim(0, 50)

ggplot(sum_cells_sec, aes(x = div_index)) + 
  geom_bar(aes(y = fraction), stat = "identity") +
  labs(title = "Fraction of cells at each division index",
       subtitle = "Corrected for increase in cell numbers with each division",
       x = "Division index",
       y = "Fraction") +
  theme_clean() + th + xlim(0, 50)

#-------------------------------------------------------------------------------
# Tertiary response
#-------------------------------------------------------------------------------

ggplot(sum_cells_ter, aes(x = div_index)) + 
  geom_bar(aes(y = log10(nr_cells)), stat = "identity") +
  labs(title = "Total number of cells for each division index",
       x = "Division index",
       y = "log10(number of cells)") +
  theme_clean() + th + xlim(0, 50)

ggplot(sum_cells_ter, aes(x = div_index)) + 
  geom_bar(aes(y = fraction), stat = "identity") +
  labs(title = "Fraction of cells at each division index",
       subtitle = "Corrected for increase in cell numbers with each division",
       x = "Division index",
       y = "Fraction") +
  theme_clean() + th + xlim(0, 50)

#-------------------------------------------------------------------------------

exp_divs_sec <- sec_parameters$bp * (sec_parameters$t_run - (sec_parameters$t_burst * sec_parameters$nr_burst_divs))
div_index_sec <- data.frame(matrix(nrow = nrow(sec_parameters), ncol = n))
div_cells_sec <- data.frame(matrix(nrow = nrow(sec_parameters), ncol = n))
colnames(div_index_sec) <- seq(0,n-1)
colnames(div_cells_sec) <- seq(0,n-1)

for (fam in seq(nrow(sec_parameters))) {
  dist_divs_sec <- dpois(x = seq(0, n), lambda = exp_divs_sec[fam])
  nr_cells_sec <- dist_divs_sec * 2^seq(0, n)
  
  dist_cells_bursts <- c(rep(0, sec_parameters$div_counter[fam]), nr_cells_sec)
  div_cells_sec[fam,] <- dist_cells_bursts[0:n]
  
  dist_divs_bursts <- c(rep(0, sec_parameters$div_counter[fam]), dist_divs_sec)
  div_index_sec[fam,] <- dist_divs_bursts[0:n]
}

sum_cells_sec <- data.frame(div_index = seq(0, n-1), 
                        nr_cells = colSums(div_cells_sec, na.rm = TRUE), 
                        fraction = colSums(div_index_sec, na.rm = TRUE)/nrow(sec_parameters))

ggplot(sum_cells_sec, aes(x = div_index)) + 
  geom_bar(aes(y = log10(nr_cells)), stat = "identity") +
  labs(title = "Total number of cells for each division index",
       x = "Division index",
       y = "log10(number of cells)") +
  theme_clean() + th + xlim(0, 50)

ggplot(sum_cells_sec, aes(x = div_index)) + 
  geom_bar(aes(y = fraction), stat = "identity") +
  labs(title = "Fraction of cells at each division index",
       subtitle = "Corrected for increase in cell numbers with each division",
       x = "Division index",
       y = "Fraction") +
  theme_clean() + th + xlim(0, 50)

sum(sum_cells$nr_cells)
sum(sum_cells$nr_cells[0:5])