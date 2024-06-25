# Open this, to get the division index distribution for each family

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#FUNCTIONS
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

get_exp_divs <- function(parameters){
  exp_divs <- parameters$bp * (parameters$t_run - (parameters$t_burst * parameters$nr_burst_divs))
  return(exp_divs)
}

generate_div_index_df <- function(parameters, n = 100){
  exp_divs <- get_exp_divs(parameters)
  
  div_index <- data.frame(matrix(nrow = nrow(parameters), ncol = n))
  colnames(div_index) <- seq(0,n-1)
  
  for (fam in seq(nrow(parameters))) {
    dist_divs <- dpois(x = seq(0, n), lambda = exp_divs[fam])
    dist_divs_bursts <- c(rep(0, parameters$div_counter[fam]), dist_divs)
    div_index[fam,] <- dist_divs_bursts[0:n]
  }
  
  div_index <- cbind(fam_nr = parameters$fam_nr, div_index)
  
  return(div_index)
}

generate_div_nr_cells_df <- function(parameters, n = 100){
  exp_divs <- get_exp_divs(parameters)
  
  div_cells <- data.frame(matrix(nrow = nrow(parameters), ncol = n))
  colnames(div_cells) <- seq(0,n-1)
  
  for (fam in seq(nrow(parameters))) {
    dist_divs <- dpois(x = seq(0, n), lambda = exp_divs[fam])
    nr_cells <- dist_divs * 2^seq(0, n)
    
    dist_cells_bursts <- c(rep(0, parameters$div_counter[fam]), nr_cells)
    div_cells[fam,] <- dist_cells_bursts[0:n]
  }
  
  div_cells <- cbind(fam_nr = parameters$fam_nr, div_cells)
  
  return(div_cells)
}

wide_to_long <- function(div_df, names_values) {
  div_df %<>% pivot_longer(cols = 2:ncol(div_df), 
                                names_to = "div_index", 
                                values_to = eval(names_values))
  div_df$div_index <- as.numeric(div_df$div_index)
  return(div_df)
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Calculate the fraction for each division index

div_index_prim <- generate_div_index_df(prim_parameters)
sum_index_prim <- data.frame(div_index = seq(0, n - 1),
                             fraction = colSums(div_index_prim[,2:ncol(div_index_prim)], na.rm = TRUE)/nrow(prim_parameters))

div_index_prim_grouped <- div_index_prim %>%
  wide_to_long("fraction") %>%
  group_by(fam_nr, div_index) %>%
  summarise(fraction = sum(fraction)/n())

# Plots division index per family
ggplot(div_index_prim_grouped, aes(x = div_index, 
                                   y = fraction,
                                   col = as.factor(fam_nr), 
                                   group = fam_nr,
                                   alpha = 0.5)) +
  labs(title = "Fraction of cells at each division index per family",
       subtitle = "Corrected for increase in cell numbers with each division",
       x = "Division index",
       y = "Fraction") +
  geom_line() +
  theme_clean() + th +
  theme(legend.position = "none") + xlim(0, 25)

# Plots total distribution of division indices
ggplot(sum_index_prim, aes(x = div_index)) + 
  geom_bar(aes(y = fraction), stat = "identity") +
  labs(title = "Fraction of cells at each division index",
       subtitle = "Corrected for increase in cell numbers with each division",
       x = "Division index",
       y = "Fraction") +
  theme_clean() + th + xlim(0, 25)

#-------------------------------------------------------------------------------
# Calculate the absolute number of cells for each division index

div_cells_prim <- generate_div_nr_cells_df(prim_parameters)

sum_cells_prim <- data.frame(div_index = seq(0, n - 1),
                             nr_cells = colSums(div_cells_prim[,2:ncol(div_cells_prim)], na.rm = TRUE))

div_cells_prim_grouped <- div_cells_prim %>%
  wide_to_long("nr_cells") %>%
  group_by(fam_nr, div_index) %>%
  summarise(nr_cells = sum(nr_cells))

# Plots division index per family
ggplot(div_cells_prim_grouped, aes(x = div_index, 
                                   y = nr_cells,
                                   col = as.factor(fam_nr), 
                                   group = fam_nr,
                                   alpha = 0.5)) +
  labs(title = "Number of cells at each division index per family",
       x = "Division index",
       y = "Number of cells") +
  geom_line() +
  theme_clean() + th +
  theme(legend.position = "none") + xlim(0, 50)

# Plots total distribution of division indices
ggplot(sum_cells_prim, aes(x = div_index)) + 
  geom_bar(aes(y = log10(nr_cells)), stat = "identity") +
  labs(title = "Number of cells at each division index",
       x = "Division index",
       y = "Number of cells (log-scale)") +
  theme_clean() + th + xlim(0, 50)

#-------------------------------------------------------------------------------
# Plots the expected divisions for each family in a violin plot
# Not really useful because it is hard to interpret for runs with a lot of families

exp_divs_prim <- get_exp_divs(prim_parameters)
exp_divs_prim <- data.frame(fam_nr = prim_parameters$fam_nr, exp_divs = (exp_divs_prim + as.numeric(prim_parameters$div_counter)))

ggplot(exp_divs_prim, aes(x = exp_divs, y = as.factor(fam_nr), 
                          fill = as.factor(fam_nr))) +
  geom_violin() +
  geom_point(shape = 21) +
  theme(legend.position = "none")

#-------------------------------------------------------------------------------
