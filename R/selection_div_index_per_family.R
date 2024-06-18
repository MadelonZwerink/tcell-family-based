# Open this, to get the division index distribution for each family

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

div_index_prim <- generate_div_index_df(prim_parameters)
sum_index_prim <- data.frame(div_index = seq(0, n-1),
                             fraction = colSums(div_index_prim, na.rm = TRUE)/nrow(prim_parameters))

div_index_prim_long <- wide_to_long(div_index_prim, "fraction")

div_index_prim_grouped <- div_index_prim_long %>%
  group_by(fam_nr, div_index) %>%
  summarise(fraction = sum(fraction)/n())

ggplot(div_index_prim_grouped, aes(x = div_index, 
                                   y = fraction,
                                   col = as.factor(fam_nr), 
                                   group = fam_nr,
                                   alpha = 0.5)) +
  geom_line() +
  theme(legend.position = "none") + xlim(0, 25)

#-------------------------------------------------------------------------------

exp_divs_prim <- get_exp_divs(prim_parameters)
exp_divs_prim <- data.frame(fam_nr = prim_parameters$fam_nr, exp_divs = (exp_divs_prim + as.numeric(prim_parameters$div_counter)))

ggplot(exp_divs_prim, aes(x = exp_divs, y = as.factor(fam_nr), 
                          fill = as.factor(fam_nr))) +
  geom_violin() +
  geom_point(shape = 21) +
  theme(legend.position = "none")

#-------------------------------------------------------------------------------

sum_cells <- data.frame(div_index = seq(0, n-1), 
                        nr_cells = colSums(div_cells, na.rm = TRUE), 
                        fraction = colSums(div_index, na.rm = TRUE)/nrow(parameters))
