#-------------------------------------------------------------------------------

#FUNCTIONS RELATED TO GETTING STATISTICS

#-------------------------------------------------------------------------------

get_total_response_stats <- function(response){
  max_cells <- max(response$cells)
  time_max <- response$time[response$cells == max_cells]
  # Assume that the duration of the table is long enough so that towards the end
  # there's only Q cells left
  Q_cells <- round(response$cells[response$time == max(response$time)])
  
  stats <- data.frame(max_cells = as.numeric(max_cells), 
                      time_max = as.numeric(time_max), 
                      Q_cells = as.numeric(Q_cells))
  
  return(stats)
}

#-------------------------------------------------------------------------------

get_max_fam_stats <- function(max_cells, table = T){
  # Calculate Spearman correlation coefficients
  cor_prim_sec <- cor.test(max_cells$cells_prim, max_cells$cells_sec, method = "spearman")
  cor_sec_ter <- cor.test(max_cells$cells_sec, max_cells$cells_ter, method = "spearman")
  
  # Calculate mean and median fam size
  mean_prim <- mean(max_cells$cells_prim)
  mean_sec <- mean(max_cells$cells_sec)
  mean_ter <- mean(max_cells$cells_ter)
  median_prim <- median(max_cells$cells_prim)
  median_sec <- median(max_cells$cells_sec)
  median_ter <- median(max_cells$cells_ter)
  if(table == T){
    stats <- data.frame(mean = c(mean_prim, mean_sec, mean_ter),
                        median = c(median_prim, median_sec, median_ter),
                        correlation = c(cor_prim_sec$estimate, 0, cor_sec_ter$estimate),
                        p_value = c(cor_prim_sec$p.value, 0, cor_sec_ter$p.value))}
  else if(table == F){
    stats <- data.frame(mean_prim, mean_sec, mean_ter, 
                        median_prim, median_sec, median_ter,
                        cor_prim_sec = cor_prim_sec$estimate, p_value_prim_sec = cor_prim_sec$p.value,
                        cor_sec_ter = cor_sec_ter$estimate, p_value_sec_ter = cor_sec_ter$p.value)}
  
  return(stats)
}

#-------------------------------------------------------------------------------

get_famsize_stats <- function(df_famsizes){
  stat_mean <- round(mean(df_famsizes$famsize[df_famsizes$famsize > 1]), digits = 0)
  stat_sd <- round(sd(df_famsizes$famsize[df_famsizes$famsize > 1]), digits = 3)
  stat_median <- round(median(df_famsizes$famsize[df_famsizes$famsize > 1]), digits = 0)
  top_5perc_families <- nrow(df_famsizes) * 0.05
  stat_disparity <- round(sum(df_famsizes$fraction[1:top_5perc_families]), digits = 3)
  
  return(c(stat_mean, stat_median, stat_sd, stat_disparity))
}

#-------------------------------------------------------------------------------

get_famsize_stats_multidays <- function(famsizes_table){
  famsizes_stats <- lapply(5:8, function(day) {
    c(day, get_famsize_stats(famsizes_table[[day-4]]))
  })
  famsizes_stats <- data.frame(do.call(rbind, famsizes_stats))
  colnames(famsizes_stats) <- c("timepoint", "mean", "median", "sd", "disparity")
  
  return(famsizes_stats)
}

#-------------------------------------------------------------------------------

get_Q_famsize_stats <- function(max_fam_table){
  prim_cor <- cor.test(max_fam_table$cells_prim, max_fam_table$Q_cells, method = "spearman")
  sec_cor <- cor.test(max_fam_table$cells_sec, max_fam_table$Q_cells, method = "spearman")
  ter_cor <- cor.test(max_fam_table$cells_ter, max_fam_table$Q_cells, method = "spearman")
  
  correlations <- data.frame(r = c(prim_cor$estimate, sec_cor$estimate, ter_cor$estimate),
                             p_value = c(prim_cor$p.value, sec_cor$p.value, ter_cor$p.value))
  
  return(correlations)
}

#-------------------------------------------------------------------------------

get_Q_prim_famsize_stats <- function(Q_famsize_table){
  prim_cor <- cor.test(Q_famsize_table$famsize, Q_famsize_table$Q_cells, method = "spearman")
  
  return(prim_cor)
}

#-------------------------------------------------------------------------------