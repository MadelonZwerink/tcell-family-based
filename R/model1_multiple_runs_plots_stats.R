# Output from cluster:

## DATA
### df_famsizes
### total_resp
### max_fam_table
### famsizes_freq_table

## STATISTICS
### famsize_stats
### max_fam_stats_long
### total_response_stats

#-------------------------------------------------------------------------------
# What do I need to make?

## PLOTS
### Response plot: small line for each simulation with mean, embed in panelfig.
### Famsize distribution plot: famsize distributions with conf int (use old code)

# STATS
### Average correlations Q cells and famsize
### Average correlations prim-sec and sec-ter
### Average number of Q cells
### Average response size

#-------------------------------------------------------------------------------
# Provide information specific for run
folder <- "./results/model1/iterative_runs/model1/"
run_name <- "st3_multiple_"

source("./R/functions_multiple_sims.R") 

# Load data
data_objects <- c("df_famsizes", "total_resp", "max_fam_table", "famsizes_freq_table",
                  "famsizes_stats", "max_fam_stats_long", "prim_response_stats",
                  "sec_response_stats", "ter_response_stats")
for(data_object in data_objects){
  assign(data_object, readRDS(paste0(folder, run_name, data_object, ".rds"))) 
}

#-------------------------------------------------------------------------------
# Response plot

total_resp_combined <- total_resp %>%
  imap_dfr(~mutate(.x, source = .y)) %>%
  arrange(source, time)

row.names(total_resp_combined) <- NULL

total_resp_mean <- total_resp_combined %>%
  group_by(time) %>%
  summarise(mean_log_cells = mean(log_cells),
            sd_log_cells = sd(log_cells),
            lower_bound = mean(log_cells) - 1.96 * sd(log_cells)/sqrt(n()),
            upper_bound = mean(log_cells) + 1.96 * sd(log_cells)/sqrt(n()),
            .groups = "drop")

response_plot_multi <- ggplot(data = total_resp_combined, aes(x = time, y = log_cells, group = source)) +
  geom_line(size = 0.01) +
  geom_line(data = total_resp_mean, aes(x = time, y = mean_log_cells), 
            color = "#9f2a63", size = 0.2, inherit.aes = F) +  # Add mean line
  geom_ribbon(data = total_resp_mean, aes(x = time, 
                                          ymin = lower_bound, 
                                          ymax = upper_bound), 
              fill = "#9f2a63", alpha = 0.5, inherit.aes = F) +
#  geom_ribbon(data = total_resp_mean, aes(x = time, ymin = log10(mean_cells)))
  labs(x = "Time (days)",
       y = "Nr. of cells (log-scale)") + 
  theme_clean() + 
  th

#-------------------------------------------------------------------------------
# Response stats

processed_total_resp_stats <- data.frame(rbind(get_processed_stats(prim_response_stats),
                                               get_processed_stats(sec_response_stats),
                                               get_processed_stats(ter_response_stats)))

#-------------------------------------------------------------------------------
# Famsize distribution plot
# Use famsizes_freq_table, which is a list of dataframes that contains a dataframe
# for each day. Each row represents a simulation and the columns indicate the fraction 
# of cells in each bin

processed_famsizes_freq_stats <- lapply(famsizes_freq_table, get_processed_stats)

famsize_dist_v_plot_multi <- 
  plot_v_grid_famsize_dist(frequencies_table = processed_famsizes_freq_stats, 
                           multi_plot = T,
                           x_axis_max = rep(19, 4))

ggsave(paste0(folder, run_name, "_famsize_dist_v_plot_multi.pdf"), plot = famsize_dist_v_plot_multi, 
       width = 5, height = 8, units = "in")

ggsave(paste0(folder, run_name, "_famsize_dist_v_plot_multi.jpg"), plot = famsize_dist_v_plot_multi, 
       width = 5, height = 8, units = "in")

#-------------------------------------------------------------------------------
# Famsize statistics
famsizes_stats <- lapply(famsizes_stats, convert_cols_to_numeric)
processed_famsize_stats <- lapply(famsizes_stats, get_processed_stats)

#-------------------------------------------------------------------------------
# Maxfam stats
max_fam_stats_long <- convert_cols_to_numeric(max_fam_stats_long)

max_fam_boxplot_cor_data <- max_fam_stats_long %>%
  select(cor_prim_sec, cor_sec_ter) %>%
  pivot_longer(cols = all_of(c("cor_prim_sec", "cor_sec_ter")),
               names_to = "metric",
               values_to = "value")

max_fam_boxplot_pvalue_data <- max_fam_stats_long %>%
  select(p_value_prim_sec, p_value_sec_ter) %>%
  mutate(p_value_sec_ter = p_value_sec_ter * 10^14) %>%
  pivot_longer(cols = all_of(c("p_value_prim_sec", "p_value_sec_ter")),
               names_to = "metric",
               values_to = "value")
  
boxplot_correlations <- ggplot(max_fam_boxplot_cor_data) + 
  geom_boxplot(aes(y = value, x = metric)) +
  theme_clean() + th +
  labs(
    y = "Spearman r",
    x = element_blank()) +
  scale_x_discrete(labels = c("prim vs sec", "sec vs ter"))
  
coeff <- 10^14
boxplot_p_values <- ggplot(max_fam_boxplot_pvalue_data) + 
  geom_boxplot(aes(y = value, x = metric)) +
  theme_clean() + th +
  labs(
    y = "P-value",
    x = element_blank()) +
  scale_x_discrete(labels = c("prim vs sec", "sec vs ter"))  +
  scale_y_continuous(
    # Features of the first axis
    name = "P-value",
    # Add a second axis and specify its features
    sec.axis = sec_axis(trans=~./coeff, name="P-value", 
                        breaks = c(0, 5*10^-15, 10^-14))
  )

processed_max_fam_stats <- get_processed_stats(max_fam_stats_long)

#-------------------------------------------------------------------------------
# Q famsize stats
Q_famsize_stats <- sapply(max_fam_table, get_Q_famsize_stats_multi)
Q_famsize_stats <- as.data.frame(t(Q_famsize_stats))
colnames(Q_famsize_stats) <- c("prim_cor", "sec_cor", "ter_cor", 
                               "prim_pvalue", "sec_pvalue", "ter_pvalue")
processed_Q_famsize_stats <- get_processed_stats(Q_famsize_stats)
#-------------------------------------------------------------------------------

output_data <- c("processed_total_resp_stats", "processed_max_fam_stats", 
                 "processed_famsizes_freq_stats", "processed_famsize_stats", 
                 "processed_Q_famsize_stats")

sink(file = paste0(folder, run_name, "_output.txt"))

for(output_stats in output_data){
  print(eval(output_stats))
  print(get(eval(output_stats)))
}

sink(file = NULL)

ggsave(paste0(folder, run_name, "_famsize_correlations_boxplot.jpg"), 
       plot = boxplot_correlations, 
       width = 4, height = 6, units = "in")

ggsave(paste0(folder, run_name, "_famsize_cor_pvalues_boxplot.jpg"), 
       plot = boxplot_p_values, 
       width = 4, height = 6, units = "in")


