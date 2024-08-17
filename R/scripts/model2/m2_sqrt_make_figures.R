
folder <- "./results/model2/multiple_runs/"
run_name <- "m2_sqrt_100"

source("./R/functions_multiple_sims.R") 

#-------------------------------------------------------------------------------
# Load data
#-------------------------------------------------------------------------------
#DATA FROM ITERATIVE RUN:

data_objects <- c("df_famsizes", 
                  "total_resp", 
                  "max_fam_table", 
                  "famsizes_freq_table",
                  "famsizes_stats", 
                  "max_fam_stats_long", 
                  "prim_response_stats",
                  "sec_response_stats", 
                  "ter_response_stats")
for(data_object in data_objects){
  assign(data_object, readRDS(paste0(folder, run_name, "_", data_object, ".rds"))) 
}

#DATA FROM SINGLE RUN:
seed <- 4321
families <- 500
nr_of_families <- 500

set.seed(seed)

bp_rule <- 'sqrt(q*11.13)'
dp_rule <- 0.5 #default
rq_rule <- 0.5 #default
t_start_dist <- 'rlnorm(nr_of_families, meanlog = 1.4, sdlog = 0.3)'
t_run_rule <- 'sqrt(q*11.13) + (nr_burst_divs[i]*t_burst[i])'
nr_burst_divs <- 'sample(c(2,3,4), nr_of_families, replace = TRUE)'
quality_dist <- "rbeta(nr_of_families, shape1 = 2.165, shape2 = 3.295)"
quality_noise <- FALSE
uniform_fam <- TRUE
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

sec_parameters <- pick_parameters(response_nr = 2,
                                  prev_parameters = prim_parameters,
                                  bp_rule = bp_rule,
                                  dp_rule = dp_rule,
                                  rq_rule = rq_rule,
                                  t_start_dist = t_start_dist,
                                  t_run_rule = t_run_rule,
                                  nr_burst_divs = nr_burst_divs,
                                  quality_dist = quality_dist,
                                  quality_noise = quality_noise,
                                  q_noise_dist = q_noise_dist,
                                  uniform_fam = uniform_fam,
                                  ASD = ASD,
                                  burst_time = burst_time,
                                  max_run_time = max_run_time,
                                  min_t_start = min_t_start)

ter_parameters <- pick_parameters(response_nr = 3,
                                  prev_parameters = sec_parameters,
                                  bp_rule = bp_rule,
                                  dp_rule = dp_rule,
                                  rq_rule = rq_rule,
                                  t_start_dist = t_start_dist,
                                  t_run_rule = t_run_rule,
                                  nr_burst_divs = nr_burst_divs,
                                  quality_dist = quality_dist,
                                  quality_noise = quality_noise,
                                  q_noise_dist = q_noise_dist,
                                  uniform_fam = uniform_fam,
                                  ASD = ASD,
                                  burst_time = burst_time,
                                  max_run_time = max_run_time,
                                  min_t_start = min_t_start)

max_fam_table_single <- generate_max_fam_table(prim_parameters, 
                                               sec_parameters, 
                                               ter_parameters, 
                                               type = "time", 
                                               timepoint = 7)

#-------------------------------------------------------------------------------
#Table 2: distributions of family sizes
#-------------------------------------------------------------------------------
famsizes_stats <- lapply(famsizes_stats, convert_cols_to_numeric)
processed_famsize_stats <- lapply(famsizes_stats, get_processed_stats)

sink(file = paste0(folder, run_name, "_T2_famsizes_distribution.txt"))

print("Table 2: famsize distributions for d5-8")
print(processed_famsize_stats)

sink(file = NULL)

#-------------------------------------------------------------------------------
#Figure 3: distributions of family sizes
#-------------------------------------------------------------------------------
processed_famsizes_freq_stats <- lapply(famsizes_freq_table, get_processed_stats)

famsize_dist_plot <-
  plot_grid_famsize_dist(frequencies_table = processed_famsizes_freq_stats, 
                         multi_plot = T,
                         x_axis_max = rep(19, 4))

ggsave(paste0(folder, run_name, "_3_famsize_dist_plot.jpg"), plot = famsize_dist_plot, 
       width = 5, height = 8, units = "in")

#-------------------------------------------------------------------------------
#Figure 4: disparity of family sizes
#-------------------------------------------------------------------------------
#A
largest_fam <- as.data.frame(t(sapply(max_fam_table, get_largest_fam))) %>%
  pivot_longer(cols = everything(),
               names_to = "metric",
               values_to = "value")

largest_fam_plot <- ggplot(largest_fam, aes(y = log10(value), x = metric)) + 
  geom_jitter() +
  stat_summary(fun=mean, geom="crossbar", colour="#9f2a63") +
  theme_clean() + th +
  labs(x = element_blank(),
       y = "Family size (cell number, 10log-scale)")  +
  scale_x_discrete(labels = c("Largest family", "Median family")) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 6)) +
  theme(plot.margin = margin(l = 1, b = 0.8, t = 0.6, unit = "cm"))

#B
famsizes_prim <- as.data.frame(sapply(max_fam_table, get_famsize_fraction))
# Each column is a simulation, each row is a family, they are ordered from big to small

twoperc_resp <- colSums(famsizes_prim[0:(0.02*nrow(famsizes_prim)),])
fiveperc_resp <- colSums(famsizes_prim[0:(0.05*nrow(famsizes_prim)),])
tenperc_resp <- colSums(famsizes_prim[0:(0.1*nrow(famsizes_prim)),])
twentyperc_resp <- colSums(famsizes_prim[0:(0.2*nrow(famsizes_prim)),])
fourtyperc_resp <- colSums(famsizes_prim[0:(0.4*nrow(famsizes_prim)),])
hundredperc_resp <- colSums(famsizes_prim)

cum_size_top_fams <- data.frame("2" = twoperc_resp,
                                "5" = fiveperc_resp, 
                                "10" = tenperc_resp,
                                "20" = twentyperc_resp, 
                                "40" = fourtyperc_resp,
                                "100" = hundredperc_resp)

cum_size_plot_data <- data.frame("families" = as.factor(c(2, 5, 10, 20, 40, 100)),
                                 "response" = colMeans(cum_size_top_fams),
                                 "sd" = sapply(cum_size_top_fams, sd))

barplot_cumsize_fam <- ggplot(cum_size_plot_data, aes(x = families, y = response)) +
  geom_col(fill = "#9f2a63", width = 0.5, col = "black") + 
  geom_errorbar(aes(ymin = response - sd, ymax = response + sd), 
                position = "dodge", width = 0.25, linewidth = 0.8) +
  theme_clean() + th +
  scale_y_continuous(limits = c(0, 100),
                     n.breaks = 6,
                     expand = c(0, 0)) +
  labs(x = "Familes (%), ordered by size",
       y = "Response (%)")


panel_famsize_disparity <- plot_grid(plotlist = list(largest_fam_plot, barplot_cumsize_fam),
                                     nrow = 1,
                                     rel_widths = c(0.4, 0.6),
                                     labels = c("A", "B"))

ggsave(paste0(folder, run_name, "_4_panel_famsize_disparity.jpg"), plot = panel_famsize_disparity, 
       width = 7, height = 3, units = "in")

sink(file = paste0(folder, run_name, "_4_famsize_disparity.txt"))

print("PANEL A")
print("Mean largest family size:")
mean(largest_fam$value[largest_fam$metric == "largest_fam"])
print("SD:")
sd(largest_fam$value[largest_fam$metric == "largest_fam"])

print("Mean median family size:")
mean(largest_fam$value[largest_fam$metric == "median_fam"])
print("SD:")
sd(largest_fam$value[largest_fam$metric == "median_fam"])

print("PANEL B")
cum_size_plot_data

sink(file = NULL)

#-------------------------------------------------------------------------------
#Figure 5: panel with correlations across responses
#-------------------------------------------------------------------------------
# Panel A
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
  #geom_line(size = 0.05) +
  geom_line(data = total_resp_mean, aes(x = time, y = mean_log_cells), 
            color = "black", size = 0.5, inherit.aes = F) +  # Add mean line
  geom_ribbon(data = total_resp_mean, aes(x = time, 
                                          ymin = mean_log_cells - sd_log_cells, 
                                          ymax = mean_log_cells + sd_log_cells), 
              fill = "#9f2a63", alpha = 0.5, inherit.aes = F) +
  #  geom_ribbon(data = total_resp_mean, aes(x = time, ymin = log10(mean_cells)))
  labs(x = "Time (days)",
       y = "Nr. of cells (log-scale)") + 
  theme_clean() + 
  th

# Panel B
max_fam_stats_long <- convert_cols_to_numeric(max_fam_stats_long)

max_fam_boxplot_cor_data <- max_fam_stats_long %>%
  select(cor_prim_sec, cor_sec_ter) %>%
  pivot_longer(cols = all_of(c("cor_prim_sec", "cor_sec_ter")),
               names_to = "metric",
               values_to = "value")

boxplot_correlations <- ggplot(max_fam_boxplot_cor_data, aes(x = metric, y = value)) + 
  geom_boxplot() +
  geom_signif(comparisons = list(c(1, 2)), map_signif_level = T) +
  theme_clean() + th +
  labs(
    y = "Spearman r",
    x = element_blank()) +
  scale_x_discrete(labels = c("1vs2", "2vs3")) +
  ylim(-0.05, 0.6) + 
  theme(plot.margin = margin(t = 0.5, r = 0.5, l = 0.5, b = 0.7, unit = "cm"))

# Panel C and D
cor_prim_sec_nrQ_plot <- plot_prim_sec_max_fam(max_fam_table_single, label_burst_divs = "shape")
cor_sec_ter_nrQ_plot <-  plot_sec_ter_max_fam(max_fam_table_single, label_burst_divs = "shape")
legend_cor_nrQ_plot <- get_legend(cor_prim_sec_nrQ_plot)
cor_prim_sec_nrQ_plot <- cor_prim_sec_nrQ_plot + theme(legend.position = "none")
cor_sec_ter_nrQ_plot <- cor_sec_ter_nrQ_plot + theme(legend.position = "none")

gt <- arrangeGrob(grobs = list(legend_cor_nrQ_plot, 
                               response_plot_multi,
                               cor_prim_sec_nrQ_plot, 
                               cor_sec_ter_nrQ_plot,
                               boxplot_correlations), 
                  nrow = 8, 
                  layout_matrix = rbind(c(2,2,2,5), 
                                        c(2,2,2,5),
                                        c(2,2,2,5),  
                                        c(3,3,4,4), 
                                        c(3,3,4,4), 
                                        c(3,3,4,4), 
                                        c(3,3,4,4), 
                                        c(1,1,1,1)),
                  padding = unit(0, "line"))

panel_plot <- as_ggplot(gt) +             # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C", "D"), 
                  size = 13,
                  x = c(0, 0.75, 0, 0.5), 
                  y = c(1, 1, 0.63, 0.63))

ggsave(paste0(folder, run_name, "_5_panel_plot_cor.jpg"), plot = panel_plot, 
       width = 6.5, height = 6.5, units = "in")

#-------------------------------------------------------------------------------
#Table 3: response characteristics
#-------------------------------------------------------------------------------
processed_total_resp_stats <- data.frame(rbind(get_processed_stats(prim_response_stats),
                                               get_processed_stats(sec_response_stats),
                                               get_processed_stats(ter_response_stats)))

sink(file = paste0(folder, run_name, "_T3_response_characteristics.txt"))

print("Table 3: characteristics of the response")
print(processed_total_resp_stats)

sink(file = NULL)

#-------------------------------------------------------------------------------
#Fig. 6: boxplot correlation number of Q cells primary response and family size across responses
#-------------------------------------------------------------------------------
Q_famsize_boxplot <- plot_Q_famsize_boxplots(max_fam_table_single)
Q_famsize_stats <- get_Q_famsize_stats(max_fam_table_single)

ggsave(paste0(folder, run_name, "_6_Q_famsize_boxplot.png"), plot = Q_famsize_boxplot, 
       width = 1950, height = 3093, units = "px")

sink(file = paste0(folder, run_name, "_6_Q_famsize_boxplot.txt"))

print("Correlation between nr. Q cells (primary resp.) and family size")
print(Q_famsize_stats)

sink(file = NULL)

