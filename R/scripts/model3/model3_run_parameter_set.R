#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------

source("./R/speed_discrete_burst.R") 
source("./R/model3_functions.R")

set.seed(seed)

print("Generating parameters...")

#-------------------------------------------------------------------------------

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

div_index_prim <- generate_div_index_df(prim_parameters)

recruited_cells_prim_parameters <- generate_recruited_cells(div_index_prim, 
                                                            recruitment_mean,
                                                            recruitment_sd)
recruited_cells_prim <- generate_processed_recruited_cells(recruited_cells_prim_parameters)

#-------------------------------------------------------------------------------

sec_parameters <- pick_parameters(response_nr = 2,
                                  prev_parameters = recruited_cells_prim,
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

div_index_sec <- generate_div_index_df(sec_parameters)

recruited_cells_sec_parameters <- generate_recruited_cells(div_index_sec, 
                                                           recruitment_mean,
                                                           recruitment_sd)
recruited_cells_sec <- generate_processed_recruited_cells(recruited_cells_sec_parameters)


#-------------------------------------------------------------------------------

ter_parameters <- pick_parameters(response_nr = 3,
                                  prev_parameters = recruited_cells_sec,
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

div_index_ter <- generate_div_index_df(ter_parameters)

recruited_cells_ter_parameters <- generate_recruited_cells(div_index_ter, 
                                                           recruitment_mean,
                                                           recruitment_sd)
recruited_cells_ter <- generate_processed_recruited_cells(recruited_cells_ter_parameters)


#-------------------------------------------------------------------------------

print("Calculating family sizes and other data...")

df_famsizes <- generate_famsize_table_multidays(prim_parameters)

Q_famsize_table <- generate_recruited_famsize_table(prim_parameters, recruited_cells_prim)

max_fam_table <- generate_max_fam_table(prim_parameters, 
                                        sec_parameters, 
                                        ter_parameters, 
                                        type = "time", 
                                        timepoint = 7,
                                        model = "noQ",
                                        recruited_cells_prim = recruited_cells_prim,
                                        recruited_cells_sec = recruited_cells_sec,
                                        recruited_cells_ter = recruited_cells_ter)

#-------------------------------------------------------------------------------

prim_resp <- generate_total_response_table(prim_parameters, 0, 30, 2,
                                           model = "noQ", recruited_cells = recruited_cells_prim_parameters)

sec_resp <- generate_total_response_table(sec_parameters, 0, 30, 2,
                                          model = "noQ", recruited_cells = recruited_cells_sec_parameters)
sec_resp$time <- sec_resp$time + 50 # secondary response starts at day 50

ter_resp <- generate_total_response_table(ter_parameters, 0, 30, 2,
                                          model = "noQ", recruited_cells = recruited_cells_ter_parameters)
ter_resp$time <- ter_resp$time + 100 # tertiary response starts at day 100

total_resp <- rbind(prim_resp, sec_resp, ter_resp)

#-------------------------------------------------------------------------------

print("Making plots...")

response_plot <- plot_response(total_resp)

famsize_dist_plot <- plot_grid_famsize_dist(famsizes_table = df_famsizes)

Q_famsize_col_plot <- plot_Q_famsize(Q_famsize_table, label_burst_divs = "col", show_legend = T, linear_model = F)

Q_famsize_all_resp_plot <- plot_Q_famsize_boxplots(max_fam_table, scale_Q_max = max(max_fam_table$Q_cells))

cor_prim_sec_nrQ_plot <- plot_prim_sec_max_fam(max_fam_table, label_burst_divs = "shape", scale_Q_max = max(max_fam_table$Q_cells))
cor_sec_ter_nrQ_plot <-  plot_sec_ter_max_fam(max_fam_table, label_burst_divs = "shape", scale_Q_max = max(max_fam_table$Q_cells))

#-------------------------------------------------------------------------------

legend_cor_nrQ_plot <- get_legend(cor_prim_sec_nrQ_plot)

cor_prim_sec_nrQ_plot <- cor_prim_sec_nrQ_plot + theme(legend.position = "none")
cor_sec_ter_nrQ_plot <- cor_sec_ter_nrQ_plot + theme(legend.position = "none")

gt_combi <- arrangeGrob(grobs = list(legend_cor_nrQ_plot, 
                                     response_plot, 
                                     Q_famsize_col_plot, 
                                     cor_prim_sec_nrQ_plot, 
                                     cor_sec_ter_nrQ_plot), 
                        nrow = 8, layout_matrix = rbind(c(2,2,3,3), 
                                                        c(2,2,3,3),
                                                        c(2,2,3,3),  
                                                        c(4,4,5,5), 
                                                        c(4,4,5,5), 
                                                        c(4,4,5,5), 
                                                        c(4,4,5,5), 
                                                        c(1,1,1,1)),
                        padding = unit(0, "line"))

panel_plots_combi <- as_ggplot(gt_combi) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C", "D"), size = 13,
                  x = c(0, 0.5, 0, 0.5), y = c(1, 1, 0.63, 0.63))

#-------------------------------------------------------------------------------

largest_fam <- get_largest_fam(max_fam_table)
largest_fam <- data.frame("metric" = names(largest_fam),
                          "value" = largest_fam)

largest_fam_plot <- ggplot(largest_fam, aes(y = log10(value), x = metric)) + 
  geom_col(fill = "#9f2a63", width = 0.5, col = "black") +
  ylim(0,6) +
  theme_clean() + th +
  labs(x = element_blank(),
       y = "Family size (cell number, 10log-scale)")  +
  scale_x_discrete(labels = c("Largest family", "Median family")) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.margin = margin(l = 1, b = 0.8, t = 0.6, unit = "cm"))

famsizes_prim <- get_famsize_fraction(max_fam_table)

cum_size_top_fams <- c(sum(famsizes_prim[0:(0.02*length(famsizes_prim))]),
                       sum(famsizes_prim[0:(0.05*length(famsizes_prim))]), 
                       sum(famsizes_prim[0:(0.1*length(famsizes_prim))]),
                       sum(famsizes_prim[0:(0.2*length(famsizes_prim))]), 
                       sum(famsizes_prim[0:(0.4*length(famsizes_prim))]),
                       sum(famsizes_prim[0:(1*length(famsizes_prim))]))

cum_size_plot_data <- data.frame("families" = as.factor(c(2, 5, 10, 20, 40, 100)),
                                 "response" = cum_size_top_fams)

barplot_cumsize_fam <- ggplot(cum_size_plot_data, aes(x = families, y = response)) +
  geom_col(fill = "#9f2a63", width = 0.5, col = "black") +
  theme_clean() + th +
  scale_y_continuous(limits = c(0, 100),
                     n.breaks = 6,
                     expand = c(0, 0)) +
  labs(x = "Familes (%), ordered by size",
       y = "Response (%)")

panel_famsize_dist_gerlach <- plot_grid(plotlist = list(largest_fam_plot, barplot_cumsize_fam),
                                        nrow = 1,
                                        rel_widths = c(0.4, 0.6),
                                        labels = c("A", "B"))

#-------------------------------------------------------------------------------

dd_total_burst <- (prim_parameters$bp * (prim_parameters$t_run - prim_parameters$nr_burst_divs * prim_parameters$t_burst)) + prim_parameters$nr_burst_divs

plot_dd_total_burst <- ggplot(as.data.frame(dd_total_burst), aes(x = dd_total_burst))  + 
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, fill = "#9f2a63", col = "black") + 
  theme_clean() + th +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Division destiny (burst included)",
       y = "Density",
       title = "All branches")

dd_total <- prim_parameters$bp * (prim_parameters$t_run - prim_parameters$nr_burst_divs * prim_parameters$t_burst)

plot_dd_total <- ggplot(as.data.frame(dd_total), aes(x = dd_total))  + 
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, fill = "#9f2a63", col = "black") + 
  theme_clean() + th +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Division destiny",
       y = "Density",
       title = "All branches")

panel_dd <- plot_grid(plotlist = list(plot_dd_total_burst, plot_dd_total),
                      labels = c("A", "B"))  

#-------------------------------------------------------------------------------

print("Calculating statistics...")

prim_resp_stats <- get_total_response_stats(prim_resp)
sec_resp_stats <- get_total_response_stats(sec_resp)
ter_resp_stats <- get_total_response_stats(ter_resp)
response_stats <- rbind(prim_resp_stats, sec_resp_stats, ter_resp_stats)

Q_prim_famsize_stats <- get_Q_prim_famsize_stats(Q_famsize_table)

Q_famsize_stats <- get_Q_famsize_stats(max_fam_table)

famsize_dist_stats <- get_famsize_stats_multidays(df_famsizes)

max_fam_stats <- get_max_fam_stats(max_fam_table)

#-------------------------------------------------------------------------------

print("Saving plots and stats...")

ggsave(paste0(folder, run_name, "_panel_plot_combi.pdf"), plot = panel_plots_combi, 
       width = 6.5, height = 6.5, units = "in")

ggsave(paste0(folder, run_name, "_panel_plot_combi.png"), plot = panel_plots_combi, 
       width = 6.5, height = 6.5, units = "in")

ggsave(paste0(folder, run_name, "_famsize_dist_plot.png"), plot = famsize_dist_plot, 
       width = 5, height = 8, units = "in")

ggsave(paste0(folder, run_name, "_Q_famsize.pdf"), plot = Q_famsize_col_plot, 
       width = 2500, height = 1800, units = "px")

ggsave(paste0(folder, run_name, "_Q_famsize_boxplot.pdf"), plot = Q_famsize_all_resp_plot, 
       width = 1950, height = 3093, units = "px")

ggsave(paste0(folder, run_name, "_Q_famsize_boxplot.png"), plot = Q_famsize_all_resp_plot, 
       width = 1950, height = 3093, units = "px")

ggsave(paste0(folder, run_name, "_cor_prim_sec_plot.pdf"), plot = cor_prim_sec_nrQ_plot, 
       width = 1500, height = 1300, units = "px")

ggsave(paste0(folder, run_name, "_cor_sec_ter_plot.pdf"), plot = cor_sec_ter_nrQ_plot, 
       width = 1500, height = 1300, units = "px")

ggsave(paste0(folder, run_name, "_response.pdf"), plot = response_plot, 
       width = 2500, height = 1800, units = "px")

ggsave(paste0(folder, run_name, "_panel_dd.jpg"), plot = panel_dd, 
       width = 7, height = 3, units = "in")

ggsave(paste0(folder, run_name, "_panel_famsize_dist.jpg"), plot = panel_famsize_dist_gerlach, 
       width = 7, height = 3, units = "in")

#-------------------------------------------------------------------------------

if(!is.null(quality_dist)){  
  plot_parameters <- ggplot(prim_parameters[prim_parameters$cell_type == "P",]) +
    geom_point(aes(x = t_run, y = bp, col = quality)) +
    theme_clean() + th +
    labs(x = "Run time (days)",
         y = "Proliferation rate (divisions/day)") +
    guides(colour = guide_colourbar(show.limits = TRUE, 
                                    title.position = "top",
                                    barwidth = 10,
                                    barheight = 0.5)) +
    labs(color = "Quality") +
    scale_color_viridis_c(option = "inferno", direction = -1)
  
  if(exists("min_q") == T){
    plot_parameters <-  plot_parameters +     
      scale_color_viridis_c(option = "inferno", direction = -1,
                            limits = c(min_q, max_q),
                            labels = c(min_q, (min_q + max_q)/2, max_q),
                            breaks = c(min_q, (min_q + max_q)/2, max_q)) 
    
    panel_plots_quality <- generate_quality_effect_panelplot(min_q, med_q, max_q)
    
    ggsave(paste0(folder, run_name, "_panel_plot_quality.png"), plot = panel_plots_quality, 
           width = 2800, height = 1000, units = "px")
  }
  
  
  ggsave(paste0(folder, run_name, "_plot_parameters.png"), plot = plot_parameters, 
         width = 1500, height = 2000, units = "px")
} else {
  plot_parameters <- ggplot(prim_parameters[prim_parameters$cell_type == "P",]) +
    geom_point(aes(x = t_run, y = bp)) +
    theme_clean() + th +
    labs(x = "Run time (days)",
         y = "Proliferation rate (divisions/day)")
  
  ggsave(paste0(folder, run_name, "_plot_parameters.png"), plot = plot_parameters, 
         width = 1500, height = 2000, units = "px")
}

#-------------------------------------------------------------------------------

sink(file = paste0(folder, run_name, "_output.txt"))

print(paste("Run_name:", run_name))
print(paste("Seed:", seed))
print(paste("Nr. of families:", families))

print(paste("Distribution proliferation rates:", bp_rule))
print(paste("Death rate proliferation cells:", dp_rule))
print(paste("Chance to become Q cell:", rq_rule))
print(paste("Distribution start times:", t_start_dist))
print(paste("Minumum start time:", min_t_start))
print(paste("Distribution of run times:", t_run_rule))
print(paste("Nr of burst divisons rule:", nr_burst_divs))
print(paste("Quality distribution:", quality_dist))
print(paste("Quality noise:", quality_noise))
print(paste("Quality noise distribution:", q_noise_dist))
print(paste("Uniform family:", uniform_fam))
print(paste("ASD:", ASD))
print(paste("Time per burst division:", burst_time))
print(paste("Forced stop time:", max_run_time))

print("Famsize stats")
print(famsize_dist_stats)

print("Response stats")
print(response_stats)

print("Correlation nr of Q cells and family size (day 7)")
print(Q_famsize_stats)

print("Correlation nr of Q cells and maximum family size in primary response")
print(Q_prim_famsize_stats)

print("Family sizes on day 7 per response and correlation prim-sec and sec-ter family size")
print(max_fam_stats)

sink(file = NULL)

#-------------------------------------------------------------------------------

print(paste(run_name, "- run finished, results saved"))



