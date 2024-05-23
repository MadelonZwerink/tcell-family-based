
#-------------------------------------------------------------------------------

set.seed(seed)

prim_parameters <- pick_parameters(nr_of_families = families,
                                   response_nr = 1,
                                   bp_rule = bp_rule,
                                   rq_rule = rq_rule,
                                   t_run_rule = t_run_rule,
                                   t_start_dist = t_start_dist,
                                   nr_burst_divs = nr_burst_divs,
                                   min_t_start = min_t_start)

sec_parameters <- pick_parameters(response_nr = 2,
                                  prev_parameters = prim_parameters,
                                  bp_rule = bp_rule,
                                  rq_rule = rq_rule,
                                  t_run_rule = t_run_rule,
                                  t_start_dist = t_start_dist,
                                  nr_burst_divs = nr_burst_divs,
                                  min_t_start = min_t_start)

ter_parameters <- pick_parameters(response_nr = 3,
                                  prev_parameters = sec_parameters,
                                  bp_rule = bp_rule,
                                  rq_rule = rq_rule,
                                  t_run_rule = t_run_rule,
                                  t_start_dist = t_start_dist,
                                  nr_burst_divs = nr_burst_divs,
                                  min_t_start = min_t_start)

#-------------------------------------------------------------------------------

df_famsizes <- generate_famsize_table_multidays(prim_parameters)

famsize_dist_h_plot <- plot_h_grid_famsize_dist(df_famsizes)

famsize_dist_sq_plot <- plot_sq_grid_famsize_dist(df_famsizes)

famsize_dist_stats <- get_famsize_stats_multidays(df_famsizes)

ggsave(paste0(folder, run_name, "_famsize_dist_h_plot.pdf"), plot = famsize_dist_h_plot, 
       width = 2500, height = 500, units = "px")

ggsave(paste0(folder, run_name, "_famsize_dist_sq_plot.pdf"), plot = famsize_dist_sq_plot, 
       width = 2500, height = 2200, units = "px")

#-------------------------------------------------------------------------------

Q_famsize_table <- generate_Q_famsize_table(prim_parameters)

Q_famsize_plot <- plot_Q_famsize(Q_famsize_table)

Q_famsize_stats <- get_Q_famsize_stats(Q_famsize_table)

ggsave(paste0(folder, run_name, "_Q_famsize.pdf"), plot = Q_famsize_plot, 
       width = 2500, height = 1800, units = "px")

#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------

max_fam_table <- generate_max_fam_table(prim_parameters, 
                       sec_parameters, 
                       ter_parameters, 
                       type = "time", 
                       timepoint = 7)

cor_prim_sec_plot <- plot_prim_sec_max_fam(max_fam_table)

cor_sec_ter_plot <-  plot_sec_ter_max_fam(max_fam_table)

max_fam_stats <- get_max_fam_stats(max_fam_table)

ggsave(paste0(folder, run_name, "_cor_prim_sec_plot.pdf"), plot = cor_prim_sec_plot, 
       width = 1500, height = 1300, units = "px")

ggsave(paste0(folder, run_name, "_cor_sec_ter_plot.pdf"), plot = cor_sec_ter_plot, 
       width = 1500, height = 1300, units = "px")

#-------------------------------------------------------------------------------

prim_resp <- generate_total_response_table(prim_parameters, 0, 50, 0.1)
prim_resp_stats <- get_total_response_stats(prim_resp)

sec_resp <- generate_total_response_table(sec_parameters, 0, 50, 0.1)
sec_resp$time <- sec_resp$time + 50 # secondary response starts at day 50
sec_resp_stats <- get_total_response_stats(sec_resp)

ter_resp <- generate_total_response_table(ter_parameters, 0, 50, 0.1)
ter_resp$time <- ter_resp$time + 100 # tertiary response starts at day 100
ter_resp_stats <- get_total_response_stats(ter_resp)

response_stats <- rbind(prim_resp_stats, sec_resp_stats, ter_resp_stats)
total_resp <- rbind(prim_resp, sec_resp, ter_resp)

response_plot <- plot_response(total_resp)

ggsave(paste0(folder, run_name, "_response.pdf"), plot = response_plot, 
       width = 2500, height = 1800, units = "px")

#-------------------------------------------------------------------------------

legend_cor_plot <- get_legend(cor_prim_sec_plot)

cor_prim_sec_plot <- cor_prim_sec_plot + theme(legend.position = "none")
cor_sec_ter_plot <- cor_sec_ter_plot + theme(legend.position = "none")

gt <- arrangeGrob(grobs = list(legend_cor_plot, response_plot, Q_famsize_plot, cor_prim_sec_plot, cor_sec_ter_plot), 
                  nrow = 8, layout_matrix = rbind(c(2,2,3,3), 
                                                  c(2,2,3,3),
                                                  c(2,2,3,3),  
                                                  c(4,4,5,5), 
                                                  c(4,4,5,5), 
                                                  c(4,4,5,5), 
                                                  c(4,4,5,5), 
                                                  c(1,1,1,1)),
                  padding = unit(0, "line"))

panel_plots <- as_ggplot(gt) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C", "D"), size = 13,
                  x = c(0, 0.5, 0, 0.5), y = c(1, 1, 0.63, 0.63))
panel_plots

#-------------------------------------------------------------------------------

sink(file = paste0(folder, run_name, "_output.txt"))

print(paste("Seed:", seed))
print(paste("Nr. of families:", families))
print(paste("Distribution start times:", t_start_dist, "\nWith minimum:", min_t_start))
print(paste("Nr of burst divisons rule:", nr_burst_divs))
print(paste("Distribution proliferation rates:", bp_rule))
print(paste("Distribution of stop times:", t_run_rule))
print(paste("Chance to become Q cell:", rq_rule))

print("Famsize stats")
famsize_dist_stats

print("Correlation nr of Q cells and family size (primary response)")
Q_famsize_stats

print("Family sizes on day 7 per response and correlation prim-sec and sec-ter family size")
max_fam_stats

sink(file = NULL)

#-------------------------------------------------------------------------------


