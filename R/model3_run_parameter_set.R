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

prim_resp <- generate_total_response_table(prim_parameters, 0, 50, 1,
                                           model = "noQ", recruited_cells = recruited_cells_prim_parameters)

sec_resp <- generate_total_response_table(sec_parameters, 0, 50, 1,
                                          model = "noQ", recruited_cells = recruited_cells_sec_parameters)
sec_resp$time <- sec_resp$time + 50 # secondary response starts at day 50

ter_resp <- generate_total_response_table(ter_parameters, 0, 50, 1,
                                          model = "noQ", recruited_cells = recruited_cells_ter_parameters)
ter_resp$time <- ter_resp$time + 100 # tertiary response starts at day 100

total_resp <- rbind(prim_resp, sec_resp, ter_resp)

#-------------------------------------------------------------------------------

print("Making plots...")

response_plot <- plot_response(total_resp)

famsize_dist_v_plot <- plot_v_grid_famsize_dist(df_famsizes)

Q_famsize_col_plot <- plot_Q_famsize(Q_famsize_table, label_burst_divs = "col", show_legend = T, linear_model = F)

Q_famsize_all_resp_plot <- plot_Q_famsize_boxplots(max_fam_table)

cor_prim_sec_nrQ_plot <- plot_prim_sec_max_fam(max_fam_table, label_burst_divs = "shape")
cor_sec_ter_nrQ_plot <-  plot_sec_ter_max_fam(max_fam_table, label_burst_divs = "shape")

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

ggsave(paste0(folder, run_name, "_famsize_dist_v_plot.pdf"), plot = famsize_dist_v_plot, 
       width = 5, height = 8, units = "in")

ggsave(paste0(folder, run_name, "_famsize_dist_v_plot.png"), plot = famsize_dist_v_plot, 
       width = 5, height = 8, units = "in")

ggsave(paste0(folder, run_name, "_Q_famsize.pdf"), plot = Q_famsize_col_plot, 
       width = 2500, height = 1800, units = "px")

ggsave(paste0(folder, run_name, "_Q_famsize_boxplot.pdf"), plot = Q_famsize_all_resp_plot, 
       width = 1950, height = 3093, units = "px")

ggsave(paste0(folder, run_name, "_Q_famsize_boxplot.png"), plot = Q_famsize_all_resp_plot, 
       width = 1950, height = 3093, units = "px")

ggsave(paste0(folder, run_name, "_cor_prim_sec_plot.pdf"), plot = cor_prim_sec_plot, 
       width = 1500, height = 1300, units = "px")

ggsave(paste0(folder, run_name, "_cor_sec_ter_plot.pdf"), plot = cor_sec_ter_plot, 
       width = 1500, height = 1300, units = "px")

ggsave(paste0(folder, run_name, "_response.pdf"), plot = response_plot, 
       width = 2500, height = 1800, units = "px")

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



