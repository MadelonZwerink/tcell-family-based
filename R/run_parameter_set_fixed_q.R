
get_previous_q_per_fam <- function(parameters){
  qualities <- unique(parameters[,c("fam_nr", "quality")])

  q_cells <- parameters %>% 
    filter(cell_type == "Q") %>%
    group_by(fam_nr) %>%
    summarise(q_cells = n()) %>%
    right_join(qualities, by = "fam_nr")
  
  q_cells[is.na(q_cells)] <- 0
  
  return(q_cells)
}

generate_fixed_q_vector <- function(parameters){
  q_cells <- get_previous_q_per_fam(parameters)
  
  q_vector <- rep(q_cells$quality, q_cells$q_cells)
  
  return(q_vector)
}

#-------------------------------------------------------------------------------
source("./R/speed_discrete_burst.R") 

set.seed(seed)

print("Generating parameters...")

if(exists("quality_noise") == F){
  quality_noise <- FALSE
}
if(exists("resolution_resp_plot") == F){
  resolution_resp_plot <- 0.1
}
if(exists("uniform_fam") == F){
  uniform_fam <- FALSE
}
if(exists("q_noise_dist") == F){
  q_noise_dist <- NULL
}

prim_parameters <- pick_parameters(bp_rule = bp_rule,
                                   dp_rule = dp_rule,
                                   rq_rule = rq_rule,
                                   t_start_dist = t_start_dist,
                                   t_run_rule = t_run_rule,
                                   nr_of_families = nr_of_families,
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
                                  quality_dist = quality_dist_2,
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
                                  quality_dist = quality_dist_2,
                                  quality_noise = quality_noise,
                                  q_noise_dist = q_noise_dist,
                                  uniform_fam = uniform_fam,
                                  ASD = ASD,
                                  burst_time = burst_time,
                                  max_run_time = max_run_time,
                                  min_t_start = min_t_start)

#-------------------------------------------------------------------------------

print("Calculating family sizes and other data...")

df_famsizes <- generate_famsize_table_multidays(prim_parameters)
#df_famsizes_sec  <- generate_famsize_table_multidays(sec_parameters)
#df_famsizes_ter <- generate_famsize_table_multidays(ter_parameters)

Q_famsize_table <- generate_Q_famsize_table(prim_parameters)

max_fam_table <- generate_max_fam_table(prim_parameters, 
                                        sec_parameters, 
                                        ter_parameters, 
                                        type = "time", 
                                        timepoint = 7)

#-------------------------------------------------------------------------------

prim_resp <- generate_total_response_table(prim_parameters, 0, 50, resolution_resp_plot)

sec_resp <- generate_total_response_table(sec_parameters, 0, 50, resolution_resp_plot)
sec_resp$time <- sec_resp$time + 50 # secondary response starts at day 50

ter_resp <- generate_total_response_table(ter_parameters, 0, 50, resolution_resp_plot)
ter_resp$time <- ter_resp$time + 100 # tertiary response starts at day 100

total_resp <- rbind(prim_resp, sec_resp, ter_resp)

#-------------------------------------------------------------------------------

print("Making plots...")

response_plot <- plot_response(total_resp)

famsize_dist_plot <- plot_grid_famsize_dist(df_famsizes)

#sec_famsize_dist_plot <- plot_v_grid_famsize_dist(df_famsizes_sec)
#ter_famsize_dist_plot <- plot_v_grid_famsize_dist(df_famsizes_ter)

Q_famsize_shape_plot <- plot_Q_famsize(Q_famsize_table, label_burst_divs = "shape")
Q_famsize_col_plot <- plot_Q_famsize(Q_famsize_table, label_burst_divs = "col", show_legend = T, linear_model = F)

Q_famsize_all_resp_plot <- plot_Q_famsize_boxplots(max_fam_table)

cor_prim_sec_nrQ_plot <- plot_prim_sec_max_fam(max_fam_table, label_burst_divs = "shape")
cor_sec_ter_nrQ_plot <-  plot_sec_ter_max_fam(max_fam_table, label_burst_divs = "shape")

cor_prim_sec_col_plot <- plot_prim_sec_max_fam(max_fam_table, label_burst_divs = "col")
cor_sec_ter_col_plot <-  plot_sec_ter_max_fam(max_fam_table, label_burst_divs = "col")

#-------------------------------------------------------------------------------

legend_cor_col_plot <- get_legend(cor_prim_sec_col_plot)

cor_prim_sec_col_plot <- cor_prim_sec_col_plot + theme(legend.position = "none")
cor_sec_ter_col_plot <- cor_sec_ter_col_plot + theme(legend.position = "none")
Q_famsize_col_nolegend_plot <- Q_famsize_col_plot + theme(legend.position = "none")

gt_col <- arrangeGrob(grobs = list(legend_cor_col_plot, 
                                   response_plot, 
                                   Q_famsize_col_nolegend_plot, 
                                   cor_prim_sec_col_plot, 
                                   cor_sec_ter_col_plot), 
                      nrow = 8, layout_matrix = rbind(c(2,2,3,3), 
                                                      c(2,2,3,3),
                                                      c(2,2,3,3),  
                                                      c(4,4,5,5), 
                                                      c(4,4,5,5), 
                                                      c(4,4,5,5), 
                                                      c(4,4,5,5), 
                                                      c(1,1,1,1)),
                      padding = unit(0, "line"))

panel_plots_col <- as_ggplot(gt_col) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C", "D"), size = 13,
                  x = c(0, 0.5, 0, 0.5), y = c(1, 1, 0.63, 0.63))

#-------------------------------------------------------------------------------

legend_cor_nrQ_plot <- get_legend(cor_prim_sec_nrQ_plot)

cor_prim_sec_nrQ_plot <- cor_prim_sec_nrQ_plot + theme(legend.position = "none")
cor_sec_ter_nrQ_plot <- cor_sec_ter_nrQ_plot + theme(legend.position = "none")

gt_nrQ <- arrangeGrob(grobs = list(legend_cor_nrQ_plot, 
                                   response_plot, 
                                   Q_famsize_shape_plot, 
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

panel_plots_nrQ <- as_ggplot(gt_nrQ) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C", "D"), size = 13,
                  x = c(0, 0.5, 0, 0.5), y = c(1, 1, 0.63, 0.63))

#-------------------------------------------------------------------------------

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
dd_P_burst <- dd_total_burst[which(prim_parameters$cell_type == "P")]

plot_dd_total_burst <- ggplot(as.data.frame(dd_total_burst), aes(x = dd_total_burst))  + 
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, fill = "#9f2a63", col = "black") + 
  theme_clean() + th +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Division destiny (burst included)",
       y = "Density",
       title = "All branches")

plot_dd_P_burst <- ggplot(as.data.frame(dd_P_burst), aes(x = dd_P_burst)) + 
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, fill = "#9f2a63", col = "black") + 
  theme_clean() + th +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Division destiny (burst included)",
       y = "Density",
       title = "Proliferating branches")

panel_dd_inc_burst <- plot_grid(plotlist = list(plot_dd_total_burst, plot_dd_P_burst),
                                labels = c("A", "B"))  


dd_total <- prim_parameters$bp * (prim_parameters$t_run - prim_parameters$nr_burst_divs * prim_parameters$t_burst)
dd_P <- dd_total[which(prim_parameters$cell_type == "P")]

plot_dd_total <- ggplot(as.data.frame(dd_total), aes(x = dd_total))  + 
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, fill = "#9f2a63", col = "black") + 
  theme_clean() + th +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Division destiny",
       y = "Density",
       title = "All branches")

plot_dd_P <- ggplot(as.data.frame(dd_P), aes(x = dd_P)) + 
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, fill = "#9f2a63", col = "black") + 
  theme_clean() + th +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Division destiny",
       y = "Density",
       title = "Proliferating branches")

panel_dd <- plot_grid(plotlist = list(plot_dd_total, plot_dd_P),
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

ggsave(paste0(folder, run_name, "_panel_plot_nrQ.pdf"), plot = panel_plots_nrQ, 
       width = 6.5, height = 6.5, units = "in")

ggsave(paste0(folder, run_name, "_panel_plot_col.pdf"), plot = panel_plots_col, 
       width = 6.5, height = 6.5, units = "in")

ggsave(paste0(folder, run_name, "_panel_plot_combi.pdf"), plot = panel_plots_combi, 
       width = 6.5, height = 6.5, units = "in")

ggsave(paste0(folder, run_name, "_panel_plot_combi.jpg"), plot = panel_plots_combi, 
       width = 6.5, height = 6.5, units = "in")

ggsave(paste0(folder, run_name, "_famsize_dist_plot.pdf"), plot = famsize_dist_plot, 
       width = 5, height = 8, units = "in")

ggsave(paste0(folder, run_name, "_famsize_dist_plot.jpg"), plot = famsize_dist_plot, 
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

ggsave(paste0(folder, run_name, "_largest_fam_plot.jpg"), plot = largest_fam_plot, 
       width = 4, height = 6, units = "in")

ggsave(paste0(folder, run_name, "_cumsize_fam_barplot.jpg"), plot = barplot_cumsize_fam, 
       width = 4, height = 3, units = "in")

ggsave(paste0(folder, run_name, "_panel_famsize_dist.jpg"), plot = panel_famsize_dist_gerlach, 
       width = 7, height = 3, units = "in")

ggsave(paste0(folder, run_name, "_panel_dd.jpg"), plot = panel_dd, 
       width = 7, height = 3, units = "in")

ggsave(paste0(folder, run_name, "_panel_dd_inc_burst.jpg"), plot = panel_dd_inc_burst, 
       width = 7, height = 3, units = "in")

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

