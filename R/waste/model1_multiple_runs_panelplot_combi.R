
if(Sys.info()[[4]]=="LAPTOP-3RJSLMKV") {
  print("Working from R-project")
  setwd("~/BiBC/Major_internship/tcell-family-based")
} else {setwd('~/Documents')
  .libPaths("./R")}

source("./R/functions_multiple_sims.R") 

# Define parameters
nr_sims <- 100
folder <- "./results/model1/iterative_runs/model1/"
run_name <- "m1_multiple_"
seed <- 4321
families <- 500

set.seed(seed)

bp_rule <- 'runif(1, min = 0, max = 3.5)'
dp_rule <- 0.5 # default
rq_rule <- 0.5 # default
t_start_dist <- 'rlnorm(nr_of_families, meanlog = 1.4, sdlog = 0.3)'
t_run_rule <- 'runif(1, min = 0.6, max = 4)' 
nr_burst_divs <- 'sample(c(2,3,4), nr_of_families, replace = TRUE)'
quality_dist <- NULL
ASD <- FALSE
burst_time <- 0.15 # default
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

sec_parameters <- pick_parameters(response_nr = 2,
                                  prev_parameters = prim_parameters,
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

ter_parameters <- pick_parameters(response_nr = 3,
                                  prev_parameters = sec_parameters,
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

#-------------------------------------------------------------------------------

Q_famsize_table <- generate_Q_famsize_table(prim_parameters)
Q_famsize_col_plot <- plot_Q_famsize(Q_famsize_table, label_burst_divs = "col", show_legend = T, linear_model = F)

max_fam_table <- generate_max_fam_table(prim_parameters, 
                                        sec_parameters, 
                                        ter_parameters, 
                                        type = "time", 
                                        timepoint = 7)
cor_prim_sec_nrQ_plot <- plot_prim_sec_max_fam(max_fam_table, label_burst_divs = "shape")
cor_sec_ter_nrQ_plot <-  plot_sec_ter_max_fam(max_fam_table, label_burst_divs = "shape")
legend_cor_nrQ_plot <- get_legend(cor_prim_sec_nrQ_plot)
cor_prim_sec_nrQ_plot <- cor_prim_sec_nrQ_plot + theme(legend.position = "none")
cor_sec_ter_nrQ_plot <- cor_sec_ter_nrQ_plot + theme(legend.position = "none")

#-------------------------------------------------------------------------------

gt_combi <- arrangeGrob(grobs = list(legend_cor_nrQ_plot, 
                                     response_plot_multi, 
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

ggsave(paste0(folder, run_name, "_panel_plot_combi.pdf"), plot = panel_plots_combi, 
       width = 6.5, height = 6.5, units = "in")

ggsave(paste0(folder, run_name, "_panel_plot_combi.jpg"), plot = panel_plots_combi, 
       width = 6.5, height = 6.5, units = "in")
