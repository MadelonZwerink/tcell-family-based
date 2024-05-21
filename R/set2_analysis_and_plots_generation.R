folder <- "./results/"
run_name <- "set2new"
seed <- 11
families = 500
bp_rule <- 'runif(1, min = 0, max = 3.5)'
t_run_rule <- 'runif(1, min = 0.5, max = 4)'
t_start_dist <- 'rlnorm(nr_of_families, meanlog = 1.4, sdlog = 0.3)'
nr_burst_divs <- 'sample(c(2,3,4), nr_of_families, replace = TRUE)'
rq_rule <- 0.5
min_t_start <- 2.5

set.seed(seed)

prim_parameters <- pick_parameters(nr_of_families = families,
                                   response_nr = 1,
                                   bp_rule = bp_rule,
                                   rq_rule = rq_rule,
                                   t_run_rule = t_run_rule,
                                   t_start_dist = t_start_dist,
                                   nr_burst_divs = nr_burst_divs,
                                   min_t_start = min_t_start)

famsize_dist <- plot_grid_famsize_dist(prim_parameters)
ggsave(paste0(folder, run_name, "_famsize_dist.pdf"), plot = famsize_dist[[1]], 
       width = 2500, height = 1800, units = "px")

Q_famsize <- plot_Q_famsize(prim_parameters)
ggsave(paste0(folder, run_name, "_Q_famsize.pdf"), plot = Q_famsize[[1]], 
       width = 2500, height = 1800, units = "px")

prim_resp <- plot_response(prim_parameters, 0, 50, 0.1)
ggsave(paste0(folder, run_name, "_response.pdf"), plot = prim_resp[[2]], 
       width = 2500, height = 1800, units = "px")

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

max_fam_plot <- generate_max_fam_plot(prim_parameters = prim_parameters, 
                                      sec_parameters = sec_parameters, 
                                      ter_parameters = ter_parameters,
                                      type = "time",
                                      time = 7,
                                      relative = F,
                                      show_title = T)

ggsave(paste0(folder, run_name, "_max_fam_plot.pdf"), plot = max_fam_plot[[2]], 
       width = 2500, height = 1300, units = "px")

# Improve this to make a panel for report!
grid.arrange(grid.arrange(prim_resp[[2]], Q_famsize[[1]], nrow = 1), max_fam_plot[[2]])

sink(file = paste0(folder, run_name, "_output.txt"))

print(paste("Seed:", seed))
print(paste("Nr. of families:", families))
print(paste("Distribution start times:", t_start_dist, "\nWith minimum:", min_t_start))
print(paste("Nr of burst divisons rule:", nr_burst_divs))
print(paste("Distribution proliferation rates:", bp_rule))
print(paste("Distribution of stop times:", t_run_rule))
print(paste("Chance to become Q cell:", rq_rule))

print("Famsize stats")
famsize_dist[[2]]

print("Correlation nr of Q cells and family size (primary response)")
Q_famsize[[2]]

print("Family sizes on day 7 per response")
max_fam_plot[[3]]

print("Correlation family size prim and sec")
max_fam_plot[[4]]

print("Correlation family size sec and ter")
max_fam_plot[[5]]

sink(file = NULL)

#-------------------------------------------------------------------------------


