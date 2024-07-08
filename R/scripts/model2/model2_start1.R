folder <- "./results/model2/"
run_name <- "start1"
seed <- 11
families = 500
quality_dist <- 'plnorm(t_start, meanlog=1.4, sdlog=0.3, lower.tail=F)'

nr_burst_divs <- 'ceiling(quality*3) + 1'
bp_rule <- 'runif(1, min = 0, max = 3.5)'
t_run_rule <- 'runif(1, min = 0.5, max = 4)'

t_start_dist <- 'rlnorm(nr_of_families, meanlog = 1.4, sdlog = 0.3)'
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
                                   min_t_start = min_t_start,
                                   quality_dist = quality_dist)

famsize_dist <- plot_grid_famsize_dist(prim_parameters)
ggsave(paste0(folder, run_name, "_famsize_dist.png"), plot = famsize_dist[[1]])

Q_famsize <- plot_Q_famsize(prim_parameters)
ggsave(paste0(folder, run_name, "_Q_famsize.png"), plot = Q_famsize[[1]])

#prim_resp <- plot_response(prim_parameters, 0, 50, 0.1)
#ggsave(paste0(folder, run_name, "_response.png"), plot = prim_resp[[2]])

sec_parameters <- pick_parameters(response_nr = 2,
                                  prev_parameters = prim_parameters,
                                  bp_rule = bp_rule,
                                  rq_rule = rq_rule,
                                  t_run_rule = t_run_rule,
                                  t_start_dist = t_start_dist,
                                  nr_burst_divs = nr_burst_divs,
                                  min_t_start = min_t_start,
                                  quality_dist = quality_dist)

ter_parameters <- pick_parameters(response_nr = 3,
                                  prev_parameters = sec_parameters,
                                  bp_rule = bp_rule,
                                  rq_rule = rq_rule,
                                  t_run_rule = t_run_rule,
                                  t_start_dist = t_start_dist,
                                  nr_burst_divs = nr_burst_divs,
                                  min_t_start = min_t_start,
                                  quality_dist = quality_dist)

max_fam_plot <- generate_max_fam_plot(prim_parameters = prim_parameters, 
                                      sec_parameters = sec_parameters, 
                                      ter_parameters = ter_parameters,
                                      type = "time",
                                      time = 7,
                                      relative = F,
                                      show_title = T)

ggsave(paste0(folder, run_name, "_max_fam_plot.png"), plot = max_fam_plot[[2]])

sink(file = paste0(folder, run_name, "_output.txt"))

print(paste("Seed:", seed))
print(paste("Nr. of families:", families))
print(paste("Distribution start times:", t_start_dist, "\nWith minimum:", min_t_start))
print(paste("Nr of burst divisons rule:", nr_burst_divs))
print(paste("Distribution proliferation rates:", bp_rule))
print(paste("Distribution of stop times:", t_run_rule))
print(paste("Chance to become Q cell:", rq_rule))
print(paste("Qualtiy distribution:", quality_dist))

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
