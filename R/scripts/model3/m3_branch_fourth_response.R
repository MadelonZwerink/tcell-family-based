
source("./R/speed_discrete_burst.R") 
source("./R/model3_functions.R")

folder <- "./results/model3/new_cluster/"
run_name <- "m3_DD_branch_2"

data_objects <- c("recruited_cells_prim",
                  "recruited_cells_sec",
                  "recruited_cells_ter",
                  "prim_parameters",
                  "sec_parameters",
                  "ter_parameters")
for(data_object in data_objects){
  assign(data_object, readRDS(paste0(folder, run_name, data_object, ".rds"))) 
}

bp_rule <- 'sqrt(q)'
dp_rule <- 0.5 #default
rq_rule <- 1 #only P cells will form
t_start_dist <- 'rlnorm(nr_of_families, meanlog = 1.4, sdlog = 0.3)'
t_run_rule <- 'sqrt(q) + (nr_burst_divs[i]*t_burst[i])' 
nr_burst_divs <- 'sample(c(2,3,4), nr_of_families, replace = TRUE)'
quality_dist <- "12.7 * rbeta(nr_of_families, shape1 = 0.227, shape2 = 1.625)"
quality_noise <- TRUE
q_noise_dist <- "12.7 * rbeta(cells[i], shape1 = 0.227, shape2 = 1.625)"
uniform_fam <- FALSE
ASD <- FALSE
burst_time <- 0.15 #default
max_run_time <- NULL
min_t_start <- 0

recruitment_mean <- 8
recruitment_sd <- 2

four_parameters <- pick_parameters(response_nr = 4,
                                  prev_parameters = recruited_cells_ter,
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

div_index_four <- generate_div_index_df(four_parameters)

recruited_cells_four_parameters <- generate_recruited_cells(div_index_four, 
                                                           recruitment_mean,
                                                           recruitment_sd)
recruited_cells_four <- generate_processed_recruited_cells(recruited_cells_four_parameters)
