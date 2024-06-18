# this will iteratively do multiple runs and calculate statistics

nr_sims <- 2
df_famsizes <- list()
total_resp <- list()
max_fam_table <- list()

#-------------------------------------------------------------------------------

folder <- "./results/model1/"
run_name <- "st3_multiple"
seed <- 4321
families = 2

bp_rule <- 'runif(1, min = 0, max = 3.5)'
dp_rule <- 0.5 #default
rq_rule <- 0.5 #default
t_start_dist <- 'rlnorm(nr_of_families, meanlog = 1.4, sdlog = 0.3)'
t_run_rule <- 'runif(1, min = 0.6, max = 4)' 
# was 0.5 in st2_nomax, but 0.6 ensures that all families can at least to their 4 burst 
# divisions, because 0.15*4 = 0.6
nr_burst_divs <- 'sample(c(2,3,4), nr_of_families, replace = TRUE)'
quality_dist <- NULL
ASD <- FALSE
burst_time <- 0.15 #default
max_run_time <- NULL
min_t_start <- 2.5

#-------------------------------------------------------------------------------

for (sim in seq(nr_sims)){
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

  #-----------------------------------------------------------------------------
  # distribution family sizes primary response
  
  df_famsizes[[sim]] <- generate_famsize_table_multidays(prim_parameters)
  
  #-----------------------------------------------------------------------------
  # total response table
  
  prim_resp <- generate_total_response_table(prim_parameters, 0, 50, 0.1)
  
  sec_resp <- generate_total_response_table(sec_parameters, 0, 50, 0.1)
  sec_resp$time <- sec_resp$time + 50 # secondary response starts at day 50
  
  ter_resp <- generate_total_response_table(ter_parameters, 0, 50, 0.1)
  ter_resp$time <- ter_resp$time + 100 # tertiary response starts at day 100
  
  total_resp[[sim]] <- rbind(prim_resp, sec_resp, ter_resp)
  
  #-----------------------------------------------------------------------------
  # max fam table
  
  max_fam_table[[sim]] <- generate_max_fam_table(prim_parameters, 
                                          sec_parameters, 
                                          ter_parameters, 
                                          type = "time", 
                                          timepoint = 7)
  
  #-----------------------------------------------------------------------------
  print(sim)
  
}

