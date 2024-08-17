#!/usr/bin/env Rscript

if(Sys.info()[[4]]=="LAPTOP-3RJSLMKV") {
  print("Working from R-project")
  setwd("~/BiBC/Major_internship/tcell-family-based")
} else {setwd('~/Documents')
  .libPaths("./R")}

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

source("./R/functions_multiple_sims.R") 

# Function to generate parameters and responses
generate_parameters_multi <- function() {
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
                                     uniform_fam = uniform_fam,
                                     q_noise_dist = q_noise_dist,
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
                                    uniform_fam = uniform_fam,
                                    q_noise_dist = q_noise_dist,
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
                                    uniform_fam = uniform_fam,
                                    q_noise_dist = q_noise_dist,
                                    ASD = ASD,
                                    burst_time = burst_time,
                                    max_run_time = max_run_time,
                                    min_t_start = min_t_start)
  
  return(list(prim_parameters = prim_parameters,
              sec_parameters = sec_parameters,
              ter_parameters = ter_parameters))
}


# Define parameters
nr_sims <- 10
folder <- "./results/model4/"
run_name <- "m4_"
seed <- 4321
nr_of_families <- 500
families <- 500
set.seed(seed)

bp_rule <- 'sqrt(q*11.13)'
dp_rule <- 0.5 #default
rq_rule <- 0.5 #default
t_start_dist <- 'rlnorm(nr_of_families, meanlog = 1.4, sdlog = 0.3)'
t_run_rule <- 'sqrt(q*11.13) + (nr_burst_divs[i]*t_burst[i])'
nr_burst_divs <- 'sample(c(2,3,4), nr_of_families, replace = TRUE)'
quality_dist <- "rbeta(nr_of_families, shape1 = 2.165, shape2 = 3.295)"
quality_dist_2 <- 'generate_fixed_q_vector(prev_parameters)'
quality_noise <- FALSE
uniform_fam <- TRUE
ASD <- FALSE
burst_time <- 0.15 #default
max_run_time <- NULL
min_t_start <- 0
resolution_resp_plot <- 0.1

# Initialize lists
df_famsizes <- vector("list", nr_sims)
total_resp <- vector("list", nr_sims)
max_fam_table <- vector("list", nr_sims)
all_parameters_list <- vector("list", nr_sims)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Run simulations
lapply(seq_len(nr_sims), generate_response_multi)

# Calculate statistics
print("Calculating statistics...")

famsize_stats <- get_famsize_stats_multidays_multi(df_famsizes)
max_fam_stats <- sapply(max_fam_table, get_max_fam_stats, table = F)
max_fam_stats_long <- t(max_fam_stats)

prim_response_stats <- convert_response_stats_part(total_resp, max_time = 50)
sec_response_stats <- convert_response_stats_part(total_resp, max_time = 100)
ter_response_stats <- convert_response_stats_part(total_resp, max_time = 150)

#-------------------------------------------------------------------------------
print("Saving data...")

saveRDS(df_famsizes, file = paste0(folder, run_name, "df_famsizes.rds"))
saveRDS(total_resp, file = paste0(folder, run_name, "total_resp.rds"))
saveRDS(max_fam_table, file = paste0(folder, run_name, "max_fam_table.rds"))

saveRDS(famsize_stats, file = paste0(folder, run_name, "famsizes_stats.rds"))
saveRDS(max_fam_stats_long, file = paste0(folder, run_name, "max_fam_stats_long.rds"))
saveRDS(prim_response_stats, file = paste0(folder, run_name, "prim_response_stats.rds"))
saveRDS(sec_response_stats, file = paste0(folder, run_name, "sec_response_stats.rds"))
saveRDS(ter_response_stats, file = paste0(folder, run_name, "ter_response_stats.rds"))

famsizes_freq_table <- generate_freq_famsize_table_multidays_multi(df_famsizes)
saveRDS(famsizes_freq_table, file = paste0(folder, run_name, "famsizes_freq_table.rds"))


