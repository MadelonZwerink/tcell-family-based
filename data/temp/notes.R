# notes script
#-------------------------------------------------------------------------------

famsizes_table_agg <- aggregate(. ~ fam_nr + fam_nr_2 + fam_nr_3, data = famsizes_table, sum) 

#-------------------------------------------------------------------------------

# Make the same plots as the toy model to check whether the dynamics are the
# same, and they are
prim <- table(prim_parameters[,c("fam_nr", "cell_type")])
sec <- table(sec_parameters[,c("fam_nr", "cell_type")])
ter <- table(ter_parameters[,c("fam_nr", "cell_type")])
all <- cbind(prim, sec, ter)
all <- data.frame(all)
colnames(all) <- c("P1", "Q1", "P2", "Q2", "P3", "Q3")
plot(all$Q1, all$Q2)
plot(all$Q2, all$Q3)
plot(all$P1, all$P2)
plot(all$P2, all$P3)
plot(all$Q1, all$P1)

all <- data.frame(prim)
# Old parameter values (16/4/2024)
sec_parameters <- pick_parameters(bp_rule = 'runif(1, min = 0.5, max = 3)',
                                  dp_rule = 0.4,
                                  rq_rule = 0.5,
                                  prev_parameters = prim_parameters,
                                  response_nr = 2,
                                  t_run_rule = 'runif(1, min = 1, max = 4)',
                                  t_start_dist = 'rlnorm(nr_of_families, meanlog = 1.35, sdlog = 0.2)',
                                  nr_burst_divs = 'sample(c(1,2,3,4), nr_of_families, replace = TRUE)')



