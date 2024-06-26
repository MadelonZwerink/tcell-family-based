library(grid)

# This is to select Q cells for the next response

# div_cells_prim contains the number of cells in each division index for all
# branches (each branch is a row). We can use this dataframe to recruit cells
# with a probability based on their division index

#-------------------------------------------------------------------------------
# First, all cells have to be converted to whole numbers, because we can't 
# recruit half a cell

div_cells_prim_rounded <- div_cells_prim

for (r in 1:nrow(div_cells_prim_rounded)){
  for (c in 1:n){
    numb <- runif(1)
    dec <- div_cells_prim_rounded[r, c+1] %% 1
    if(numb <= dec){
      div_cells_prim_rounded[r, c+1] <- 1 + as.integer(div_cells_prim_rounded[r, c+1])
    }
    if(numb >= dec){
      div_cells_prim_rounded[r, c+1] <- 0 + as.integer(div_cells_prim_rounded[r, c+1])
    }
  }
}

#-------------------------------------------------------------------------------
# Now, we can select each cell with a certain probability depending on 
# division index
recruitment_mean <- 5
recruitment_sd <- 1

recruitment_prob <- pnorm(0:n, mean = recruitment_mean, 
                          sd = recruitment_sd, lower.tail = FALSE)
prob_df <- as.data.frame(cbind(div_index = seq(0,n), prob = recruitment_prob))
plot_prob_dist <- ggplot(prob_df, aes(x = div_index, y = prob)) + geom_point() + xlim(0,25) +
  labs(x = "Division index",
       y = "Recruitment probability") +
  theme_clean() + th

#-------------------------------------------------------------------------------
# Check the number of Q cells and the division index of Q cells for different 
# recruitment probability distributions

# Select cells depending on probability
selected_Q_cells <- data.frame(matrix(0, nrow = nrow(div_cells_prim_rounded), 
                                      ncol = ncol(div_cells_prim_rounded)))
colnames(selected_Q_cells) <- colnames(div_cells_prim_rounded)
selected_Q_cells$fam_nr <- div_cells_prim_rounded$fam_nr

for (r in 1:nrow(selected_Q_cells)){
  for (c in 1:n){
    numb <- runif(div_cells_prim_rounded[r, c+1])
    Q_cells <- sum(as.numeric(numb <= recruitment_prob[c]))
    selected_Q_cells[r, c+1] <- Q_cells
  }
}

selected_Q_cells_grouped <- selected_Q_cells %>%
  wide_to_long("Q_cells") %>%
  group_by(fam_nr, div_index) %>%
  summarise(Q_cells = sum(Q_cells))

Q_cells_per_div_index <- data.frame(Q_cells = colSums(selected_Q_cells)[-1], 
                                    div_index = seq(0, n - 1))

Q_cells_per_family <- selected_Q_cells_grouped %>%
  group_by(fam_nr) %>%
  summarise(Q_cells = sum(Q_cells))

plot_dist_per_fam <- ggplot(Q_cells_per_family) + 
  geom_boxplot(aes(y = Q_cells))  +
  labs(title = "Q cells per family",
       x = element_blank(),
       y = "Number of Q cells") +
  theme_clean() + th + theme(axis.text.x = element_blank())

# Calculate total number of Q cells and mean div index
total_nr_Q_cells <- sum(selected_Q_cells_grouped$Q_cells)
mean_div_index <- round(sum(Q_cells_per_div_index$Q_cells * Q_cells_per_div_index$div_index)/total_nr_Q_cells, digits = 2)

# This plots the number of Q cells that are formed as a function of div index
plot_Q_cells_total <- ggplot(Q_cells_per_div_index, aes(x = div_index, y = Q_cells)) +
  geom_bar(stat = "identity") + xlim(0,25) +
  labs(title = "Total number of Q cells",
       x = "Division index",
       y = "Number of Q cells") +
  theme_clean() + th

plot_Q_per_fam <- ggplot(selected_Q_cells_grouped, aes(x = div_index, 
                                   y = Q_cells,
                                   col = as.factor(fam_nr), 
                                   group = fam_nr,
                                   alpha = 0.5)) +
  labs(title = "Number of Q cells per family",
       x = "Division index",
       y = "Number of Q cells") +
  geom_line() +
  theme_clean() + th +
  theme(legend.position = "none") + xlim(0, 25)

stats <- textGrob(label = paste("Recruitment function:\nMean:", recruitment_mean,
                                "sd:", recruitment_sd,
                                "\nTotal cells:", total_nr_Q_cells,
                                "\nMean div. index:", mean_div_index),
                 x = unit(0.10, "npc"),
                 hjust = 0)

gt_recruitment_prob <- arrangeGrob(grobs = list(plot_prob_dist,
                                                plot_Q_cells_total,
                                                plot_Q_per_fam,
                                                stats,
                                                plot_dist_per_fam), 
                                    nrow = 3, layout_matrix = rbind(c(1,1,5,4),
                                                                    c(2,2,3,3),
                                                                    c(2,2,3,3)),
                                    padding = unit(0, "line"))

panel_plots_prob <- as_ggplot(gt_recruitment_prob) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C", "D"), size = 13,
                  x = c(0, 0.5, 0, 0.5), y = c(1, 1, 0.66, 0.66))
panel_plots_prob

selected_Q_cells_long <- wide_to_long(selected_Q_cells, "Q_cells")
selected_Q_cells_long <- selected_Q_cells_long[selected_Q_cells_long$Q_cells != 0,]

selected_Q_cells_processed <- selected_Q_cells_long %>%
  uncount(Q_cells)
selected_Q_cells_processed$cell_type <- "Q"
colnames(selected_Q_cells_processed) <- c("fam_nr", "div_counter", "cell_type")

sec_parameters <- pick_parameters(response_nr = 2,
                                  prev_parameters = selected_Q_cells_processed,
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


