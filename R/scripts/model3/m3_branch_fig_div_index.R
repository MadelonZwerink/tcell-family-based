folder <- "./results/model1/"
run_name <- "m1_"

data_objects <- c("prim_parameters",
                  "sec_parameters",
                  "ter_parameters")
for(data_object in data_objects){
  assign(data_object, readRDS(paste0(folder, run_name, data_object, ".rds"))) 
}

folder <- "./results/model3/new_cluster/"
single_run_name <- "m3_DD_branch_2"

data_objects <- c("recruited_cells_prim",
                  "recruited_cells_sec",
                  "recruited_cells_ter")

for(data_object in data_objects){
  assign(data_object, readRDS(paste0(folder, single_run_name, data_object, ".rds"))) 
}

m3_prim <- data.frame(div_counter = recruited_cells_prim$div_counter,
                      model = 3)
m3_sec <- data.frame(div_counter = recruited_cells_sec$div_counter,
                     model = 3)
m3_ter <- data.frame(div_counter = recruited_cells_ter$div_counter,
                     model = 3)
m1_prim <- data.frame(div_counter = prim_parameters[prim_parameters$cell_type == "Q",]$div_counter,
                      model = 1)
m1_sec <- data.frame(div_counter = sec_parameters[sec_parameters$cell_type == "Q",]$div_counter,
                     model = 1)
m1_ter <- data.frame(div_counter = ter_parameters[ter_parameters$cell_type == "Q",]$div_counter,
                     model = 1)
prim_divs <- rbind(m1_prim, m3_prim)
sec_divs <- rbind(m1_sec, m3_sec)
ter_divs <- rbind(m1_ter, m3_ter)

div_prim <- ggplot(prim_divs, aes(y = div_counter, 
                                  x = factor(model), 
                                  group = factor(model))) +
  geom_boxplot() +
  geom_segment(aes(x = 1.55, y = 8, yend = 8, xend = 2.45), 
               col = "#9f2a63", linewidth = 0.9, linetype = "dashed") +
  theme_clean() + th + ylim(0, 15) + 
  labs(y = "Division index",
       x = element_blank(),
       title = "Primary response") + 
  scale_x_discrete(labels = c("1" = paste("M1\n","n =", nrow(m1_prim)), 
                              "3" = paste("M3\n","n =", nrow(m3_prim))))

div_sec <- ggplot(sec_divs, aes(y = div_counter, 
                                x = factor(model), 
                                group = factor(model))) +
  geom_boxplot() + 
  geom_segment(aes(x = 1.55, y = 8, yend = 8, xend = 2.45), 
               col = "#9f2a63", linewidth = 0.9, linetype = "dashed") +
  theme_clean() + th + ylim(0, 15) + 
  labs(y = "Division index",
       x = element_blank(),
       title = "Secondary response") + 
  scale_x_discrete(labels = c("1" = paste("M1\n","n =", nrow(m1_sec)), 
                              "3" = paste("M3\n","n =", nrow(m3_sec))))

div_ter <- ggplot(ter_divs, aes(y = div_counter, 
                                x = factor(model), 
                                group = factor(model))) +
  geom_boxplot() + 
  geom_segment(aes(x = 1.55, y = 8, yend = 8, xend = 2.45), 
               col = "#9f2a63", linewidth = 0.9, linetype = "dashed") +
  theme_clean() + th + ylim(0, 15) + 
  labs(y = "Division index",
       x = element_blank(),
       title = "Tertiary response") + 
  scale_x_discrete(labels = c("1" = paste("M1\n","n =", nrow(m1_ter)), 
                              "3" = paste("M3\n","n =", nrow(m3_ter))))

gt <- arrangeGrob(grobs = list(div_prim, div_sec, div_ter), 
                  nrow = 1, 
                  padding = unit(0, "line"))

div_plot <- as_ggplot(gt) +             # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C"), 
                  size = 15,
                  x = c(0.01, 0.34, 0.68), 
                  y = c(0.97, 0.97, 0.97))

ggsave("./results/model3/new_cluster/m3_DD_branch_7_div_plot.jpg", plot = div_plot, 
       width = 8, height = 3.5, units = "in")
