quality_prim <- prim_parameters %>% group_by(fam_nr) %>%
  summarize(quality = mean(quality))

max_fam_table$quality <- quality_prim$quality

  
plot_prim_sec_quality <- ggplot(data = max_fam_table, 
               aes(x = log10(cells_prim), y = log10(cells_sec))) +
  geom_point(aes(col = 11.13*quality, shape = as.factor(nr_burst_divs)), 
             size = 1, alpha = 0.75) +
  scale_color_viridis_c(option = "viridis", direction = -1,
                        limits = c(0, 12),
                        labels = seq(0, 12, 2),
                        breaks = seq(0, 12, 2),
                        alpha = 0.75) +
  guides(colour = guide_colourbar(show.limits = TRUE, 
                                  title.position = "top",
                                  barwidth = 10,
                                  barheight = 0.5)) +
  labs(color = "DD") +
  labs(title = "Primary vs. Secondary Response",
       x = "Family size primary (log-scale)",
       y = "Family size secondary (log-scale)",
       shape = "Nr. of burst divisions (prim. resp.)") +
  theme_clean() +
  th + 
  theme(legend.box.background = element_rect(color="lightgrey")) +
  guides(size = "none", shape = guide_legend(title.position = "top"))  +
  scale_shape_manual(values = c(20, 17, 15)) +
  scale_x_continuous(limits = c(0, 6)) +
  scale_y_continuous(limits = c(0, 6))

ggsave(paste0(folder, run_name, "_plot_prim_sec_quality.pdf"), plot = plot_prim_sec_quality, 
       width = 1500, height = 1300, units = "px")

#-------------------------------------------------------------------------------

  
plot_sec_ter_quality <- ggplot(data = max_fam_table, 
               aes(x = log10(cells_sec), y = log10(cells_ter)))  +
  geom_point(aes(col = 11.13*quality, shape = as.factor(nr_burst_divs)), 
             size = 1, alpha = 0.75) +
  scale_color_viridis_c(option = "viridis", direction = -1,
                        limits = c(0, 12),
                        labels = seq(0, 12, 2),
                        breaks = seq(0, 12, 2),
                        alpha = 0.75) +
  guides(colour = guide_colourbar(show.limits = TRUE, 
                                  title.position = "top",
                                  barwidth = 10,
                                  barheight = 0.5)) +
  labs(color = "DD") +
  labs(
    title = "Secondary vs. Tertiary Response",
    x = "Family size secondary (log-scale)",
    y = "Family size tertiary (log-scale)") +
  theme_clean() +
  th +  
  guides(size = "none", color = "none", 
         shape = "none") +
  scale_shape_manual(values = c(20, 17, 15)) +
  scale_x_continuous(limits = c(0, 6)) +
  scale_y_continuous(limits = c(0, 6))

ggsave(paste0(folder, run_name, "_plot_sec_ter_quality.pdf"), plot = plot_sec_ter_quality, 
       width = 1500, height = 1300, units = "px")
  

