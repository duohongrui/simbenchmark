###--------------------------------------------------------------------------###
###                         Supplementary Figure 7
###--------------------------------------------------------------------------###
library(tidyverse)
library(ggrepel)
library(ggpmisc)
overall_data <- readRDS("Chunk8-Data Analysis/overall_data.rds")
method_class_colors <- c(RColorBrewer::brewer.pal(12, "Set3")[4],
                         "#94D3C7","#AAAAE0","#ADCE65","#EAC5E3")
p <- ggplot(overall_data, aes(x = accuracy, y = functionality)) +
  geom_smooth(method = 'lm', linewidth = 0.6, color = "red") +
  geom_point(size = 0.8) +
  stat_cor(method = "pearson", size = 2) +
  geom_label_repel(aes(label = id,
                       fill = Category),
                   color = "black",
                   size = 1,
                   box.padding = unit(0.4, "lines"),
                   point.padding = unit(0.1, "lines"),
                   label.padding = unit(0.1, "lines"),
                   segment.colour = "black",
                   max.overlaps = 50,
                   segment.size = 0.2) +
  theme_bw() +
  theme(legend.position = c(0.1, 0.86),
        axis.text = element_text(size = 4, color = "black"),
        axis.title = element_text(size = 5),
        legend.title = element_text(size = 4),
        legend.text = element_text(size = 4),
        legend.background = element_blank(),
        panel.grid = element_blank(),
        legend.key.size = unit(0.2, 'cm')) +
  scale_fill_manual(values = method_class_colors) +
  xlab("Accuracy scores of accuracy") +
  ylab("Accuracy scores of functionality")

ggsave(plot = p,
       filename = "../sim-article/figures/Supp_Fig_7.pdf",
       width = 7,
       height = 7,
       units = "cm")
