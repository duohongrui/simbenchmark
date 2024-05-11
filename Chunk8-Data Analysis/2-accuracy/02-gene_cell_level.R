################################################################################
##########################   Supplementary Figure 4  ###########################
################################################################################
library(tidyverse)
library(ggpubr)
method_class_colors <- c(RColorBrewer::brewer.pal(12, "Set3")[4],
                         "#94D3C7","#AAAAE0","#ADCE65","#EAC5E3")
overall_data <- readRDS("Chunk8-Data Analysis/overall_data.rds")

cell_level <- c("LS","FZC","CCC","TMM","ELS","FCO","RLZ")
gene_level <- c("ME","SD","CV","FZG","FGO","RMS","RMZ", "RDM")

cell_data <- overall_data %>% 
  select(all_of(c("id", cell_level))) %>% 
  pivot_longer(cols = all_of(cell_level), names_to = "property", values_to = "score") %>% 
  mutate(
    level = "cell level"
  )

gene_data <- overall_data %>% 
  select(all_of(c("id", gene_level))) %>% 
  pivot_longer(cols = all_of(gene_level), names_to = "property", values_to = "score") %>% 
  mutate(
    level = "gene level"
  )

level <- cell_data %>% 
  rbind(gene_data)


### row metadata
row_meta <- data.frame(
  `class` = c(rep("class1", 10),
              rep("class2", 15),
              rep("class3", 9),
              rep("class4", 9),
              rep("class5", 6)),
  `color` = c(rep(method_class_colors[1], 10),
              rep(method_class_colors[2], 15),
              rep(method_class_colors[3], 9),
              rep(method_class_colors[4], 9),
              rep(method_class_colors[5], 6))
)
rownames(row_meta) <- overall_data$id

level$id <- factor(level$id, levels = overall_data$id)
odds <- seq(1, 49, 2)
rect <- level[odds, ]
sample_size <- level %>%
  drop_na() %>% 
  group_by(id, level) %>%
  summarize(num = n()) %>% 
  mutate(y = 0)
p <- ggplot(level, aes(x = id, y = score, color = level))+
  geom_rect(data = rect,
            xmin = odds - 0.5,
            xmax = odds + 0.5,
            ymin = -Inf,
            ymax = +Inf,
            fill = "gray60",
            alpha = 0.3,
            inherit.aes = F) +
  geom_boxplot(outlier.alpha = 0,
               width = 1,
               size = 0.1) +
  geom_jitter(aes(group = level),
              pch = 19,
              size = 0.05,
              position = position_jitterdodge(dodge.width=0.9)) +
  stat_compare_means(mapping = aes(group = level),
                     method = "wilcox.test",
                     label = "p.signif",
                     hide.ns = TRUE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45,
                               vjust = 1,
                               hjust = 1,
                               size = 5.5,
                               color = "black"),
    axis.text.y = element_text(size = 4, color = "black"),
    axis.title.y = element_text(size = 5),
    axis.ticks = element_line(linewidth = 0.05, lineend = 0.01, color = "black"),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 4),
    legend.position = c(0.5, -0.25),
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.box = "horizontal",
    legend.key = element_rect(fill = NA),
    legend.key.size = unit(0.3, 'cm'),
    panel.grid = element_blank()
  )

ggsave(plot = p,
       filename = "../sim-article/figures/Supp_Fig_3.pdf",
       width = 18,
       height = 8,
       units = "cm")
