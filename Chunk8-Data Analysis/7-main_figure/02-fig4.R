library(stringr)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(patchwork)
library(ggpubr)
library(EnvStats)

################################################################################
############   Figure 4 --- Accuracy Scores in Diff Data Types  ################
################################################################################
overall_data <- readRDS("./Chunk8-Data Analysis/overall_data.rds") %>% 
  select(1:67)
colnames(overall_data)[40:67] <- str_split(colnames(overall_data)[40:67], "_", simplify = TRUE)[, 2]
method_class_colors <- c(RColorBrewer::brewer.pal(12, "Set3")[4],
                         "#94D3C7","#AAAAE0","#ADCE65","#EAC5E3")
platforms <- c("MARS-seq",
               "10X Genomics",
               "Smart-seq2",
               "Mix sources1",
               "CEL-seq",
               "Fluidigm C1",
               "CEL-seq2",
               "Mix sources2",
               "Drop-seq",
               "inDrop",
               "STRT-seq",
               "Microwell-seq",
               "Smart-seq",
               "ST",
               "HDST",
               "10X Visium",
               "Slide-seq",
               "Slide-seqV2",
               "seqFISH",
               "seqFISH+",
               "osmFISH",
               "sci-Space",
               "MERFISH",
               "Stereo-seq")
###--------------------------------------------------------------------------###
###                           Supplementary Figure 2
###--------------------------------------------------------------------------###
supp_fig2a <- overall_data %>% 
  select(c(1,2, 17:24)) %>% 
  pivot_longer(cols = 3:ncol(.), names_to = "metric", values_to = "value")
supp_fig2a$Category <- factor(supp_fig2a$Category, levels = paste0("Class ", 1:5))
supp_fig2a$metric <- factor(supp_fig2a$metric, levels = c("MAD", "MAE", "RMSE", "KS", "OV", "bhattacharyya", "KDE", "multiKS"))
supp_f2a <- ggboxplot(data = supp_fig2a,
                      x = "Category",
                      y = "value",
                      color = "Category",
                      ggtheme = theme_pubr(),
                      size = 0.2,
                      alpha = 0.6,
                      width = 0.8,
                      bxp.errorbar = FALSE,
                      outlier.size = 0,
                      add = "jitter",
                      add.params = list(size = 0.05))+
  ylab("Metric values") +
  scale_color_manual(values = method_class_colors) +
  scale_fill_manual(values = method_class_colors) +
  facet_wrap(.~ metric, ncol = 4, strip.position = "top") +
  stat_n_text(size = 1, vjust = 0, color = "black") +
  theme(axis.text = element_text(size = 4),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title = element_text(size = 5),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.length = unit(0.05, 'cm'),
        axis.line = element_line(linewidth = 0.2),
        legend.text = element_text(size = 4),
        legend.position = c(0.5, 1.14),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.title = element_text(size = 4),
        legend.key.size = unit(0.5, 'cm'),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 4),
        strip.text.x = element_text(margin = margin(-0.01,0,-0.01,0, "cm")),
        plot.margin = unit(c(0.5,0.1,0.1,0.1), 'cm')) +
  ylim(0, 1) +
  scale_y_continuous(labels = c("0", "0.8"), breaks = c(0, 0.8))


###--------------------------------------------------------------------------###
###                                    supp_fig2b
###--------------------------------------------------------------------------###
supp_fig2b <- overall_data %>% 
  select(c(1,2, 25:39)) %>% 
  pivot_longer(cols = 3:ncol(.), names_to = "property", values_to = "value")
supp_fig2b$Category <- factor(supp_fig2b$Category, levels = paste0("Class ", 1:5))
supp_fig2b$property <- factor(supp_fig2b$property)

supp_f2b <- ggboxplot(data = supp_fig2b,
                  x = "Category",
                  y = "value",
                  color = "Category",
                  ggtheme = theme_pubr(),
                  size = 0.2,
                  alpha = 0.6,
                  width = 0.8,
                  bxp.errorbar = FALSE,
                  outlier.size = 0,
                  add = "jitter",
                  add.params = list(size = 0.05))+
  ylab("Metric values") +
  scale_color_manual(values = method_class_colors) +
  scale_fill_manual(values = method_class_colors) +
  facet_wrap(.~ property, ncol = 4, strip.position = "top") +
  stat_n_text(size = 1, vjust = 0, color = "black") +
  theme(axis.text = element_text(size = 4),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title = element_text(size = 5),
        axis.title.x = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.length = unit(0.05, 'cm'),
        axis.line = element_line(linewidth = 0.2),
        legend.position = "none",
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 4),
        strip.text.x = element_text(margin = margin(-0.01,0,-0.01,0, "cm")),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), 'cm')) +
  ylim(0, 1) +
  scale_y_continuous(labels = c("0", "0.8"), breaks = c(0, 0.8))

ggsave(plot = supp_f2a / supp_f2b + plot_layout(height = c(1,2)),
       filename = "../sim-article/figures/Supp_Fig_2.pdf",
       width = 10,
       height = 10,
       units = "cm")


###--------------------------------------------------------------------------###
###                                    Fig4a
###--------------------------------------------------------------------------###
technique_colors <- c("#7DAEF4", "#E0709E")
method_class_colors <- c(RColorBrewer::brewer.pal(12, "Set3")[4],
                         "#94D3C7","#AAAAE0","#ADCE65","#EAC5E3")

arrange_by_group <- function(tibble){
  groups <- unique(tibble %>% pull("Category"))
  result <- tibble()
  for(i in groups){
    tmp <- tibble %>% 
      filter(Category == i) %>% 
      arrange(desc(accuracy))
    result <- rbind(result, tmp)
  }
  return(result)
}
per_platform_score <- overall_data[, c("id", "Category", "accuracy", platforms)]
per_platform_score <- arrange_by_group(per_platform_score) %>% 
  select(-2,-3) %>% 
  tibble::column_to_rownames("id")

#### heatmap
library(ComplexHeatmap)

### colume metadata
meta_data <- data.frame(
  `Technique Platform` = c(rep("scRNA-seq", 13), rep("ST technology", 11))
)
rownames(meta_data) <- colnames(per_platform_score)

### row metadata
row_meta <- data.frame(
  `class` = c(rep("class1", 15),
              rep("class2", 10),
              rep("class3", 9),
              rep("class4", 9),
              rep("class5", 6)),
  `color` = c(rep(method_class_colors[1], 15),
              rep(method_class_colors[2], 10),
              rep(method_class_colors[3], 9),
              rep(method_class_colors[4], 9),
              rep(method_class_colors[5], 6))
)
rownames(row_meta) <- rownames(per_platform_score)

col <- list(Technique.Platform = c(`scRNA-seq` = technique_colors[1],
                                   `ST technology` = technique_colors[2]),
            class = c(class1 = method_class_colors[1],
                      class2 = method_class_colors[2],
                      class3 = method_class_colors[3],
                      class4 = method_class_colors[4],
                      class5 = method_class_colors[5]))

pdf(file = "../sim-article/figures/Fig4a.pdf", width = 7, height = 8)
p4_a <- ComplexHeatmap::pheatmap(per_platform_score %>% as.matrix(),
                                 show_colnames = TRUE,
                                 show_rownames = TRUE,
                                 scale = "none",
                                 cluster_cols = FALSE,
                                 cluster_rows = FALSE,
                                 annotation_col = meta_data,
                                 annotation_colors = col,
                                 annotation_row = row_meta %>% select(1),
                                 treeheight_row = 20,
                                 border_color = "white",
                                 gaps_row = c(15, 25, 34, 43), 
                                 color = c(colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(40)),
                                 na_col = "gray80")
print(p4_a)
dev.off()


###--------------------------------------------------------------------------###
###                                    Fig4b
###--------------------------------------------------------------------------###
### fig4b
library(ggrepel)
library(ggpmisc)
fig4b_data <- overall_data %>% 
  select(all_of(c("id", "scRNA-seq data", "spatial transcriptome data", "Category"))) %>% 
  rename(., `ST technology` = `spatial transcriptome data`) %>% 
  rename(., `scRNA-seq` = `scRNA-seq data`) %>% 
  rename(., Method = id) %>% 
  drop_na()

fig4b <- ggplot(fig4b_data, aes(x = `ST technology`, y = `scRNA-seq`)) +
  geom_smooth(method = 'lm', linewidth = 0.6, color = "red") +
  geom_point(size = 0.8) +
  stat_cor(method = "pearson", size = 2) +
  geom_label_repel(aes(label = Method,
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
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(size = 5),
        legend.title = element_text(size = 4),
        legend.text = element_text(size = 4),
        legend.background = element_blank(),
        panel.grid = element_blank(),
        legend.key.size = unit(0.2, 'cm')) +
  xlab("Accuracy scores on the spatial datasets") +
  ylab("Accuracy scores on the scRNA-seq datasets") +
  scale_fill_manual(values = method_class_colors)

ggsave(plot = fig4b,
       filename = "../sim-article/figures/Fig4b.pdf",
       width = 7,
       height = 7,
       units = "cm")

###--------------------------------------------------------------------------###
###                                Fig4c
###--------------------------------------------------------------------------###
library(ggrepel)
library(ggpmisc)
spatial_accuracy_data <- readRDS("Chunk8-Data Analysis/2-accuracy/spatial_accuracy_plot_data.rds")
fig4c <- fig4b_data %>% 
  select(-3) %>% 
  left_join(spatial_accuracy_data %>% select(1:2), by = "Method") %>% 
  rename(., `ST technology` = accuracy)

fig4c <- ggplot(fig4c, aes(x = `ST technology`, y = `scRNA-seq`)) +
  geom_smooth(method = 'lm', linewidth = 0.3, color = "red") +
  geom_point(size = 0.8) +
  stat_cor(method = "pearson", size = 2) +
  geom_label_repel(aes(label = Method,
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
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(size = 5),
        legend.title = element_text(size = 4),
        legend.text = element_text(size = 4),
        legend.background = element_blank(),
        panel.grid = element_blank(),
        legend.key.size = unit(0.2, 'cm')) +
  xlab("Accuracy scores on the spatial datasets \n (calculated by spatial-level metrics)") +
  ylab("Accuracy scores on the scRNA-seq datasets") +
  scale_fill_manual(values = method_class_colors)

ggsave(plot = fig4c,
       filename = "../sim-article/figures/Fig4c.pdf",
       width = 7,
       height = 7,
       units = "cm")



###--------------------------------------------------------------------------###
###                                    Fig4d
###--------------------------------------------------------------------------###
#### Fig4d
Fig4d_data <- spatial_accuracy_data
colnames(Fig4d_data)[9:19] <- str_split(colnames(Fig4d_data)[9:19], "_", simplify = TRUE)[, 2]
Fig4d_data <- overall_data %>% 
  select(all_of(c("id", platforms[1:13]))) %>% 
  rename(Method = id) %>% 
  full_join(Fig4d_data, by = "Method")
Fig4d_data <- Fig4d_data %>% 
  select(-c(15:21)) %>% 
  pivot_longer(cols = 2:ncol(.), names_to = "Platform", values_to = "platform_score")
Fig4d_data$Method <- factor(Fig4d_data$Method, levels = overall_data$id)
Fig4d_data <- Fig4d_data %>% 
  mutate(
    technology = case_when(
      Platform %in% platforms[1:13] ~ "scRNA-seq",
      Platform %in% platforms[14:24] ~ "ST technology"
    )
  )
odds <- seq(1, 49, 2)
rect <- Fig4d_data[odds, ]
fig4d <- ggplot(Fig4d_data, aes(x = Method, y = platform_score, color = technology))+
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
  geom_jitter(aes(group = technology),
              pch = 19,
              size = 0.05,
              position = position_jitterdodge(dodge.width=0.9)) +
  stat_compare_means(mapping = aes(group = technology),
                     method = "wilcox.test",
                     label = "p.signif",
                     hide.ns = TRUE) +
  scale_color_manual(values = technique_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45,
                               vjust = 1,
                               hjust = 1,
                               size = 4.5,
                               color = "black"),
    axis.text.y = element_text(size = 4, color = "black"),
    axis.title.y = element_text(size = 5),
    axis.ticks = element_line(linewidth = 0.05, lineend = 0.01, color = "black"),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 4),
    legend.position = c(0.5, -0.5),
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.box = "horizontal",
    legend.key = element_rect(fill = NA),
    legend.key.size = unit(0.3, 'cm'),
    panel.grid = element_blank()
  )

ggsave(plot = fig4d,
       filename = "../sim-article/figures/Fig4d.pdf",
       width = 18,
       height = 6,
       units = "cm")

###--------------------------------------------------------------------------###
###                             Supplementary Fig 6
###--------------------------------------------------------------------------###
accuracy <- readRDS("Chunk8-Data Analysis/2-accuracy/accuracy_long_data.rds")
data_info <- openxlsx::read.xlsx("./Chunk1-Data preparation/evaluation_datasets.xlsx", rowNames = TRUE, sep.names = " ") %>% 
  mutate(
    Data = str_split(`Dataset ID`, "_", simplify = TRUE)[, 1]
  ) %>% 
  select(Data, Platform, `Quantification Strategy`)

accuracy <- accuracy %>% 
  full_join(data_info, by = "Data") %>% 
  relocate(Platform, .after = "Data") %>% 
  mutate(
    Platform = case_when(
      Platform == "Smart-seq2\r\n10X Genomics" ~ "Mix sources1",
      Platform == "CEL-seq\r\nCEL-seq2" ~ "Mix sources2",
      TRUE ~ Platform
    ),
    Type = case_when(
      Data %in% paste0("data", 1:101) ~ "scRNA-seq data",
      Data %in% paste0("data", 101:152) ~ "spatial transcriptome data"
    )
  ) %>% 
  relocate(Type, .after = "Platform")

supp6_data <- accuracy %>% 
  group_by(Method, metric, Platform, Type, Data, `Quantification Strategy`) %>% 
  summarise(
    value = mean(value, na.rm = TRUE),
  ) %>% 
  group_by(Method, metric, `Quantification Strategy`, Type, Data) %>% 
  summarise(
    value = mean(value, na.rm = TRUE),
  ) %>% 
  filter(Type == "scRNA-seq data") %>% 
  group_by(Method, `Quantification Strategy`, Data) %>% 
  summarise(
    value = mean(value, na.rm = TRUE),
  ) %>% 
  drop_na()

supp6_data$Method <- factor(supp6_data$Method, levels = overall_data$id)
odds <- seq(1, 47, 2)
rect <- supp6_data[odds, ]

sample_size <- supp6_data %>%
  drop_na() %>% 
  group_by(Method, `Quantification Strategy`) %>%
  summarize(num = n()) %>% 
  mutate(y = 0)

supp6 <- ggplot(supp6_data, aes(x = Method, y = value, color = `Quantification Strategy`))+
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
  geom_jitter(aes(group = `Quantification Strategy`),
              pch = 19,
              size = 0.05,
              position = position_jitterdodge(dodge.width=0.9)) +
  geom_text(data = sample_size, aes(x = Method, y = y, label = num), size = 1) +
  stat_compare_means(mapping = aes(group = `Quantification Strategy`),
                     method = "wilcox.test",
                     label = "p.signif",
                     hide.ns = TRUE) +
  # scale_color_manual(values = technique_colors) +
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
    legend.position = c(0.5, -0.3),
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.box = "horizontal",
    legend.key = element_rect(fill = NA),
    legend.key.size = unit(0.3, 'cm'),
    panel.grid = element_blank()
  ) +
  ylim(0, 1)

ggsave(plot = supp6,
       filename = "../sim-article/figures/Supp_Fig_6.pdf",
       width = 18,
       height = 8,
       units = "cm")
