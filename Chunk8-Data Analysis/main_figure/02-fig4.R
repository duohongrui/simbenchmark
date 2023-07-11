library(stringr)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(patchwork)
library(ggpubr)

################################################################################
############   Figure 4 --- Accuracy Scores in Diff Data Types  ################
################################################################################
overall_data <- readRDS("./Chunk8-Data Analysis/overall_data.rds")
acc_data <- readRDS("./Chunk8-Data Analysis/accuracy/accuracy_long_data.rds")
data_info <- openxlsx::read.xlsx("./Chunk1-Data preparation/evaluation_datasets.xlsx", rowNames = TRUE) %>% 
  mutate(
    Data = str_split(Dataset.ID, "_", simplify = TRUE)[, 1]
  ) %>% 
  select(Data, Platform)
platforms <- c("MARS-Seq",
               "10X Genomics",
               "Smart-seq2",
               "Mix sources1",
               "CEL-seq",
               "Fluidigm C1",
               "CEL-seq2",
               "Mix sources2",
               "Drop-seq",
               "inDrop",
               "Microwell-seq",
               "Smart-seq",
               "ST",
               "HDST",
               "10X Visium",
               "Slide-Seq",
               "Slide-SeqV2",
               "seqFISH",
               "seqFISH+",
               "osmFISH",
               "sci-Space",
               "MERFISH",
               "Stereo-Seq")
### add platform and type to accuracy table
acc_data <- acc_data %>% 
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

####################### Platform
technique_colors <- RColorBrewer::brewer.pal(8, "Set2")[1:2]
method_class_colors <- c("#F97369", "#69ACD1", "#F9A24B", "#ADDD52", "#F9B4DA")

platform_data <- acc_data %>% 
  group_by(Platform, metric, Method) %>% 
  summarise(
    value = mean(value, na.rm = TRUE)
  ) %>% 
  ungroup()

### per_platform_score
per_platform_score <- platform_data %>% 
  group_by(Platform, Method) %>% 
  summarise(
    value = mean(value, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = everything(), names_from = "Platform", values_from = "value") %>% 
  column_to_rownames(var = "Method")
per_platform_score <- per_platform_score[, platforms]

#### heatmap
library(ComplexHeatmap)
library(pheatmap)
per_platform_score <- per_platform_score[overall_data$id, ]

### colume metadata
meta_data <- data.frame(
  `Technique Platform` = c(rep("scRNA-seq", 12), rep("spatial transcriptomics technique", 11))
)
rownames(meta_data) <- colnames(per_platform_score)

### row metadata
row_meta <- data.frame(
  `class` = c(rep("class1", 10),
              rep("class2", 13),
              rep("class3", 9),
              rep("class4", 7),
              rep("class5", 6)),
  `color` = c(rep(method_class_colors[1], 10),
              rep(method_class_colors[2], 13),
              rep(method_class_colors[3], 9),
              rep(method_class_colors[4], 7),
              rep(method_class_colors[5], 6))
)
rownames(row_meta) <- rownames(per_platform_score)

col <- list(Technique.Platform = c(`scRNA-seq` = technique_colors[1],
                                   `spatial transcriptomics technique` = technique_colors[2]),
            class = c(class1 = method_class_colors[1],
                      class2 = method_class_colors[2],
                      class3 = method_class_colors[3],
                      class4 = method_class_colors[4],
                      class5 = method_class_colors[5]))

pdf(file = "../sim-article/figures/Fig4-a.pdf", width = 8, height = 8)
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
                                 gaps_row = c(10, 23, 32, 39), 
                                 color = c(colorRampPalette(c("#19547b","#ffd89b"))(40)),
                                 na_col = "grey80",
                                 right_annotation = rowAnnotation(bar1 = anno_barplot(apply(per_platform_score[, 1:12], 1, mean, na.rm = TRUE),
                                                                                      gp = gpar(fill = row_meta$color)),
                                                                  bar2 = anno_barplot(apply(per_platform_score[, 13:23], 1, mean, na.rm = TRUE),
                                                                                      gp = gpar(fill = row_meta$color))))
print(p4_a)
dev.off()


#### Fig4-S
colors <- RColorBrewer::brewer.pal(12, "Set3")[4:8]
methods <- openxlsx::read.xlsx("./Chunk1-Data preparation/methods.xlsx")

### normalize every metrics to [0, 1]
platform_p4s_data <- platform_data %>% 
  group_by(Platform, Method) %>% 
  summarise(
    platform_score = mean(value, na.rm = TRUE)
  )%>%
  ungroup() %>% 
  full_join(methods %>% select(1:2), by = "Method")
platform_p4s_data$Method <- factor(platform_p4s_data$Method, levels = overall_data$id)
P4_s <- ggbarplot(data = platform_p4s_data,
                  x = "Method",
                  y = "platform_score",
                  color = "Category",
                  fill = "Category",
                  ggtheme = theme_pubr(),
                  size = 1,
                  alpha = 0.6)+
  ylab("Metric values") +
  scale_color_manual(values = method_class_colors) +
  scale_fill_manual(values = colors) +
  facet_grid(Platform ~.) +
  ylim(0, 1) +
  scale_y_continuous(labels = c("0", "0.8"), breaks = c(0, 0.8)) +
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = "top",
        legend.background = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = NA, color = "black"),
        strip.text = element_text(size = 8, angle = 45)) +
  expand_limits(x = c(0, 48))

ggsave(plot = P4_s, filename = "../sim-article/figures/Fig4-S.pdf", width = 15, height = 14, units = "in")


#### Fig4b boxplot
Fig4b_data <- platform_p4s_data
Fig4b_data$Category <- factor(Fig4b_data$Category, levels = paste0("Class ", 1:5))
P4_b <- ggboxplot(data = Fig4b_data,
                  x = "Category",
                  y = "platform_score",
                  color = "Category",
                  ggtheme = theme_pubr(),
                  size = 0.2,
                  alpha = 0.6,
                  width = 0.8,
                  bxp.errorbar = FALSE,
                  outlier.size = 0,
                  add = "jitter",
                  add.params = list(size = 0.15))+
  ylab("Metric values") +
  scale_color_manual(values = c("#F97369", "#69ACD1", "#F9A24B", "#ADDD52", "#F9B4DA")) +
  scale_fill_manual(values = colors) +
  facet_wrap(.~ Platform, ncol = 4, strip.position = "top") +
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = "top",
        legend.background = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = "white", color = "white")) +
  ylim(0, 0.9) +
  scale_y_continuous(labels = c("0.1", "0.8"), breaks = c(0.1, 0.8))

ggsave(plot = P4_b, filename = "../sim-article/figures/Fig4-b.pdf", width = 6, height = 8, units = "in")


#### Fig4c boxplot
Fig4c_data <- platform_p4s_data %>% 
  mutate(
    technology = case_when(
      Platform %in% platforms[1:12] ~ "scRNA-seq",
      Platform %in% platforms[13:23] ~ "ST technology"
    )
  )
odds <- seq(1, 45, 1)
rect <- Fig4c_data[odds, ]
fig4c <- ggplot(Fig4c_data, aes(x = Method, y= platform_score, fill = technology))+
  geom_rect(data = rect,
            xmin = odds - 0.5,
            xmax = odds + 0.5,
            ymin = -Inf,
            ymax = +Inf,
            fill = row_meta$color,
            alpha = 0.3,
            inherit.aes = F) +
  geom_boxplot(outlier.alpha = 0,
               width = 0.5,
               size = 0.1) +
  geom_jitter(data = Fig4c_data,
              aes(x = Method, y = platform_score),
              pch = 19,
              size = 0.2) +
  theme_bw() +
  scale_fill_manual(values = technique_colors) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.box = "horizontal",
    legend.key = element_rect(fill = NA),
    panel.grid = element_blank(),
    axis.title.x = element_blank()
  ) +
  stat_compare_means(data = Fig4c_data,
                     mapping = aes(group = technology),
                     method = "wilcox.test",
                     label = "p.signif",
                     hide.ns = TRUE)

ggsave(plot = fig4c, filename = "../sim-article/figures/Fig4-c.pdf", width = 18, height = 6, units = "in")


### fig4d
library(ggrepel)
fig4d_data <- Fig4c_data %>% 
  group_by(Method, technology) %>% 
  summarise(
    technology_score = mean(platform_score, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = everything(),
              names_from = technology,
              values_from = technology_score) %>% 
  right_join(methods %>% select(1:2), by = "Method") %>% 
  drop_na()

cor_value <- sprintf("%0.2f", cor(fig4d_data$`ST technology`, fig4d_data$`scRNA-seq`))

fig4d <- ggplot(fig4d_data, aes(x = `ST technology`, y = `scRNA-seq`)) +
  geom_point(size = 3) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_label_repel(aes(label = Method,
                       fill = Category),
                   color = "white",
                   size = 4,
                   box.padding = unit(0.6, "lines"),
                   point.padding = unit(0.8, "lines"),
                   segment.colour = "black",
                   max.overlaps = 50) +
  annotate(geom = "text",
           x = 0.25,
           y = 0.7,
           label = paste0("Cor: ", cor_value),
           size = 5) +
  theme_bw() +
  theme(legend.position = c(0.1, 0.8),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 10),
        panel.grid = element_blank()) +
  xlab("Accuracy scores on the spatial datasets") +
  ylab("Accuracy scores on the scRNA-seq datasets") +
  scale_fill_manual(values = method_class_colors)

ggsave(plot = fig4d, filename = "../sim-article/figures/Fig4-d.pdf", width = 9, height = 9, units = "in")


### fig4e
fig4e_data <- Fig4c_data
fig4e <- ggplot(fig4e_data, aes(x = Category, y = platform_score, fill = technology)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(data = fig4e_data,
              aes(x = Category, y = platform_score),
              pch = 19,
              size = 0.2) +
  theme_bw() +
  scale_fill_manual(values = technique_colors) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.box = "horizontal",
    legend.key = element_rect(fill = NA),
    panel.grid = element_blank(),
    axis.title.x = element_blank()
  ) +
  stat_compare_means(data = fig4e_data,
                     mapping = aes(group = technology),
                     method = "wilcox.test",
                     label = "p.signif",
                     hide.ns = TRUE)

ggsave(plot = fig4e, filename = "../sim-article/figures/Fig4-e.pdf", width = 8, height = 8, units = "in")


