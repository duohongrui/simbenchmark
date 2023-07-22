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
overall_data <- readRDS("./Chunk8-Data Analysis/overall_data.rds") %>% 
  select(1:62)
colnames(overall_data)[38:62] <- str_split(colnames(overall_data)[38:62], "_", simplify = TRUE)[, 2]
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
###--------------------------------------------------------------------------###
###                                    Fig4c
###--------------------------------------------------------------------------###
technique_colors <- RColorBrewer::brewer.pal(8, "Set2")[1:2]
technique_colors <- c("#86b3d3", "#f5926e")
method_class_colors <- c("#F97369", "#b4b9f9", "#fec79c", "#92C694", "#F9B4DA")
method_class_colors <- RColorBrewer::brewer.pal(12, "Set3")[4:8]

platform_data <- overall_data %>% 
  select(all_of(c("id", "Category", platforms)))

#### heatmap
library(ComplexHeatmap)

### colume metadata
meta_data <- data.frame(
  `Technique Platform` = c(rep("scRNA-seq", 12), rep("ST technology", 11))
)
rownames(meta_data) <- colnames(platform_data)[3:ncol(platform_data)]

### row metadata
row_meta <- data.frame(
  `class` = c(rep("class1", 10),
              rep("class2", 15),
              rep("class3", 10),
              rep("class4", 8),
              rep("class5", 6)),
  `color` = c(rep(method_class_colors[1], 10),
              rep(method_class_colors[2], 15),
              rep(method_class_colors[3], 10),
              rep(method_class_colors[4], 8),
              rep(method_class_colors[5], 6))
)
rownames(row_meta) <- rownames(platform_data)

col <- list(Technique.Platform = c(`scRNA-seq` = technique_colors[1],
                                   `ST technology` = technique_colors[2]),
            class = c(class1 = method_class_colors[1],
                      class2 = method_class_colors[2],
                      class3 = method_class_colors[3],
                      class4 = method_class_colors[4],
                      class5 = method_class_colors[5]))

pdf(file = "../sim-article/figures/Fig4c.pdf", width = 7, height = 7)
p4_c <- ComplexHeatmap::pheatmap(platform_data %>% column_to_rownames("id") %>% select(-1) %>% as.matrix(),
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
                                 gaps_row = c(10, 25, 35, 43), 
                                 color = c(colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(40)),
                                 na_col = "gray80",
                                 right_annotation = rowAnnotation(bar1 = anno_barplot(overall_data$`scRNA-seq data`,
                                                                                      gp = gpar(fill = row_meta$color)),
                                                                  bar2 = anno_barplot(overall_data$`spatial transcriptome data`,
                                                                                      gp = gpar(fill = row_meta$color))))
print(p4_c)
dev.off()




###--------------------------------------------------------------------------###
###                                    Fig4a
###--------------------------------------------------------------------------###
fig4a <- overall_data %>% 
  select(c(1,2, 15:22)) %>% 
  pivot_longer(cols = 3:ncol(.), names_to = "metric", values_to = "value")
fig4a$Category <- factor(fig4a$Category, levels = paste0("Class ", 1:5))
fig4a$metric <- factor(fig4a$metric, levels = c("MAD", "MAE", "RMSE", "KS", "OV", "bhattacharyya", "KDE", "multiKS"))
p4_a <- ggboxplot(data = fig4a,
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

ggsave(plot = p4_a,
       filename = "../sim-article/figures/Fig4a.pdf",
       width = 8,
       height = 3.5,
       units = "cm")


###--------------------------------------------------------------------------###
###                                    Fig4b
###--------------------------------------------------------------------------###
fig4b <- overall_data %>% 
  select(c(1,2, 23:37)) %>% 
  pivot_longer(cols = 3:ncol(.), names_to = "property", values_to = "value")
fig4b$Category <- factor(fig4b$Category, levels = paste0("Class ", 1:5))
fig4b$property <- factor(fig4b$property)
p4_b <- ggboxplot(data = fig4b,
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

ggsave(plot = p4_a / p4_b,
       filename = "../sim-article/figures/Fig4b.pdf",
       width = 8,
       height = 9,
       units = "cm")

ggsave(plot = p4_a / p4_b + plot_layout(height = c(1,2)),
       filename = "../sim-article/figures/Fig4ab.pdf",
       width = 9,
       height = 12.5,
       units = "cm")

#### Fig4-S
colors <- RColorBrewer::brewer.pal(12, "Set3")[4:8]
methods <- openxlsx::read.xlsx("./Chunk1-Data preparation/methods.xlsx")

platform_p4s_data <- platform_data %>% 
  pivot_longer(cols = 3:ncol(.), names_to = "Platform", values_to = "platform_score") %>% 
  rename(., Method = id)

platform_p4s_data$Method <- factor(platform_p4s_data$Method, levels = overall_data$id)
### top 5 method be labelled
top_5_methods <- platform_p4s_data %>% 
  group_by(Platform) %>% 
  top_n(., 5, platform_score) %>% 
  arrange(desc(platform_score)) %>% 
  mutate(
    n = 1:n()
  ) %>% 
  ungroup() %>% 
  mutate(
    n = case_when(
      n == 6 ~ 5,
      TRUE ~ n
    )
  )

### bottom 5 method be labelled
bottom_5_methods <- platform_p4s_data %>% 
  mutate(
    platform_score = 1 - platform_score
  ) %>% 
  group_by(Platform) %>% 
  top_n(., 5, platform_score) %>% 
  arrange(desc(platform_score)) %>% 
  mutate(
    n = 1:n()
  ) %>% 
  ungroup() %>% 
  mutate(
    n = case_when(
      n == 6 ~ 5,
      TRUE ~ n
    ),
    n = -n
  )


P4_s <- ggbarplot(data = platform_p4s_data,
                  x = "Method",
                  y = "platform_score",
                  color = "Category",
                  fill = "Category",
                  ggtheme = theme_pubr(),
                  size = 1,
                  alpha = 0.6)+
  geom_text(data = top_5_methods,
            mapping = aes(x = Method, y = platform_score, label = n), nudge_y = -0.29) +
  geom_text(data = bottom_5_methods,
            mapping = aes(x = Method, y = platform_score, label = n), nudge_y = -0.4) +
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
        strip.text.y = element_text(size = 8, angle = 0))

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
odds <- seq(1, 49, 1)
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
library(ggpmisc)
fig4d_data <- overall_data %>% 
  select(all_of(c("id", "scRNA-seq data", "spatial transcriptome data", "Category"))) %>% 
  rename(., `ST technology` = `spatial transcriptome data`) %>% 
  rename(., `scRNA-seq` = `scRNA-seq data`) %>% 
  rename(., Method = id) %>% 
  drop_na()

cor_value <- sprintf("%0.2f", cor(fig4d_data$`ST technology`, fig4d_data$`scRNA-seq`))

fig4d <- ggplot(fig4d_data, aes(x = `ST technology`, y = `scRNA-seq`)) +
  geom_smooth(method = 'lm', linewidth = 1, color = "red") +
  geom_point(size = 3) +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p.value.label"), sep = "*\", \"*"),
               formula = y ~ x,
               parse = TRUE,
               size = 4.5,
               label.x = 0.1,
               label.y = 0.96) +
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
  theme(legend.position = c(0.1, 0.7),
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


