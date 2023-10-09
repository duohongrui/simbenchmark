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
  select(1:66)
colnames(overall_data)[39:66] <- str_split(colnames(overall_data)[39:66], "_", simplify = TRUE)[, 2]
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
###                                    Fig4a
###--------------------------------------------------------------------------###
fig4a <- overall_data %>% 
  select(c(1,2, 16:23)) %>% 
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
###                                    Fig4b
###--------------------------------------------------------------------------###
fig4b <- overall_data %>% 
  select(c(1,2, 24:38)) %>% 
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

ggsave(plot = p4_a / p4_b + plot_layout(height = c(1,2)),
       filename = "../sim-article/figures/Fig4ab.pdf",
       width = 10,
       height = 10,
       units = "cm")

###--------------------------------------------------------------------------###
###                         Fig4c,d (03-accuracy_model.R)
###--------------------------------------------------------------------------###


###--------------------------------------------------------------------------###
###                                    Fig4e
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
rownames(row_meta) <- rownames(per_platform_score)

col <- list(Technique.Platform = c(`scRNA-seq` = technique_colors[1],
                                   `ST technology` = technique_colors[2]),
            class = c(class1 = method_class_colors[1],
                      class2 = method_class_colors[2],
                      class3 = method_class_colors[3],
                      class4 = method_class_colors[4],
                      class5 = method_class_colors[5]))

pdf(file = "../sim-article/figures/Fig4e.pdf", width = 7, height = 8)
p4_e <- ComplexHeatmap::pheatmap(per_platform_score %>% as.matrix(),
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
                                 gaps_row = c(10, 25, 34, 43), 
                                 color = c(colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(40)),
                                 na_col = "gray80")
print(p4_e)
dev.off()


###--------------------------------------------------------------------------###
###                                    Fig4f
###--------------------------------------------------------------------------###
### fig4f
library(ggrepel)
library(ggpmisc)
fig4f_data <- overall_data %>% 
  select(all_of(c("id", "scRNA-seq data", "spatial transcriptome data", "Category"))) %>% 
  rename(., `ST technology` = `spatial transcriptome data`) %>% 
  rename(., `scRNA-seq` = `scRNA-seq data`) %>% 
  rename(., Method = id) %>% 
  drop_na()

fig4f <- ggplot(fig4f_data, aes(x = `ST technology`, y = `scRNA-seq`)) +
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

ggsave(plot = fig4f,
       filename = "../sim-article/figures/Fig4f.pdf",
       width = 7,
       height = 7,
       units = "cm")

###--------------------------------------------------------------------------###
###                                    Fig4g
###--------------------------------------------------------------------------###
### fig4f
fig4g_data <- overall_data %>% 
  select(all_of(c("id", "Category", platforms))) %>% 
  pivot_longer(cols = 3:ncol(.), names_to = "Platform", values_to = "platform_score") %>% 
  rename(., Method = id)
fig4g_data$Method <- factor(fig4g_data$Method, levels = overall_data$id)
fig4g_data <- fig4g_data %>% 
  mutate(
    technology = case_when(
      Platform %in% platforms[1:13] ~ "scRNA-seq",
      Platform %in% platforms[14:24] ~ "ST technology"
    )
  )
sample_size <- fig4g_data %>%
  drop_na() %>% 
  group_by(technology, Category) %>%
  summarize(num = n()) %>% 
  mutate(y = 0.1)
fig4g <- ggplot(fig4g_data, aes(x = Category, y = platform_score, color = technology)) +
  geom_boxplot(outlier.alpha = 0,
               width = 1,
               size = 0.1) +
  geom_jitter(data = fig4g_data,
              aes(x = Category, y = platform_score, color = technology),
              pch = 19,
              size = 0.08,
              position = position_jitterdodge(dodge.width = 0.9)) +
  geom_text(data = sample_size, aes(x = Category, y = y, label = num), size = 1) +
  theme_bw() +
  scale_color_manual(values = technique_colors) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 5, color = "black"),
    axis.text.y = element_text(size = 5, color = "black"),
    axis.title.y = element_text(size = 6.5),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5),
    legend.position = "bottom",
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.box = "horizontal",
    legend.key = element_rect(fill = NA),
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks = element_line(linewidth = 0.3, lineend = 0.5)
  ) +
  stat_compare_means(data = fig4g_data,
                     mapping = aes(group = technology),
                     method = "wilcox.test",
                     label = "p.signif",
                     hide.ns = TRUE)+
  ylab("Accuracy score")

ggsave(plot = fig4g,
       filename = "../sim-article/figures/Fig4g.pdf",
       width = 6,
       height = 8,
       units = "cm")


###--------------------------------------------------------------------------###
###                                    Fig4h
###--------------------------------------------------------------------------###
#### Fig4d boxplot
fig4h_data <- fig4g_data
odds <- seq(1, 49, 2)
rect <- fig4h_data[odds, ]
sample_size <- fig4h_data %>%
  drop_na() %>% 
  group_by(Method, technology) %>%
  summarize(num = n()) %>% 
  mutate(y = 0.1)
fig4h <- ggplot(fig4h_data, aes(x = Method, y = platform_score, color = technology))+
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
  geom_text(data = sample_size, aes(x = Method, y = y, label = num), size = 1) +
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

ggsave(plot = fig4h,
       filename = "../sim-article/figures/Fig4h.pdf",
       width = 18,
       height = 5,
       units = "cm")



###--------------------------------------------------------------------------###
###                             Supplementary Fig. 5
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

supp5_data <- accuracy %>% 
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

supp5_data$Method <- factor(supp5_data$Method, levels = overall_data$id)
odds <- seq(1, 47, 2)
rect <- supp5_data[odds, ]

sample_size <- supp5_data %>%
  drop_na() %>% 
  group_by(Method, `Quantification Strategy`) %>%
  summarize(num = n()) %>% 
  mutate(y = 0)

supp5 <- ggplot(supp5_data, aes(x = Method, y = value, color = `Quantification Strategy`))+
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

ggsave(plot = supp5,
       filename = "../sim-article/figures/Supp_Fig_5.pdf",
       width = 18,
       height = 8,
       units = "cm")
