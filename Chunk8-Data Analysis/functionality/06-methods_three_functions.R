source("./Chunk8-Data Analysis/functionality/utils_functions.R")

overall_data <- readRDS("./Chunk8-Data Analysis/overall_data.rds")

function_data <- overall_data %>% 
  select(c(1:3, 69:121)) %>% 
  slice(1:43)
colnames(function_data)[32:56] <- str_split(colnames(function_data)[32:56], "_", simplify = TRUE)[, 2]

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

####################### Platform
color <- c("#F1947C", "#6CB9B0")
method_class_colors <- c("#F97369", "#69ACD1", "#F9A24B", "#ADDD52", "#F9B4DA")

per_platform_score <- function_data[, c("id", platforms)] %>% 
  column_to_rownames("id")

#### heatmap
library(ComplexHeatmap)

### colume metadata
meta_data <- data.frame(
  `Technique Platform` = c(rep("scRNA-seq", 12), rep("ST technology", 11))
)
rownames(meta_data) <- colnames(per_platform_score)

### row metadata
row_meta <- data.frame(
  `class` = c(rep("class1", 10),
              rep("class2", 15),
              rep("class3", 10),
              rep("class4", 8)),
  `color` = c(rep(method_class_colors[1], 10),
              rep(method_class_colors[2], 15),
              rep(method_class_colors[3], 10),
              rep(method_class_colors[4], 8))
)
rownames(row_meta) <- rownames(per_platform_score)

col <- list(Technique.Platform = c(`scRNA-seq` = color[1],
                                   `ST technology` = color[2]),
            class = c(class1 = method_class_colors[1],
                      class2 = method_class_colors[2],
                      class3 = method_class_colors[3],
                      class4 = method_class_colors[4]))

pdf(file = "../sim-article/figures/Fig5-a1.pdf", width = 10, height = 8)
p5_a <- ComplexHeatmap::pheatmap(per_platform_score %>% as.matrix(),
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
                                 right_annotation = rowAnnotation(bar1 = anno_barplot(function_data$`scRNA-seq data`[1:43],
                                                                                      gp = gpar(fill = row_meta$color)),
                                                                  bar2 = anno_barplot(function_data$`spatial transcriptome data`[1:43],
                                                                                      gp = gpar(fill = row_meta$color)))
                                 )
print(p5_a)
dev.off()



#### Fig5b boxplot
methods <- openxlsx::read.xlsx("./Chunk1-Data preparation/methods.xlsx")

### normalize every metrics to [0, 1]
platform_p5s_data <- function_data %>% 
  rename(., Method = id) %>% 
  select(c(1:2, 34:56)) %>% 
  pivot_longer(cols = 3:25, names_to = "Platform", values_to = "platform_score")
  
platform_p5s_data$Method <- factor(platform_p5s_data$Method, levels = overall_data$id)
platform_p5s_data$Category <- factor(platform_p5s_data$Category, levels = paste0("Class ", 1:5))
P5_b <- ggboxplot(data = platform_p5s_data,
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
  scale_color_manual(values = method_class_colors) +
  scale_fill_manual(values = colors) +
  facet_wrap(.~ Platform, ncol = 4, strip.position = "top") +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10),
        legend.position = "top",
        legend.background = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text.x = element_text(size = 8))

ggsave(plot = P5_b, filename = "../sim-article/figures/Fig5-b2.pdf", width = 6, height = 8, units = "in")
