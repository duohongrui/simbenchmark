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
technique_colors <- c("#7DAEF4", "#E4B083")
method_class_colors <- RColorBrewer::brewer.pal(12, "Set3")[c(4,1,3,7,8)]

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

col <- list(Technique.Platform = c(`scRNA-seq` = technique_colors[1],
                                   `ST technology` = technique_colors[2]),
            class = c(class1 = method_class_colors[1],
                      class2 = method_class_colors[2],
                      class3 = method_class_colors[3],
                      class4 = method_class_colors[4]))

pdf(file = "../sim-article/figures/Fig5a.pdf", width = 7.6, height = 8)
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
                                 gaps_row = c(10, 25, 35, 43), 
                                 color = c(colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(40)),
                                 na_col = "grey80")
print(p5_a)
dev.off()



#### Fig5b boxplot
methods <- openxlsx::read.xlsx("./Chunk1-Data preparation/methods.xlsx")

### normalize every metrics to [0, 1]
platform_p5s_data <- function_data %>% 
  rename(., Method = id) %>% 
  select(c(1:2, 34:56)) %>% 
  pivot_longer(cols = 3:25, names_to = "Platform", values_to = "platform_score") %>% 
  filter(Platform != "MERFISH", Platform != "Microwell-seq", Platform != "Stereo-Seq")
  
platform_p5s_data$Method <- factor(platform_p5s_data$Method, levels = overall_data$id)
platform_p5s_data$Category <- factor(platform_p5s_data$Category, levels = paste0("Class ", 1:5))
platform_p5s_data$Platform <- factor(platform_p5s_data$Platform, levels = platforms)
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
  scale_fill_manual(values = method_class_colors) +
  facet_wrap(.~ Platform, ncol = 4, strip.position = "top") +
  theme(axis.text = element_text(size = 4),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title = element_text(size = 5),
        axis.title.x = element_blank(),
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
  ylab("Platform scores")

ggsave(plot = P5_b,
       filename = "../sim-article/figures/Fig5b.pdf",
       width = 10,
       height = 9,
       units = "cm")
