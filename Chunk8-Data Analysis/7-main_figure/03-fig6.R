library(ggpubr)
source("./Chunk8-Data Analysis/3-functionality/07-utils_functions.R")

overall_data <- readRDS("./Chunk8-Data Analysis/overall_data.rds")

function_data <- overall_data %>% 
  select(c(1:3, 14, 74:129)) %>% 
  slice(1:43)
colnames(function_data)[33:60] <- str_split(colnames(function_data)[33:60], "_", simplify = TRUE)[, 2]

data_info <- openxlsx::read.xlsx("./Chunk1-Data preparation/evaluation_datasets.xlsx", rowNames = TRUE) %>% 
  mutate(
    Data = str_split(Dataset.ID, "_", simplify = TRUE)[, 1]
  ) %>% 
  select(Data, Platform)
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

technique_colors <- c("#7DAEF4", "#E0709E")
method_class_colors <- c(RColorBrewer::brewer.pal(12, "Set3")[4],
                         "#94D3C7","#AAAAE0","#ADCE65","#EAC5E3")
###--------------------------------------------------------------------------###
###                                    Fig6a
###--------------------------------------------------------------------------###
functionality <- c("Group", "DEGs", "Batch", "Trajectory")

function_bar_plot <- map(functionality, .f = function(x){
  sub_data <- function_data %>% 
    select(all_of(c("id", "Category", paste0(x, "_score")))) %>% 
    rename(., x = paste0(x, "_score")) %>% 
    drop_na() %>% 
    arrange(desc(x)) %>% 
    mutate(
      color = case_when(
        Category == "Class 1" ~ method_class_colors[1],
        Category == "Class 2" ~ method_class_colors[2],
        Category == "Class 3" ~ method_class_colors[3],
        Category == "Class 4" ~ method_class_colors[4],
      )
    )
  sub_data$id <- factor(sub_data$id, levels = rev(sub_data$id))
  sub_data$Category <- factor(sub_data$Category, levels = paste0("Class ", 1:4))
  if(x == "Trajectory"){
    color <- method_class_colors[1]
  }else if(x == "Batch"){
    color <- method_class_colors[c(2,3,4)]
  }else{
    color <- method_class_colors[-5]
  }
  dotplot <- ggdotchart(data = sub_data,
                        x = "id",
                        y = "x",
                        color = "Category",
                        dot.size = 1,
                        shape = 15,
                        size = 0.3,
                        add = "segments") +
    geom_text(aes(label = sprintf("%0.2f", x)),
              nudge_y = 0.14,
              size = 1.1) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_line(linewidth = 0.1, lineend = 0.04),
      axis.line = element_line(colour = "black", linewidth = 0.1),
      axis.title = element_blank(),
      axis.text = element_text(color = "black", size = 3.2)
    ) +
    coord_flip() +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, 1.1),
                       breaks = seq(0,1,0.2),
                       labels = as.character(seq(0,1,0.2))) +
    scale_color_manual(values = color)
  
  if(x == "Group"){
    dotplot <- dotplot +
      theme(
        legend.position = c(0.9, 0.1),
        title = element_text(size = 4),
        legend.title = element_text(size = 3),
        legend.text = element_text(size = 3),
        legend.key.size = unit(0.1, "cm"),
        legend.direction = "vertical",
        legend.background = element_blank()
      )
  }else{
    dotplot <- dotplot +
      theme(
        legend.position = "none"
      )
  }
  
  dotplot
})

ggsave(plot = (function_bar_plot[[1]] | function_bar_plot[[2]]) / (function_bar_plot[[3]] | function_bar_plot[[4]]) + plot_layout(heights = c(2,1)),
       filename = "../sim-article/figures/Fig6a.pdf",
       width = 7,
       height = 9,
       units = "cm")

###--------------------------------------------------------------------------###
###                        Fig6b and Supplementary Fig 8
###--------------------------------------------------------------------------###
library(ggrepel)
library(ggpmisc)

functionality_scores <- c("Group_score", "DEGs_score", "Batch_score", "Trajectory_score")
combn_func <- combn(functionality_scores,2)

fig6b <- map(1:4, function(x){
  x1 <- combn_func[, x][1]
  x2 <- combn_func[, x][2]
  data <- function_data %>% 
    select(all_of(c("id", x1, x2, "Category"))) %>% 
    rename(., Method = id) %>% 
    drop_na()
  ggplot(data, aes(x = get(x1), y = get(x2))) +
    geom_smooth(method = 'lm', linewidth = 0.6, color = "red") +
    geom_point(size = 0.8) +
    stat_cor(method = "pearson", size = 2) +
    geom_label_repel(aes(label = Method,
                         fill = Category),
                     color = "black",
                     size = 1,
                     box.padding = unit(0.2, "lines"),
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
          legend.key.size = unit(0.2, 'cm'),
          axis.ticks = element_line(color = "black", linewidth = 0.1, lineend = 0.3),
          axis.line = element_line(color = "black", linewidth = 0.08)) +
    xlab(paste0("Scores of ", x1, " functionality")) +
    ylab(paste0("Scores of ", x2, " functionality")) +
    scale_fill_manual(values = method_class_colors)
})

ggsave(plot = patchwork::wrap_plots(fig6b, ncol = 2),
       filename = "../sim-article/figures/Fig6b.pdf",
       width = 14,
       height = 14,
       units = "cm")



###--------------------------------------------------------------------------###
###                                    Fig6c
###--------------------------------------------------------------------------###
arrange_by_group <- function(tibble){
  groups <- unique(tibble %>% pull("Category"))
  result <- tibble()
  for(i in groups){
    tmp <- tibble %>% 
      filter(Category == i) %>% 
      arrange(desc(functionality))
    result <- rbind(result, tmp)
  }
  return(result)
}
per_platform_score <- function_data[, c("id", "Category", "functionality", platforms)]
per_platform_score <- arrange_by_group(per_platform_score) %>% 
  select(-c(2,3,14,15,26,27)) %>% 
  tibble::column_to_rownames("id")

#### heatmap
library(ComplexHeatmap)

### colume metadata
meta_data <- data.frame(
  `Technique Platform` = c(rep("scRNA-seq", 11), rep("ST technology", 9))
)
rownames(meta_data) <- colnames(per_platform_score)

### row metadata
row_meta <- data.frame(
  `class` = c(rep("class1", 15),
              rep("class2", 10),
              rep("class3", 9),
              rep("class4", 9)),
  `color` = c(rep(method_class_colors[1], 15),
              rep(method_class_colors[2], 10),
              rep(method_class_colors[3], 9),
              rep(method_class_colors[4], 9))
)
rownames(row_meta) <- rownames(per_platform_score)

col <- list(Technique.Platform = c(`scRNA-seq` = technique_colors[1],
                                   `ST technology` = technique_colors[2]),
            class = c(class1 = method_class_colors[1],
                      class2 = method_class_colors[2],
                      class3 = method_class_colors[3],
                      class4 = method_class_colors[4]))

pdf(file = "../sim-article/figures/fig6c.pdf", width = 7.2, height = 8.5)
p6_c <- ComplexHeatmap::pheatmap(per_platform_score %>% as.matrix(),
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
                                 gaps_row = c(15, 25, 34), 
                                 color = c(colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(40)),
                                 na_col = "grey80")
print(p6_c)
dev.off()



###--------------------------------------------------------------------------###
###                           Supplementary Fig 11
###--------------------------------------------------------------------------###
methods <- openxlsx::read.xlsx("./Chunk1-Data preparation/methods.xlsx")

### normalize every metrics to [0, 1]
platform_supp11_data <- function_data %>% 
  rename(., Method = id) %>% 
  select(c(1:2, 35:58)) %>% 
  pivot_longer(cols = 3:26, names_to = "Platform", values_to = "platform_score") %>% 
  filter(Platform != "MERFISH", Platform != "Microwell-seq", Platform != "Stereo-seq", Platform != "STRT-seq",)

platform_supp11_data$Method <- factor(platform_supp11_data$Method, levels = overall_data$id)
platform_supp11_data$Category <- factor(platform_supp11_data$Category, levels = paste0("Class ", 1:5))
platform_supp11_data$Platform <- factor(platform_supp11_data$Platform, levels = platforms)
supp_fig_11 <- ggboxplot(data = platform_supp11_data,
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
  stat_n_text(size = 1, vjust = 0, color = "black") +
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
        legend.position = c(0.5, 1.05),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.title = element_text(size = 4),
        legend.key.size = unit(0.5, 'cm'),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 4),
        strip.text.x = element_text(margin = margin(-0.01,0,-0.01,0, "cm")),
        plot.margin = unit(c(0.6,0.1,0.1,0.1), 'cm')) +
  ylab("Platform scores")

ggsave(plot = supp_fig_11,
       filename = "../sim-article/figures/Supp_Fig_11.pdf",
       width = 6,
       height = 9,
       units = "cm")

###--------------------------------------------------------------------------###
###                         Supplementary Fig 10 (Moran'C statistics)
###--------------------------------------------------------------------------###
data_list <- list.files("../spatial_autocorrelation/")
moransC_data <- map_dfr(data_list, function(x){
  tmp_data <- readRDS(file.path("../spatial_autocorrelation", x))
  tmp_data
})
spatial_auto <- moransC_data %>% 
  group_by(method, data) %>% 
  summarise(
    Moran = mean(Moran, na.rm = TRUE)
  ) %>% 
  ungroup()
spatial_auto <- overall_data %>% 
  select(id, Category) %>% 
  rename(method = id) %>% 
  right_join(spatial_auto, by = "method")

spatial_auto_rank <- spatial_auto %>% 
  group_by(method) %>% 
  summarise(
    Moran = mean(Moran, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  arrange(Moran) %>% 
  pull(method)

spatial_auto$method <- factor(spatial_auto$method, levels = spatial_auto_rank)
ggplot(spatial_auto, aes(x = method, y = Moran))+
  geom_boxplot(aes(color = Category), outlier.size = 0, width = 0.5) +
  geom_point(aes(color = Category), size = 0.8) +
  scale_color_manual(values = method_class_colors) +
  coord_flip() +
  theme_bw() +
  theme(
    axis.title = element_text(color = "black"),
    axis.title.y = element_blank(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "top"
  ) +
  ylab("Moran's C statistics")

ggsave(filename = "../sim-article/figures/Supp_Fig_10.pdf",
       width = 6,
       height = 6)


###--------------------------------------------------------------------------###
###                                    Fig6d
###--------------------------------------------------------------------------###
library(ggrepel)
library(ggpmisc)
fig6d_data <- function_data %>% 
  select(all_of(c("id", "scRNA-seq data", "spatial transcriptome data", "Category"))) %>% 
  rename(., `ST technology` = `spatial transcriptome data`) %>% 
  rename(., `scRNA-seq` = `scRNA-seq data`) %>% 
  rename(., Method = id) %>% 
  drop_na()

fig6d <- ggplot(fig6d_data, aes(x = `ST technology`, y = `scRNA-seq`)) +
  geom_smooth(method = 'lm', linewidth = 0.6, color = "red") +
  geom_point(size = 0.8) +
  stat_cor(method = "pearson", size = 2) +
  geom_label_repel(aes(label = Method,
                       fill = Category),
                   color = "black",
                   size = 1,
                   box.padding = unit(0.2, "lines"),
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
        legend.key.size = unit(0.2, 'cm'),
        axis.ticks = element_line(color = "black", linewidth = 0.1, lineend = 0.3),
        axis.line = element_line(color = "black", linewidth = 0.08)) +
  xlab("Functionality scores on the spatial datasets") +
  ylab("Functionality scores on the scRNA-seq datasets") +
  scale_fill_manual(values = method_class_colors)

ggsave(plot = fig6d,
       filename = "../sim-article/figures/Fig6d.pdf",
       width = 7,
       height = 7,
       units = "cm")



###--------------------------------------------------------------------------###
###                             Supplementary Fig 12a
###--------------------------------------------------------------------------###
source("Chunk8-Data Analysis/3-functionality/07-utils_functions.R")
techniques_data <- data_process(function_data = group_long_data,
                                functionality = "Group",
                                info = "technique")
cor_data <- techniques_data %>% 
  pivot_wider(names_from = "Type", values_from = "value")
cor_data <- cor_data %>% 
  mutate(
    across(all_of(c("scRNA-seq", "ST technology")), ~ replace_na(.x, 0))
  ) %>% 
  filter(`scRNA-seq` != 0, `ST technology` != 0)
supp_fig12a <- ggplot(cor_data, aes(x = `ST technology`, y = `scRNA-seq`)) +
  geom_smooth(method = 'lm', linewidth = 0.6, color = "red") +
  geom_point(size = 0.8) +
  stat_cor(method = "pearson", size = 2) +
  geom_label_repel(aes(label = Method,
                       fill = Category),
                   color = "black",
                   size = 1,
                   box.padding = unit(0.2, "lines"),
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
        legend.key.size = unit(0.2, 'cm'),
        axis.ticks = element_line(color = "black", linewidth = 0.1, lineend = 0.3),
        axis.line = element_line(color = "black", linewidth = 0.08)) +
  xlab("Group scores on the spatial datasets") +
  ylab("Group scores on the scRNA-seq datasets") +
  scale_fill_manual(values = method_class_colors)

ggsave(plot = supp_fig12a,
       filename = "../sim-article/figures/Supp_Fig_12a.pdf",
       width = 7,
       height = 7,
       units = "cm")

###--------------------------------------------------------------------------###
###                             Supp_Fig8+9
###--------------------------------------------------------------------------###
source("./Chunk8-Data Analysis/3-functionality/07-utils_functions.R")
group_long_data <- readRDS("Chunk8-Data Analysis/3-functionality/group_long_data.rds")

group <- function_metric_for_datasets(function_data = group_long_data,
                                      functionality = "Group")

DEGs_long_data <- readRDS("Chunk8-Data Analysis/3-functionality/DEGs_long_data.rds")
DEGs <- function_metric_for_datasets(function_data = DEGs_long_data,
                                     functionality = "DEGs")

batch_long_data <- readRDS("Chunk8-Data Analysis/3-functionality/batch_long_data.rds")
batch <- function_metric_for_datasets(function_data = batch_long_data,
                                      functionality = "Batch")

trajectory_long_data <- readRDS("Chunk8-Data Analysis/3-functionality/trajectory_long_data.rds")
traj <- function_metric_for_datasets(function_data = trajectory_long_data,
                                     functionality = "Trajectory")

###--------------------------------------------------------------------------###
###                             Supp_Fig_12bc
###--------------------------------------------------------------------------###

ggsave(plot = wrap_plots(DEGs$fig + traj$fig) + plot_annotation(tag_levels = "a"),
       filename = "../sim-article/figures/Supp_Fig_12bc.pdf",
       width = 16,
       height = 10,
       units = "cm")

###--------------------------------------------------------------------------###
###                             Supp_Fig_13
###--------------------------------------------------------------------------###

ggsave(plot = wrap_plots(group$technique_bar / DEGs$technique_bar / traj$technique_bar) +
         plot_annotation(tag_levels = "a"),
       filename = "../sim-article/figures/Supp_Fig_13.pdf")
