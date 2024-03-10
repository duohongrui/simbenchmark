library(tidyverse)
################################################################################
###############################    Figure 4b  ##################################
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

platform_score_per_model <- overall_data[, c("id", "Model Category", "accuracy", platforms)]
platform_score_per_model <- platform_score_per_model %>% 
  group_by(`Model Category`) %>% 
  summarise(
    across(all_of(colnames(platform_score_per_model)[3:27]), ~ mean(.x, na.rm = TRUE))
  ) %>% 
  select(-2) %>% 
  tibble::column_to_rownames("Model Category")

technique_colors <- c("#7DAEF4", "#E0709E")
method_class_colors <- c(RColorBrewer::brewer.pal(12, "Set3")[4],
                         "#94D3C7","#AAAAE0","#ADCE65","#EAC5E3")
meta_data <- data.frame(
  `Technique Platform` = c(rep("scRNA-seq", 13), rep("ST technology", 11))
)
rownames(meta_data) <- colnames(platform_score_per_model)
col <- list(Technique.Platform = c(`scRNA-seq` = technique_colors[1],
                                   `ST technology` = technique_colors[2]))


pdf(file = "../sim-article/figures/Fig4b.pdf", width = 7, height = 3.2)
p4_b <- ComplexHeatmap::pheatmap(platform_score_per_model %>% as.matrix(),
                                 show_colnames = TRUE,
                                 show_rownames = TRUE,
                                 scale = "none",
                                 cluster_cols = FALSE,
                                 cluster_rows = TRUE,
                                 annotation_col = meta_data,
                                 annotation_colors = col,
                                 treeheight_row = 20,
                                 cutree_rows = 4,
                                 border_color = "white",
                                 color = c(colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(40)),
                                 na_col = "gray80")
print(p4_b)
dev.off()

################################################################################
###############################    Figure 4c  ##################################
################################################################################
accuracy_per_model <- overall_data[, c("id", "Model Category", "accuracy")]

ggplot(accuracy_per_model, aes(fill = `Model Category`,
                               x = reorder(`Model Category`, accuracy, decreasing = TRUE),
                               y = accuracy))+
  geom_bar(stat = "summary", position = position_dodge(1), color = "black", alpha = 0.4) +
  stat_summary(fun.data = 'mean_se',
               geom = "errorbar",
               colour = "black",
               width = 0.2,
               position = position_dodge(1)) +
  geom_point(size = 3.5,
             color = "black",
             shape = 21,
             stroke = 0.01,
             show.legend = FALSE) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(name = "Paired", n = 12)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  ylab("Accuracy scores") +
  xlab("Model Category")
ggsave(filename = "../sim-article/figures/Fig4c_revised.pdf")




################################################################################
###############################    Figure 4d  ##################################
################################################################################
### read data
accuracy_data <- readRDS("./Chunk8-Data Analysis/2-accuracy/accuracy_long_data.rds")
method <- openxlsx::read.xlsx("Chunk1-Data preparation/methods.xlsx", sheet = 1, colNames = TRUE, sep.names = " ") %>% 
  select(-2)
### exclude MFA as it do not use zero inflated model by default
accuracy_data <- accuracy_data %>% 
  left_join(method %>% select(Method, `Model Category`), by = "Method") %>% 
  filter(Method != "MFA")

# accuracy_data <- accuracy_data %>% 
#   left_join(method %>% select(Method, `Model Category`), by = "Method")

### property
accuracy_summary_per_property <- accuracy_data %>% 
  group_by(Method, property, Data, `Model Category`) %>% 
  summarise(
    value = mean(value, na.rm = TRUE)
  ) %>% 
  group_by(Method, property, `Model Category`) %>% 
  summarise(
    value = mean(value, na.rm = TRUE)
  ) %>% 
  group_by(property, `Model Category`) %>% 
  summarise(
    value = mean(value, na.rm = TRUE)
  )

property <- accuracy_summary_per_property %>% 
  pivot_wider(names_from = "property", values_from = "value") %>% 
  column_to_rownames(var = "Model Category")

pdf(file = "../sim-article/figures/Fig4c.pdf", width = 5, height = 6.3)
ComplexHeatmap::pheatmap(property %>% as.matrix() %>% t(),
                         show_colnames = TRUE,
                         show_rownames = TRUE,
                         scale = "none",
                         cluster_cols = TRUE,
                         cluster_rows = FALSE,
                         cutree_cols = 4,
                         # annotation_col = meta_data,
                         # annotation_colors = col,
                         # annotation_row = row_meta %>% select(1),
                         treeheight_col = 20,
                         border_color = "white",
                         color = c(colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(40)),
                         na_col = "gray80")
dev.off()



### summarize accuracy scores for each model class
accuracy_summary_per_metric <- accuracy_data %>% 
  group_by(Method, metric, `Model Category`, Data) %>% 
  summarise(
    value = mean(value, na.rm = TRUE),
  ) %>% 
  group_by(Method, metric, `Model Category`) %>% 
  summarise(
    value = mean(value, na.rm = TRUE)
  ) %>% 
  group_by(metric, `Model Category`) %>% 
  summarise(
    value = mean(value, na.rm = TRUE)
  )

model_accuracy <- accuracy_summary_per_metric %>% 
  pivot_wider(names_from = "metric", values_from = "value") %>% 
  column_to_rownames(var = "Model Category")
library(ComplexHeatmap)
pdf(file = "../sim-article/figures/Fig4d.pdf", width = 6, height = 5)
ComplexHeatmap::pheatmap(model_accuracy %>% as.matrix() %>% t(),
                         show_colnames = TRUE,
                         show_rownames = TRUE,
                         scale = "none",
                         cluster_cols = TRUE,
                         cluster_rows = FALSE,
                         cutree_cols = 4,
                         # annotation_col = meta_data,
                         # annotation_colors = col,
                         # annotation_row = row_meta %>% select(1),
                         treeheight_col = 20,
                         border_color = "white",
                         color = c(colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(40)),
                         na_col = "gray80")
dev.off()
