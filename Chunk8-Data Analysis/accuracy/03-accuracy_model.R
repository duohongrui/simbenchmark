library(tidyverse)

### read data
accuracy_data <- readRDS("./Chunk8-Data Analysis/accuracy/accuracy_long_data.rds")
method <- openxlsx::read.xlsx("Chunk1-Data preparation/methods.xlsx", sheet = 1, colNames = TRUE, sep.names = " ") %>% 
  select(-2)
# accuracy_data <- accuracy_data %>% 
#   left_join(method %>% select(Method, `Model Category`), by = "Method")

### exclude MFA as it do not use zero inflated model by default
accuracy_data <- accuracy_data %>% 
  left_join(method %>% select(Method, `Model Category`), by = "Method") %>% 
  filter(Method != "MFA")

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

pdf(file = "../sim-article/figures/Fig4c", width = 7, height = 4.5)
ComplexHeatmap::pheatmap(property %>% as.matrix(),
                         show_colnames = TRUE,
                         show_rownames = TRUE,
                         scale = "none",
                         cluster_cols = FALSE,
                         cluster_rows = TRUE,
                         cutree_rows = 4,
                         # annotation_col = meta_data,
                         # annotation_colors = col,
                         # annotation_row = row_meta %>% select(1),
                         treeheight_row = 20,
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
pdf(file = "../sim-article/figures/Fig4d", width = 7, height = 4.5)
ComplexHeatmap::pheatmap(model_accuracy %>% as.matrix(),
                         show_colnames = TRUE,
                         show_rownames = TRUE,
                         scale = "none",
                         cluster_cols = FALSE,
                         cluster_rows = TRUE,
                         cutree_rows = 4,
                         # annotation_col = meta_data,
                         # annotation_colors = col,
                         # annotation_row = row_meta %>% select(1),
                         treeheight_row = 20,
                         border_color = "white",
                         color = c(colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(40)),
                         na_col = "gray80")
dev.off()
