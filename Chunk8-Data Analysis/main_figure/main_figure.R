library(stringr)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(funkyheatmap)


################################################################################
#######################   Figure 2 --- detail heatmap  #########################
################################################################################
accuracy <- readRDS("Chunk8-Data Analysis/accuracy/accuracy_long_data.rds")
functionality <- readRDS("Chunk8-Data Analysis/functionality/functionality_data.rds")
usability <- readRDS("Chunk8-Data Analysis/usability/usability.rds")

### summarize accuracy score per metric
accuracy_summary_per_metric <- accuracy %>% 
  group_by(Method, metric) %>% 
  summarise(
    value = mean(value, na.rm = TRUE)
  ) %>%
  pivot_wider(., names_from = metric, values_from = value)
  
### summarize accuracy score
accuracy_score <- accuracy %>% 
  group_by(Method, metric) %>% 
  summarise(
    value = mean(value, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  group_by(Method) %>% 
  summarise(
    accuracy = mean(value, na.rm = TRUE)
  )

### add accuracy score to the metric result
accuracy_plot <- accuracy_score %>% 
  full_join(., accuracy_summary_per_metric, by = "Method") %>% 
  arrange(desc(accuracy))



### summarize functionality score per metric
functionality_summary_per_metric <- functionality %>% 
  pivot_longer(cols = 3:ncol(.), names_to = "metric", values_to = "value") %>% 
  group_by(Method, metric) %>% 
  summarise(
    value = mean(value, na.rm = TRUE)
  ) %>%
  pivot_wider(., names_from = metric, values_from = value)

### summarize functionality score
functionality_score <- functionality %>% 
  pivot_longer(cols = 3:ncol(.), names_to = "metric", values_to = "value") %>% 
  group_by(Method, metric) %>% 
  summarise(
    value = mean(value, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  group_by(Method) %>% 
  summarise(
    functionality = mean(value, na.rm = TRUE)
  )

### add functionality score to the metric result
functionality_plot <- functionality_score %>% 
  full_join(., functionality_summary_per_metric, by = "Method") %>% 
  arrange(desc(functionality))


colnames(usability) <- c("Method", "usability")

overall_data <- accuracy_plot %>% 
  inner_join(., functionality_plot, by = "Method") %>% 
  inner_join(., usability, by = "Method") %>% 
  relocate(accuracy, .after = Method) %>% 
  relocate(functionality, .after = accuracy) %>% 
  relocate(usability, .after = functionality)

overall_data <- overall_data%>% 
  mutate(
    overall = (accuracy + functionality + usability)/3
  ) %>% 
  relocate(overall, .after = Method) %>% 
  arrange(desc(overall))


overall_data <- overall_data%>% 
  mutate(
    id = Method,
  ) %>% 
  relocate(id, .before = Method) %>% 
  select(-2)

### define column information

column_info <- column_info <- tribble(
  ~ id,                  ~group,             ~name,                     ~geom,       ~palette,
  "id",                  "method",           "",                        "text",      NA,
  "overall",             "overall",          "Overall",                 "bar",       "palette1",
  "accuracy",            "overall",          "Accuracy",                "bar",       "palette2",
  "functionality",       "overall",          "Functionality",           "bar",       "palette3",
  "usability",           "overall",          "Usability",               "bar",       "palette4",
  "bhattacharyya",       "property",         "bhattacharyya distance",  "funkyrect", "palette2",
  "KDE",                 "property",         "KDE",                     "funkyrect", "palette2",
  "KS",                  "property",         "KS",                      "funkyrect", "palette2",
  "MAD",                 "property",         "MAD",                     "funkyrect", "palette2",
  "MAE",                 "property",         "MAE",                     "funkyrect", "palette2",
  "multiKS",             "property",         "multiKS",                 "funkyrect", "palette2",
  "OV",                  "property",         "Overlapping",             "funkyrect", "palette2",
  "RMSE",                "property",         "RMSE",                    "funkyrect", "palette2",
  "Accuracy",            "property",         "Accuracy",                "funkyrect", "palette2",
  "AWS_batch",           "group",            "AWS_batch",               "funkyrect", "palette3",
  "CDI",                 "group",            "CDI",                     "funkyrect", "palette3",
  "cms",                 "group",            "cms",                     "funkyrect", "palette3",
  "connectivity",        "group",            "Connectivity",            "funkyrect", "palette3",
  "DB_index",            "group",            "DB Index",                "funkyrect", "palette3",
  "distribution_score",  "group",            "Distributon Score",       "funkyrect", "palette3",
  "dunn",                "batch",            "Dunn Index",              "funkyrect", "palette4",
  "F1",                  "batch",            "F1 Score",                "funkyrect", "palette4",
  "kBET",                "batch",            "kBET",                    "funkyrect", "palette4",
  "LISI",                "batch",            "LISI",                    "funkyrect", "palette4",
  "mm",                  "batch",            "mm",                      "funkyrect", "palette4",
  "pcr",                 "batch",            "PCR",                     "funkyrect", "palette4",
  "Precision",           "DEGs",             "Precision",               "funkyrect", "palette1",
  "Recall",              "DEGs",             "Recall",                  "funkyrect", "palette1",
  "ROUGE",               "DEGs",             "Recall",                  "funkyrect", "palette1",
  "shannon_entropy",     "DEGs",             "Shannon Entropy",         "funkyrect", "palette1",
  "silhouette",          "DEGs",             "Silhouette",              "funkyrect", "palette1",
  "true_proportion",     "DEGs",             "True Proportion",         "funkyrect", "palette1"
)

### column grouping

column_groups <- tribble(
  ~ Experiment,     ~Category,                            ~group,         ~palette,
  "Method",         "",                                   "method",       "palette1",
  "Summary",        "Aggregated scores per term",         "overall",      "palette1",
  "Accuracy",       "Per metric",                         "property",     "palette2",
  "Functionality",  "Group",                              "group",        "palette3",
  "Functionality",  "Batch",                              "batch",        "palette3",
  "Functionality",  "DEGs",                               "DEGs",         "palette3"
)


### determine palettes

palettes <- tribble(
  ~palette,             ~colours,
  "palette1",           grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Greys")[-1]))(101),
  "palette2",           grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues") %>% c("#011636")))(101),
  "palette3",           grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Reds")[-8:-9]))(101),
  "palette4",           grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "YlOrBr")[-8:-9]))(101)
)


detailed_plot <- funky_heatmap(data = overall_data,
                               column_info = column_info,
                               column_groups = column_groups,
                               palettes = palettes,
                               scale_column = FALSE,
                               expand = c(xmin = 1, xmax = 4, ymin = 0, ymax = 0))

ggsave(detailed_plot, filename = "Chunk8-Data Analysis/main_figure/detailed_plot.pdf", units = "in", width = 20, height = 10)





