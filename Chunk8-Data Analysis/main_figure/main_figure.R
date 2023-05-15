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
scalability <- readRDS("Chunk8-Data Analysis/scalability/score_data.rds")
usability <- readRDS("Chunk8-Data Analysis/usability/usability_long_data.rds")
method <- openxlsx::read.xlsx("Chunk1-Data preparation/methods.xlsx", sheet = 1, colNames = TRUE)
colnames(method) <- c("Method",
                      "Category",
                      "Platform",
                      "Prior Information",
                      "Simulate Groups",
                      "Simulate DEGs",
                      "Simulate Batches",
                      "Simulate Trajectory")
method_names <- method$Method

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

functionality_summary_per_metric <- functionality_summary_per_metric[, c("Method", colnames(functionality[c(-1, -2)]))]

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


### scalability score
colnames(scalability) <- c("Method",
                           "estimation time",
                           "estimation memory",
                           "simulation time",
                           "simulation memory",
                           "time",
                           "memory",
                           "scalability")

### usability
usability <- usability %>% 
  group_by(method, category) %>% 
  summarise(
    score = mean(score, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "category", values_from = score)

usability <- usability %>% 
  mutate(
    usability = apply(.[-1], 1, mean)
  )

colnames(usability)[1] <- "Method"


### overall data
overall_data <- method %>% 
  inner_join(., accuracy_plot, by = "Method") %>% 
  inner_join(., functionality_plot, by = "Method") %>% 
  inner_join(., scalability, by = "Method") %>% 
  inner_join(., usability, by = "Method") %>% 
  relocate("accuracy", .after = `Simulate Trajectory`) %>% 
  relocate("functionality", .after = accuracy) %>% 
  relocate("scalability", .after = functionality) %>% 
  relocate("usability", .after = scalability)

overall_data <- overall_data%>% 
  mutate(
    overall = apply(.[, -c(1:8)], MARGIN = 1, function(x){mean(c(x[1], x[2], x[3], x[4]), na.rm = TRUE)})
  ) %>% 
  relocate(overall, .after = `Simulate Trajectory`)


overall_data <- overall_data%>% 
  mutate(
    id = Method,
  ) %>% 
  relocate(id, .before = Method) %>% 
  select(-2)

overall_data <- overall_data %>% 
  mutate(
    across(all_of(c("Prior Information",
                    "Simulate Groups",
                    "Simulate DEGs",
                    "Simulate Batches",
                    "Simulate Trajectory")), ~ replace_na(.x, ""))
  ) %>% 
  as_tibble()

### arrange method

arrange_by_group <- function(tibble){
  groups <- unique(tibble %>% pull("Category"))
  result <- tibble()
  for(i in groups){
    tmp <- tibble %>% 
      filter(Category == i) %>% 
      arrange(desc(overall))
    result <- rbind(result, tmp)
  }
  return(result)
}

overall_data <- arrange_by_group(tibble = overall_data)

###--------------------------------------------------------------------------###
###                            Method Summary 1
###--------------------------------------------------------------------------###

### define column information
column_info <- tribble(
  ~ id,                    ~group,             ~name,                     ~geom,       ~palette,    ~options,
  "id",                    "method",           "",                        "text",      NA,          list(hjust = 0, width = 1.8),
  "Platform",              "method_info",      "Platform",                "text",      NA,          list(width = 1.5),
  "Prior Information",     "method_info",      "Prior Information",       "text",      NA,          list(width = 3.5),
  "Simulate Groups",       "function",         "Groups",                  "text",      NA,          list(width = 1),
  "Simulate DEGs",         "function",         "DEGs",                    "text",      NA,          list(width = 1),
  "Simulate Batches",      "function",         "Batches",                 "text",      NA,          list(width = 1),
  "Simulate Trajectory",   "function",         "Trajectory",              "text",      NA,          list(width = 1),
  "overall",               "overall",          "Overall",                 "bar",       "palette1",  list(width = 4),
  "accuracy",              "overall",          "Accuracy",                "bar",       "palette2",  list(width = 4),
  "functionality",         "overall",          "Functionality",           "bar",       "palette3",  list(width = 4),
  "scalability",           "overall",          "Scalability",             "bar",       "palette5",  list(width = 4),
  "usability",             "overall",          "Usability",               "bar",       "palette4",  list(width = 4)
)

### column grouping

column_groups <- tribble(
  ~ Experiment,                ~Category,                                        ~group,         ~palette,
  "Method Characteristics",    "",                                               "method",       "palette1",
  "Method Characteristics",    "",                                               "method_info",  "palette1",
  "Method Characteristics",    "The application scenarios",                      "function",     "palette1",
  "Evaluation Summary",        "Scores of overall performace and four criteria", "overall",      "palette1"
)

### row info
row_info <- tribble(
   ~group,         ~ id,                   
   "Class 1",      "SPARSim",
   "Class 1",      "Splat",
   "Class 1",      "powsimR",
   "Class 1",      "SplatPop",
   "Class 1",      "SPsimSeq",
   "Class 2",      "Lun",
   "Class 2",      "scDesign",
   "Class 2",      "muscat",
   "Class 2",      "Lun2",
   "Class 2",      "ESCO",
   "Class 3",      "scDesign2",
   "Class 3",      "POWSC",
   "Class 3",      "hierarchicell",
   "Class 3",      "SparseDC",
   "Class 6",      "ESCO-tree",
   "Class 6",      "ESCO-traj",
   "Class 6",      "TedSim",
   "Class 6",      "SymSim",
   "Class 6",      "VeloSim"
)

### row groups
row_groups <- tribble(
  ~group,         ~ id,                   
  "Class 1",      "Class 1 Method",
  "Class 2",      "Class 2 Method",
  "Class 3",      "Class 3 Method",
  "Class 6",      "Class 6 Method",
)

### determine palettes

palettes <- tribble(
  ~palette,             ~colours,
  "palette1",           grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Greys")[-1]))(101),
  "palette2",           grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues") %>% c("#011636")))(101),
  "palette3",           grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Reds")[-9]))(101),
  "palette4",           grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "YlOrBr")[-9]))(101),
  "palette5",           grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Greens")[-9]))(101)
)

data <- overall_data %>% 
  select(-2) %>% 
  select(1:12)

detailed_plot <- funky_heatmap(data = data,
                               column_info = column_info,
                               column_groups = column_groups,
                               row_info = row_info,
                               row_groups = row_groups,
                               palettes = palettes,
                               scale_column = FALSE,
                               expand = c(xmin = 1, xmax = 2, ymin = 1, ymax = 1))

library(Cairo)
ggsave(detailed_plot,
       filename = "Chunk8-Data Analysis/main_figure/detailed_plot1.pdf",
       units = "in",
       width = 18,
       height = 13, device = cairo_pdf)


###--------------------------------------------------------------------------###
###                            Method Summary 2
###--------------------------------------------------------------------------###

data2 <- overall_data %>% 
  select(-c(2:13))

### define column information
column_info <- column_info <- tribble(
  ~ id,                  ~group,             ~name,                     ~geom,       ~palette,    ~options,
  "id",                  "method",           "",                        "text",      NA,          list(hjust = 0, width = 3),
  "KDE",                 "property",         "KDE",                     "funkyrect", "palette2",  NULL,
  "KS",                  "property",         "KS",                      "funkyrect", "palette2",  NULL,
  "MAD",                 "property",         "MAD",                     "funkyrect", "palette2",  NULL,
  "MAE",                 "property",         "MAE",                     "funkyrect", "palette2",  NULL,
  "OV",                  "property",         "Overlapping",             "funkyrect", "palette2",  NULL,
  "RMSE",                "property",         "RMSE",                    "funkyrect", "palette2",  NULL,
  "bhattacharyya",       "property",         "bhattacharyya",           "funkyrect", "palette2",  NULL,
  "multiKS",             "property",         "multiKS",                 "funkyrect", "palette2",  NULL,
  "CDI",                 "group",            "CDI",                     "funkyrect", "palette3",  NULL,
  "ROUGE",               "group",            "ROUGE",                   "funkyrect", "palette3",  NULL,
  "silhouette",          "group",            "silhouette",              "funkyrect", "palette3",  NULL,
  "dunn",                "group",            "dunn",                    "funkyrect", "palette3",  NULL,
  "connectivity",        "group",            "connectivity",            "funkyrect", "palette3",  NULL,
  "DB_index",            "group",            "DB Index",                "funkyrect", "palette3",  NULL,
  "cms",                 "batch",            "cms",                     "funkyrect", "palette3",  NULL,
  "LISI",                "batch",            "LISI",                    "funkyrect", "palette3",  NULL,
  "mm",                  "batch",            "mm",                      "funkyrect", "palette3",  NULL,
  "shannon_entropy",     "batch",            "Shannon Entropy",         "funkyrect", "palette3",  NULL,
  "kBET",                "batch",            "kBET",                    "funkyrect", "palette3",  NULL,
  "AWS_batch",           "batch",            "AWS_batch",               "funkyrect", "palette3",  NULL,
  "pcr",                 "batch",            "pcr",                     "funkyrect", "palette3",  NULL,
  "distribution_score",  "DEGs",             "Distribution Score",      "funkyrect", "palette3",  NULL,
  "true_proportion",     "DEGs",             "True Proportion",         "funkyrect", "palette3",  NULL,
  "Accuracy",            "DEGs",             "Accuracy",                "funkyrect", "palette3",  NULL,
  "Precision",           "DEGs",             "Precision",               "funkyrect", "palette3",  NULL,
  "Recall",              "DEGs",             "Recall",                  "funkyrect", "palette3",  NULL,
  "F1",                  "DEGs",             "F1 Score",                "funkyrect", "palette3",  NULL,
  "AUC",                 "DEGs",             "AUC",                     "funkyrect", "palette3",  NULL,
  "estimation time",     "scala",            "Estimation Time",         "funkyrect", "palette5",  NULL,
  "estimation memory",   "scala",            "Estimation Memory",       "funkyrect", "palette5",  NULL,
  "simulation time",     "scala",            "Simulation Time",         "funkyrect", "palette5",  NULL,
  "simulation memory",   "scala",            "Simulation Memory",       "funkyrect", "palette5",  NULL,
  "time",                "scala",            "Time",                    "funkyrect", "palette5",  NULL,
  "memory",              "scala",            "Memory",                  "funkyrect", "palette5",  NULL,
  "Avalability",         "usability",        "Avalability",             "funkyrect", "palette4",  NULL,
  "Code",                "usability",        "Code",                    "funkyrect", "palette4",  NULL,
  "Documentation",       "usability",        "Documentation",           "funkyrect", "palette4",  NULL,
  "Evaluation",          "usability",        "Evaluation",              "funkyrect", "palette4",  NULL,
  "Maintenance",         "usability",        "Maintenance",             "funkyrect", "palette4",  NULL,
  "Paper",               "usability",        "Paper",                   "funkyrect", "palette4",  NULL
)

### column grouping

column_groups <- tribble(
  ~ Experiment,     ~Category,                                 ~group,         ~palette,
  "Method",         "",                                        "method",       "palette1",
  "Accuracy",       "Per metric",                              "property",     "palette2",
  "Functionality",  "Group",                                   "group",        "palette3",
  "Functionality",  "Batch",                                   "batch",        "palette3",
  "Functionality",  "DEGs",                                    "DEGs",         "palette3",
  "Scalability",    "Time comsuming and \nmemory usage",       "scala",        "palette5",
  "Usability",      "The quality of the codes and softwares",  "usability",    "palette4"
)


### determine palettes

palettes <- tribble(
  ~palette,             ~colours,
  "palette1",           grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Greys")[-1]))(101),
  "palette2",           grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues") %>% c("#011636")))(101),
  "palette3",           grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Reds")[-9]))(101),
  "palette4",           grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "YlOrBr")[-9]))(101),
  "palette5",           grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Greens")[-9]))(101)
)

detailed_plot <- funky_heatmap(data = data2,
                               column_info = column_info,
                               column_groups = column_groups,
                               row_info = row_info,
                               row_groups = row_groups,
                               palettes = palettes,
                               scale_column = FALSE,
                               expand = c(xmin = 1, xmax = 4, ymin = 1, ymax = 1))

ggsave(detailed_plot, filename = "Chunk8-Data Analysis/main_figure/detailed_plot.pdf", units = "in", width = 28, height = 15)


column_info <- dynbenchmark_data$column_info
