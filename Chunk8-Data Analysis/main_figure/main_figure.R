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

### errors
errors <- openxlsx::read.xlsx(xlsxFile = "./Chunk8-Data Analysis/error_reason/error.xlsx", sheet = 1)
errors <- errors %>% 
  mutate(
    method = str_split(error_file, "_", simplify = TRUE)[, 1],
    data = case_when(
      str_starts(str_split(error_file, "_", simplify = TRUE)[, 2], "data") ~ str_split(error_file, "_", simplify = TRUE)[, 2],
      TRUE ~ str_split(error_file, "_", simplify = TRUE)[, 3]
    ),
    data_type = case_when(
      data %in% paste0("data", 1:101) ~ "single-cell data",
      data %in% paste0("data", 102:152) ~ "spatial data"
    )
  )

### error per method, data_type and error category
method_error <- errors %>% 
  group_by(method, category, data_type) %>% 
  summarise(
    count = n()
  )
saveRDS(method_error, file = "./Chunk8-Data Analysis/error_reason/method_error.rds")

### proportion of errors
simulated_list <- list.files("/Volumes/Elements/sim_bench/simulation_data/")
succeed_methods <- table(str_split(simulated_list, pattern = "_", simplify = TRUE)[, 1]) %>% as.data.frame()
colnames(succeed_methods) <- c("Method", "Succeed")
method_error_count <- method_error %>% 
  group_by(method) %>% 
  summarise(
    Failed = sum(count)
  )
colnames(method_error_count)[1] <- "Method"
methods_error_table <- full_join(succeed_methods, method_error_count, by = "Method") %>% 
  mutate(
    across(c(2,3), ~ replace_na(.x, 0)),
    Total = Succeed + Failed,
    error_proportion = Failed / Total,
    error = paste0(as.character(round(error_proportion, digits = 3) * 100), "%")
  )

### proportion of errors per category
error_proportion_category <- method_error %>% 
  group_by(method) %>% 
  summarise(
    error_proportion_category = count/sum(count),
    category = category
  ) %>% 
  ungroup()
category <- unique(error_proportion_category$category)
error_proportion_category <- map(methods_error_table$Method, .f = function(x){
  values <- error_proportion_category %>% 
    filter(method == x) %>% 
    pull(error_proportion_category)
  if(S4Vectors::isEmpty(values)){
    values <- rep(0, 11)
    names(values) <- category
  }else{
    names(values) <- error_proportion_category %>% 
      filter(method == x) %>% 
      pull(category)
  }
  values
})

error_proportion_category_tibble <- map_dfr(1:length(error_proportion_category), .f = function(x){
  tmp <- error_proportion_category[[x]][category]
  names(tmp) <- category
  tmp[is.na(tmp)] <- 0
  error_proportion_category[[x]] <- tmp
  a <- tibble(
    Method = methods_error_table$Method[x],
    Reason = error_proportion_category[x]
  )
  a
})

error_info <- full_join(methods_error_table, error_proportion_category_tibble, by = "Method")
saveRDS(error_info, file = "./Chunk8-Data Analysis/error_reason/error_info.rds")

### Add to accuracy
accuracy_plot <- accuracy_plot %>% 
  full_join(., error_info, by = "Method")


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
functionality_plot <- functionality_plot %>% 
  mutate(
    Group_score = apply(X = functionality_plot %>% select(3:8), MARGIN = 1, FUN = function(x){mean(x, na.rm = TRUE)}),
    DEGs_score = apply(X = functionality_plot %>% select(9:15), MARGIN = 1, FUN = function(x){mean(x, na.rm = TRUE)}),
    Batch_score = apply(X = functionality_plot %>% select(16:23), MARGIN = 1, FUN = function(x){mean(x, na.rm = TRUE)}),
    Trajectory_score = apply(X = functionality_plot %>% select(24:27), MARGIN = 1, FUN = function(x){mean(x, na.rm = TRUE)})
  )

### scalability score
colnames(scalability) <- c("Method",
                           "estimation time",
                           "estimation memory",
                           "simulation time",
                           "simulation memory",
                           "time",
                           "memory",
                           "scalability")
scalability_plot <- readRDS("./Chunk8-Data Analysis/scalability/scalability_plot_data.rds")
colnames(scalability_plot)[1] <- "Method"
scalability <- scalability %>% 
  full_join(scalability_plot, by = "Method")

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
  full_join(., accuracy_plot, by = "Method") %>% 
  full_join(., functionality_plot, by = "Method") %>% 
  full_join(., scalability, by = "Method") %>% 
  full_join(., usability, by = "Method") %>% 
  relocate("accuracy", .after = `Simulate Trajectory`) %>% 
  relocate("functionality", .after = accuracy) %>% 
  relocate("scalability", .after = functionality) %>% 
  relocate("usability", .after = scalability)

overall_data <- overall_data%>% 
  mutate(
    overall = apply(overall_data %>% select(9:12), MARGIN = 1, function(x){mean(c(x[1], x[2], x[3], x[4]), na.rm = TRUE)})
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
  "id",                    "method",           "",                        "text",      NA,          list(hjust = 0, width = 5),
  "Platform",              "method_info",      "Platform",                "text",      NA,          list(width = 1.5),
  "Prior Information",     "method_info",      "Prior Information",       "text",      NA,          list(width = 8),
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
   "Class 1",      "SCRIP-BP",
   "Class 1",      "Splat",
   "Class 1",      "SCRIP-GP-trendedBCV",
   "Class 1",      "SCRIP-GP-commonBCV",
   "Class 1",      "SCRIP-BGP-trendedBCV",
   "Class 1",      "SCRIP-BGP-commonBCV",
   "Class 1",      "powsimR",
   "Class 1",      "SplatPop",
   "Class 1",      "SPsimSeq",
   "Class 2",      "PROSSTT",
   "Class 2",      "SCRIP-paths",
   "Class 2",      "Splat-paths",
   "Class 2",      "ESCO-tree",
   "Class 2",      "SplatPop-paths",
   "Class 2",      "ESCO-traj",
   "Class 2",      "TedSim",
   "Class 2",      "phenopath",
   "Class 2",      "MFA",
   "Class 2",      "SymSim",
   "Class 2",      "dyngen",
   "Class 2",      "VeloSim",
   "Class 2",      "dyntoy",
   "Class 3",      "Lun",
   "Class 3",      "scDesign",
   "Class 3",      "Lun2",
   "Class 3",      "muscat",
   "Class 3",      "ESCO",
   "Class 3",      "scDD",
   "Class 3",      "scDesign3",
   "Class 3",      "zingeR",
   "Class 3",      "zinbwaveZinger",
   "Class 4",      "scDesign2",
   "Class 4",      "hierarchicell",
   "Class 4",      "scGAN",
   "Class 4",      "POWSC",
   "Class 4",      "SparseDC",
   "Class 4",      "SimBPDD",
   "Class 4",      "BASiCS",
   "Class 5",      "Simple",
   "Class 5",      "Kersplat",
   "Class 5",      "zinbwave",
   "Class 5",      "BEARscc",
   "Class 5",      "dropsim",
   "Class 5",      "CancerInSilico"
)


### row groups
row_groups <- tribble(
  ~group,         ~ id,                   
  "Class 1",      "Class 1 Method",
  "Class 2",      "Class 2 Method",
  "Class 3",      "Class 3 Method",
  "Class 4",      "Class 4 Method",
  "Class 5",      "Class 5 Method"
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
       width = 16,
       height = 24,
       device = cairo_pdf)


###--------------------------------------------------------------------------###
###                            Method Summary 2
###--------------------------------------------------------------------------###

data2 <- overall_data %>% 
  select(-c(2:13))

### define column information
column_info <- column_info <- tribble(
  ~ id,                                     ~group,             ~name,                     ~geom,       ~palette,         ~options,
  "id",                                     "method",           "",                        "text",      NA,               list(hjust = 0, width = 7, size = 4),
  "KDE",                                    "property",         "KDE",                     "funkyrect", "palette2",       NULL,
  "KS",                                     "property",         "KS",                      "funkyrect", "palette2",       NULL,
  "MAD",                                    "property",         "MAD",                     "funkyrect", "palette2",       NULL,
  "MAE",                                    "property",         "MAE",                     "funkyrect", "palette2",       NULL,
  "OV",                                     "property",         "Overlapping",             "funkyrect", "palette2",       NULL,
  "RMSE",                                   "property",         "RMSE",                    "funkyrect", "palette2",       NULL,
  "bhattacharyya",                          "property",         "bhattacharyya",           "funkyrect", "palette2",       NULL,
  "multiKS",                                "property",         "multiKS",                 "funkyrect", "palette2",       NULL,
  "error",                                  "error",            "Proportion of errors",    "text",      NA,               list(hjust = 0, width = 2, size = 4),
  "Reason",                                 "error",            "Reasons",                 "pie",       "error_reasons",  NULL,
  "Group_score",                            "group",            "Group Score",             "funkyrect", "palette3",       NULL,
  "CDI",                                    "group",            "CDI",                     "funkyrect", "palette3",       NULL,
  "ROUGE",                                  "group",            "ROUGE",                   "funkyrect", "palette3",       NULL,
  "silhouette",                             "group",            "silhouette",              "funkyrect", "palette3",       NULL,
  "dunn",                                   "group",            "dunn",                    "funkyrect", "palette3",       NULL,
  "connectivity",                           "group",            "connectivity",            "funkyrect", "palette3",       NULL,
  "DB_index",                               "group",            "DB Index",                "funkyrect", "palette3",       NULL,
  "DEGs_score",                             "DEGs",             "DEGs Score",              "funkyrect", "palette3",       NULL,
  "distribution_score",                     "DEGs",             "Distribution Score",      "funkyrect", "palette3",       NULL,
  "true_proportion",                        "DEGs",             "True Proportion",         "funkyrect", "palette3",       NULL,
  "Accuracy",                               "DEGs",             "Accuracy",                "funkyrect", "palette3",       NULL,
  "Precision",                              "DEGs",             "Precision",               "funkyrect", "palette3",       NULL,
  "Recall",                                 "DEGs",             "Recall",                  "funkyrect", "palette3",       NULL,
  "F1",                                     "DEGs",             "F1 Score",                "funkyrect", "palette3",       NULL,
  "AUC",                                    "DEGs",             "AUC",                     "funkyrect", "palette3",       NULL,
  "Batch_score",                            "batch",            "Batch Score",             "funkyrect", "palette3",       NULL,
  "cms",                                    "batch",            "cms",                     "funkyrect", "palette3",       NULL,
  "LISI",                                   "batch",            "LISI",                    "funkyrect", "palette3",       NULL,
  "mm",                                     "batch",            "mm",                      "funkyrect", "palette3",       NULL,
  "shannon_entropy",                        "batch",            "Shannon Entropy",         "funkyrect", "palette3",       NULL,
  "kBET",                                   "batch",            "kBET",                    "funkyrect", "palette3",       NULL,
  "AWS_batch",                              "batch",            "AWS_batch",               "funkyrect", "palette3",       NULL,
  "pcr",                                    "batch",            "pcr",                     "funkyrect", "palette3",       NULL,
  "ISI",                                    "batch",            "ISI",                     "funkyrect", "palette3",       NULL,
  "Trajectory_score",                       "trajectory",       "Trajectory Score",        "funkyrect", "palette3",       NULL,
  "HIM",                                    "trajectory",       "HIM",                     "funkyrect", "palette3",       NULL,
  "F1_branches",                            "trajectory",       "F1_branches",             "funkyrect", "palette3",       NULL,
  "F1_milestones",                          "trajectory",       "F1_milestones",           "funkyrect", "palette3",       NULL,
  "Cor_dist",                               "trajectory",       "Cor_dist",                "funkyrect", "palette3",       NULL,
  "estimation time",                        "scala",            "Estimation Time",         "funkyrect", "palette5",       NULL,
  "estimation memory",                      "scala",            "Estimation Memory",       "funkyrect", "palette5",       NULL,
  "simulation time",                        "scala",            "Simulation Time",         "funkyrect", "palette5",       NULL,
  "simulation memory",                      "scala",            "Simulation Memory",       "funkyrect", "palette5",       NULL,
  "time",                                   "scala",            "Time",                    "funkyrect", "palette5",       NULL,
  "memory",                                 "scala",            "Memory",                  "funkyrect", "palette5",       NULL,
  "scalability_100_1000_estimation_score",  "scala_rect",       "100x1000",                "rect",      "scala_rect",     NULL,
  "scalability_100_1000_estimation_score",  "scala_rect",       "",                        "text",      "white6black4",   list(label = "scalability_100_1000_estimation_text",
                                                                                                                               overlay = TRUE,
                                                                                                                               size = 2.1),
  "scalability_1000_1000_estimation_score",  "scala_rect",      "1000x1000",               "rect",      "scala_rect",     NULL,
  "scalability_1000_1000_estimation_score",  "scala_rect",      "",                        "text",      "white6black4",   list(label = "scalability_1000_1000_estimation_text",
                                                                                                                               overlay = TRUE,
                                                                                                                               size = 2.1),
  "scalability_1000_10000_estimation_score",  "scala_rect",     "1000x10000",              "rect",      "scala_rect",     NULL,
  "scalability_1000_10000_estimation_score",  "scala_rect",     "",                        "text",      "white6black4",   list(label = "scalability_1000_10000_estimation_text",
                                                                                                                               overlay = TRUE,
                                                                                                                               size = 2.1),
  "scalability_10000_1000_estimation_score",  "scala_rect",     "10000x1000",              "rect",      "scala_rect",     NULL,
  "scalability_10000_1000_estimation_score",  "scala_rect",     "",                        "text",      "white6black4",   list(label = "scalability_10000_1000_estimation_text",
                                                                                                                               overlay = TRUE,
                                                                                                                               size = 2.1),
  "scalability_100_1000_simulation_score",  "scala_rect2",      "100x1000",                "rect",      "scala_rect",     NULL,
  "scalability_100_1000_simulation_score",  "scala_rect2",      "",                        "text",      "white6black4",   list(label = "scalability_100_1000_simulation_text",
                                                                                                                               overlay = TRUE,
                                                                                                                               size = 2.1),
  "scalability_1000_1000_simulation_score",  "scala_rect2",     "1000x1000",               "rect",      "scala_rect",     NULL,
  "scalability_1000_1000_simulation_score",  "scala_rect2",     "",                        "text",      "white6black4",   list(label = "scalability_1000_1000_simulation_text",
                                                                                                                               overlay = TRUE,
                                                                                                                               size = 2.1),
  "scalability_1000_10000_simulation_score",  "scala_rect2",    "1000x10000",              "rect",      "scala_rect",     NULL,
  "scalability_1000_10000_simulation_score",  "scala_rect2",    "",                        "text",      "white6black4",   list(label = "scalability_1000_10000_simulation_text",
                                                                                                                                overlay = TRUE,
                                                                                                                                size = 2.1),
  "scalability_10000_1000_simulation_score",  "scala_rect2",    "10000x1000",              "rect",      "scala_rect",     NULL,
  "scalability_10000_1000_simulation_score",  "scala_rect2",    "",                        "text",      "white6black4",   list(label = "scalability_10000_1000_simulation_text",
                                                                                                                                overlay = TRUE,
                                                                                                                                size = 2.1),
  "Avalability",                            "usability",        "Avalability",             "funkyrect", "palette4",       NULL,
  "Code",                                   "usability",        "Code",                    "funkyrect", "palette4",       NULL,
  "Documentation",                          "usability",        "Documentation",           "funkyrect", "palette4",       NULL,
  "Evaluation",                             "usability",        "Evaluation",              "funkyrect", "palette4",       NULL,
  "Maintenance",                            "usability",        "Maintenance",             "funkyrect", "palette4",       NULL,
  "Paper",                                  "usability",        "Paper Quality",           "funkyrect", "palette4",       NULL
)

### column grouping

column_groups <- tribble(
  ~ Experiment,     ~Category,                                       ~group,         ~palette,
  "Method",         "",                                              "method",       "palette1",
  "Accuracy",       "Per metric",                                    "property",     "palette2",
  "Accuracy",       "Error",                                         "error",        "palette2",
  "Functionality",  "Group",                                         "group",        "palette3",
  "Functionality",  "DEGs",                                          "DEGs",         "palette3",
  "Functionality",  "Batch",                                         "batch",        "palette3",
  "Functionality",  "trajectory",                                    "trajectory",   "palette3",
  "Scalability",    "Scores of time consuming \nand memory usage",   "scala",        "palette5",
  "Scalability",    "Estimation \n(cell x gene)",                    "scala_rect",   "palette5",
  "Scalability",    "Simulation \n(cell x gene)",                    "scala_rect2",  "palette5",
  "Usability",      "The quality of the codes \nand softwares",      "usability",    "palette4"
)


### determine palettes
error_reasons_palettes <- RColorBrewer::brewer.pal(11, "Paired")
names(error_reasons_palettes) <- category

palettes <- tribble(
  ~palette,             ~colours,
  "palette1",           grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Greys")[-1]))(101),
  "palette2",           grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues") %>% c("#011636")))(101),
  "palette3",           grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Reds")[-9]))(101),
  "palette4",           grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "YlOrBr")[-9]))(101),
  "palette5",           grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Greens")[-9]))(101),
  "error_reasons",      error_reasons_palettes,
  "scala_rect",         grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Reds")[-9]))(101),
  "white6black4",       c("white","white","white","black","black","black","black","black","black","black")
)

detailed_plot <- funky_heatmap(data = data2,
                               column_info = column_info,
                               column_groups = column_groups,
                               row_info = row_info,
                               row_groups = row_groups,
                               palettes = palettes,
                               scale_column = FALSE,
                               expand = c(xmin = 1, xmax = 4, ymin = 1, ymax = 1))

ggsave(detailed_plot,
       filename = "Chunk8-Data Analysis/main_figure/detailed_plot.pdf",
       units = "in",
       width = 24,
       height = 28,
       device = cairo_pdf)


