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
source("Chunk8-Data Analysis/main_figure/utils.R")
accuracy <- readRDS("Chunk8-Data Analysis/accuracy/accuracy_long_data.rds")
functionality <- readRDS("Chunk8-Data Analysis/functionality/functionality_long_data.rds")
scalability <- readRDS("Chunk8-Data Analysis/scalability/score_data.rds")
usability <- readRDS("Chunk8-Data Analysis/usability/usability_long_data.rds")
method <- openxlsx::read.xlsx("Chunk1-Data preparation/methods.xlsx", sheet = 1, colNames = TRUE) %>% 
  select(-2)
colnames(method) <- c("Method",
                      "Category",
                      "Platform",
                      "Model",
                      "Prior Information",
                      "Simulate Groups",
                      "Simulate DEGs",
                      "Simulate Batches",
                      "Simulate Trajectory")
method_names <- method$Method
data_info <- openxlsx::read.xlsx("./Chunk1-Data preparation/evaluation_datasets.xlsx", rowNames = TRUE) %>% 
  mutate(
    Data = str_split(Dataset.ID, "_", simplify = TRUE)[, 1]
  ) %>% 
  select(Data, Platform)

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

accuracy_plot <- accuracy_process_function(accuracy)


############################# errors
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
succeed_methods <- succeed_methods %>% 
  rbind(data.frame("Var1" = c("scMultiSim", "scMultiSim-tree", "scDesign3-traj", "SRTsim"),
                   "Freq" = c(10, 10, 10, 10)))
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
  reframe(
    error_proportion_category = count/sum(count),
    category = category
  )
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
functionality <- functionality %>% 
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

functionality_plot <- functionality_process_function(functionality)

### scalability score
colnames(scalability) <- c("Method",
                           "estimation time",
                           "estimation memory",
                           "simulation time",
                           "simulation memory",
                           "time",
                           "memory",
                           "scalability",
                           "Cor(time_estimation)",
                           "Cor(memory_estimation)",
                           "Cor(time_simulation)",
                           "Cor(memory_simulation)")
scalability_plot <- readRDS("./Chunk8-Data Analysis/scalability/scalability_plot_data.rds")
colnames(scalability_plot)[1] <- "Method"
scalability <- scalability %>% 
  full_join(scalability_plot, by = "Method")
scalability$`Cor(time_estimation)` <- as.character(scalability$`Cor(time_estimation)`)
scalability$`Cor(memory_estimation)` <- as.character(scalability$`Cor(memory_estimation)`)
scalability$`Cor(time_simulation)` <- as.character(scalability$`Cor(time_simulation)`)
scalability$`Cor(memory_simulation)` <- as.character(scalability$`Cor(memory_simulation)`)

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
  left_join(., accuracy_plot, by = "Method") %>% 
  left_join(., functionality_plot, by = "Method") %>% 
  left_join(., scalability, by = "Method") %>% 
  left_join(., usability, by = "Method") %>% 
  relocate("accuracy", .after = `Simulate Trajectory`) %>% 
  relocate("functionality", .after = accuracy) %>% 
  relocate("scalability", .after = functionality) %>% 
  relocate("usability", .after = scalability)

overall_data <- overall_data %>% 
  mutate(
    overall = apply(overall_data %>% select(10:13), MARGIN = 1, function(x){mean(c(x[1], x[2], x[3], x[4]), na.rm = TRUE)})
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
    across(all_of(c("Model",
                    "Prior Information",
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
saveRDS(overall_data, file = "./Chunk8-Data Analysis/overall_data.rds")
###--------------------------------------------------------------------------###
###                            Fig2. Summary
###--------------------------------------------------------------------------###

### define column information
column_info <- tribble(
  ~ id,                    ~group,             ~name,                     ~geom,       ~palette,    ~options,
  "id",                    "method",           "",                        "text",      NA,          list(hjust = 0, width = 5),
  "Platform",              "method_info",      "Platform",                "text",      NA,          list(width = 1.5),
  "Model",                 "method_info",      "Model",                   "text",      NA,          list(width = 8),
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
   "Class 1",      "Splat",
   "Class 1",      "SCRIP-BP",
   "Class 1",      "SCRIP-GP-commonBCV",
   "Class 1",      "SCRIP-BGP-commonBCV",
   "Class 1",      "SCRIP-GP-trendedBCV",
   "Class 1",      "SCRIP-BGP-trendedBCV",
   "Class 1",      "powsimR",
   "Class 1",      "SplatPop",
   "Class 1",      "SPsimSeq",
   "Class 2",      "PROSSTT",
   "Class 2",      "SCRIP-paths",
   "Class 2",      "Splat-paths",
   "Class 2",      "ESCO-tree",
   "Class 2",      "SplatPop-paths",
   "Class 2",      "ESCO-traj",
   "Class 2",      "phenopath",
   "Class 2",      "TedSim",
   "Class 2",      "MFA",
   "Class 2",      "SymSim",
   "Class 2",      "dyngen",
   "Class 2",      "VeloSim",
   "Class 2",      "dyntoy",
   "Class 2",      "scMultiSim-tree",
   "Class 2",      "scDesign3-traj",
   "Class 3",      "scDesign",
   "Class 3",      "Lun",
   "Class 3",      "muscat",
   "Class 3",      "Lun2",
   "Class 3",      "ESCO",
   "Class 3",      "scDD",
   "Class 3",      "scDesign3",
   "Class 3",      "zingeR",
   "Class 3",      "zinbwaveZinger",
   "Class 4",      "scMultiSim",
   "Class 4",      "scDesign2",
   "Class 4",      "hierarchicell",
   "Class 4",      "POWSC",
   "Class 4",      "SparseDC",
   "Class 4",      "scGAN",
   "Class 4",      "BASiCS",
   "Class 4",      "SimBPDD",
   "Class 4",      "SRTsim",
   "Class 5",      "Simple",
   "Class 5",      "Kersplat",
   "Class 5",      "BEARscc",
   "Class 5",      "zinbwave",
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

summary_plot <- funky_heatmap(data = overall_data,
                              column_info = column_info,
                              column_groups = column_groups,
                              row_info = row_info,
                              row_groups = row_groups,
                              palettes = palettes,
                              scale_column = FALSE,
                              expand = c(xmin = 1, xmax = 2, ymin = 1, ymax = 1))

library(Cairo)
ggsave(summary_plot,
       filename = "/Users/duohongrui/Desktop/sim-article/figures/Fig2.pdf",
       units = "in",
       width = 17,
       height = 24,
       device = cairo_pdf)

###--------------------------------------------------------------------------###
###                              Supplementary Fig.1
###--------------------------------------------------------------------------###

### define column information
column_info <- column_info <- tribble(
  ~ id,                                     ~group,             ~name,                     ~geom,       ~palette,         ~options,
  "id",                                     "method",           "",                        "text",      NA,               list(hjust = 0,
                                                                                                                               width = 4.5,
                                                                                                                               size = 4,
                                                                                                                               legend = FALSE),
  "KDE",                                    "property",         "KDE",                     "funkyrect", "palette2",       NULL,
  "KS",                                     "property",         "KS",                      "funkyrect", "palette2",       NULL,
  "MAD",                                    "property",         "MAD",                     "funkyrect", "palette2",       NULL,
  "MAE",                                    "property",         "MAE",                     "funkyrect", "palette2",       NULL,
  "OV",                                     "property",         "Overlapping",             "funkyrect", "palette2",       NULL,
  "RMSE",                                   "property",         "RMSE",                    "funkyrect", "palette2",       NULL,
  "bhattacharyya",                          "property",         "bhattacharyya",           "funkyrect", "palette2",       NULL,
  "multiKS",                                "property",         "multiKS",                 "funkyrect", "palette2",       NULL,
  "acc_scRNA-seq data",                     "technique",        "scRNA-seq",               "funkyrect", "palette2",       NULL,
  "acc_spatial transcriptome data",         "technique",        "ST technology",           "funkyrect", "palette2",       NULL,
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
  "Trajectory_score",                       "trajectory",       "Trajectory Score",        "funkyrect", "palette3",       NULL,
  "HIM",                                    "trajectory",       "HIM",                     "funkyrect", "palette3",       NULL,
  "F1_branches",                            "trajectory",       "F1_branches",             "funkyrect", "palette3",       NULL,
  "F1_milestones",                          "trajectory",       "F1_milestones",           "funkyrect", "palette3",       NULL,
  "Cor_dist",                               "trajectory",       "Cor_dist",                "funkyrect", "palette3",       NULL,
  "func_scRNA-seq data",                    "func_technique",   "scRNA-seq",               "funkyrect", "palette3",       NULL,
  "func_spatial transcriptome data",        "func_technique",   "ST technology",           "funkyrect", "palette3",       NULL,
  "estimation time",                        "scala",            "Estimation Time",         "funkyrect", "palette5",       NULL,
  "Cor(time_estimation)",                   "scala",            "Cor(estimation time)",    "text",      NA,               list(hjust = 0, width = 1, size = 3),
  "estimation memory",                      "scala",            "Estimation Memory",       "funkyrect", "palette5",       NULL,
  "Cor(memory_estimation)",                 "scala",            "Cor(estimation memory)",  "text",      NA,               list(hjust = 0, width = 1, size = 3),
  "simulation time",                        "scala",            "Simulation Time",         "funkyrect", "palette5",       NULL,
  "Cor(time_simulation)",                   "scala",            "Cor(simulation time)",    "text",      NA,               list(hjust = 0, width = 1, size = 3),
  "simulation memory",                      "scala",            "Simulation Memory",       "funkyrect", "palette5",       NULL,
  "Cor(memory_simulation)",                 "scala",            "Cor(simulation memory)",  "text",      NA,               list(hjust = 0, width = 1, size = 3),
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
  "Paper",                                  "usability",        "Paper Quality",           "funkyrect", "palette4",       NULL,
  "error",                                  "error",            "Proportion of errors",    "text",      NA,               list(hjust = 0, width = 2, size = 4),
  "Reason",                                 "error",            "Reasons",                 "pie",       "error_reasons",  NULL
)

### column grouping

column_groups <- tribble(
  ~ Experiment,     ~Category,                                       ~group,           ~palette,
  "Method",         "",                                              "method",         "palette1",
  "Accuracy",       "Per metric",                                    "property",       "palette2",
  "Accuracy",       "Per technique",                                 "technique",      "palette2",
  "Functionality",  "Group",                                         "group",          "palette3",
  "Functionality",  "DEGs",                                          "DEGs",           "palette3",
  "Functionality",  "Batch",                                         "batch",          "palette3",
  "Functionality",  "trajectory",                                    "trajectory",     "palette3",
  "Functionality",  "Per technique",                                 "func_technique", "palette3",
  "Scalability",    "Scores of time consuming \nand memory usage",   "scala",          "palette5",
  "Scalability",    "Estimation \n(cell x gene)",                    "scala_rect",     "palette5",
  "Scalability",    "Simulation \n(cell x gene)",                    "scala_rect2",    "palette5",
  "Usability",      "The quality of the codes \nand softwares",      "usability",      "palette4",
  "Usability",      "Error",                                         "error",          "palette4",
)


### determine palettes
error_reasons_palettes <- RColorBrewer::brewer.pal(11, "Paired")
a <- overall_data[1, ] %>% pull("Reason")
category <- names(a[[1]])
rm(a)
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

detailed_plot <- funky_heatmap(data = overall_data,
                               column_info = column_info,
                               column_groups = column_groups,
                               row_info = row_info,
                               row_groups = row_groups,
                               palettes = palettes,
                               scale_column = FALSE,
                               expand = c(xmin = 1, xmax = 3.5, ymin = 0.5, ymax = 0.5),
                               col_annot_offset = 4.5,
                               removed_entries = c("SPsimseq"))

ggsave(detailed_plot,
       filename = "/Users/duohongrui/Desktop/sim-article/figures/Supp_Fig1.pdf",
       units = "in",
       width = 20,
       height = 16,
       device = cairo_pdf)


###--------------------------------------------------------------------------###
###                             Fig.3
###--------------------------------------------------------------------------###

### define column information
column_info <- column_info <- tribble(
  ~ id,                                     ~group,             ~name,                     ~geom,       ~palette,         ~options,
  "id",                                     "method",           "",                        "text",      NA,               list(hjust = 0,
                                                                                                                               width = 4,
                                                                                                                               size = 4, 
                                                                                                                               egend = FALSE),
  "KDE",                                    "property",         "KDE",                     "funkyrect", "palette2",       NULL,
  "KS",                                     "property",         "KS",                      "funkyrect", "palette2",       NULL,
  "MAD",                                    "property",         "MAD",                     "funkyrect", "palette2",       NULL,
  "MAE",                                    "property",         "MAE",                     "funkyrect", "palette2",       NULL,
  "OV",                                     "property",         "Overlapping",             "funkyrect", "palette2",       NULL,
  "RMSE",                                   "property",         "RMSE",                    "funkyrect", "palette2",       NULL,
  "bhattacharyya",                          "property",         "bhattacharyya",           "funkyrect", "palette2",       NULL,
  "multiKS",                                "property",         "multiKS",                 "funkyrect", "palette2",       NULL,
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
  "Trajectory_score",                       "trajectory",       "Trajectory Score",        "funkyrect", "palette3",       NULL,
  "HIM",                                    "trajectory",       "HIM",                     "funkyrect", "palette3",       NULL,
  "F1_branches",                            "trajectory",       "F1_branches",             "funkyrect", "palette3",       NULL,
  "F1_milestones",                          "trajectory",       "F1_milestones",           "funkyrect", "palette3",       NULL,
  "Cor_dist",                               "trajectory",       "Cor_dist",                "funkyrect", "palette3",       NULL,
  "estimation time",                        "scala",            "Estimation Time",         "funkyrect", "palette5",       NULL,
  "estimation memory",                      "scala",            "Estimation Memory",       "funkyrect", "palette5",       NULL,
  "simulation time",                        "scala",            "Simulation Time",         "funkyrect", "palette5",       NULL,
  "simulation memory",                      "scala",            "Simulation Memory",       "funkyrect", "palette5",       NULL,
  "Avalability",                            "usability",        "Avalability",             "funkyrect", "palette4",       NULL,
  "Code",                                   "usability",        "Code",                    "funkyrect", "palette4",       NULL,
  "Documentation",                          "usability",        "Documentation",           "funkyrect", "palette4",       NULL,
  "Evaluation",                             "usability",        "Self-evaluation",         "funkyrect", "palette4",       NULL,
  "Maintenance",                            "usability",        "Maintenance",             "funkyrect", "palette4",       NULL,
  "Paper",                                  "usability",        "Paper Quality",           "funkyrect", "palette4",       NULL,
  "error",                                  "error",            "Proportion of errors",    "text",      NA,               list(hjust = 0, width = 2, size = 4),
  "Reason",                                 "error",            "Reasons",                 "pie",       "error_reasons",  NULL
)

### column grouping

column_groups <- tribble(
  ~ Experiment,     ~Category,                                       ~group,           ~palette,
  "Method",         "",                                              "method",         "palette1",
  "Accuracy",       "Per metric",                                    "property",       "palette2",
  "Functionality",  "Group",                                         "group",          "palette3",
  "Functionality",  "DEGs",                                          "DEGs",           "palette3",
  "Functionality",  "Batch",                                         "batch",          "palette3",
  "Functionality",  "Trajectory",                                    "trajectory",     "palette3",
  "Scalability",    "Time consuming \n memory usage",             "scala",          "palette5",
  "Usability",      "Quality of the codes \nand softwares",          "usability",      "palette4",
  "Usability",      "Error",                                         "error",          "palette4",
)


fig2_plot <- funky_heatmap(data = overall_data,
                           column_info = column_info,
                           column_groups = column_groups,
                           row_info = row_info,
                           row_groups = row_groups,
                           palettes = palettes,
                           scale_column = FALSE,
                           expand = c(xmin = 1, xmax = 4, ymin = 0.5, ymax = 0.5),
                           col_annot_offset = 4)

ggsave(fig2_plot,
       filename = "/Users/duohongrui/Desktop/sim-article/figures/Fig3.pdf",
       units = "in",
       width = 16,
       height = 19)


###--------------------------------------------------------------------------###
###                       Supplementary Fig.2 (Accuracy)
###--------------------------------------------------------------------------###

### define column information
column_info <- column_info <- tribble(
  ~ id,                                     ~group,             ~name,                     ~geom,       ~palette,         ~options,
  "id",                                     "method",           "",                        "text",      NA,               list(hjust = 0, width = 7, size = 4, legend = FALSE),
  "KDE",                                    "metric",           "KDE",                     "funkyrect", "palette2",       NULL,
  "KS",                                     "metric",           "KS",                      "funkyrect", "palette2",       NULL,
  "MAD",                                    "metric",           "MAD",                     "funkyrect", "palette2",       NULL,
  "MAE",                                    "metric",           "MAE",                     "funkyrect", "palette2",       NULL,
  "OV",                                     "metric",           "Overlapping",             "funkyrect", "palette2",       NULL,
  "RMSE",                                   "metric",           "RMSE",                    "funkyrect", "palette2",       NULL,
  "bhattacharyya",                          "metric",           "bhattacharyya",           "funkyrect", "palette2",       NULL,
  "multiKS",                                "metric",           "multiKS",                 "funkyrect", "palette2",       NULL,
  "LS",                                     "property",         "LS",                      "funkyrect", "palette2",       NULL,
  "FZC",                                    "property",         "FZC",                     "funkyrect", "palette2",       NULL,
  "CCC",                                    "property",         "CCC",                     "funkyrect", "palette2",       NULL,
  "TMM",                                    "property",         "TMM",                     "funkyrect", "palette2",       NULL,
  "ELS",                                    "property",         "ELS",                     "funkyrect", "palette2",       NULL,
  "FCO",                                    "property",         "FCO",                     "funkyrect", "palette2",       NULL,
  "RLZ",                                    "property",         "RLZ",                     "funkyrect", "palette2",       NULL,
  "ME",                                     "property",         "ME",                      "funkyrect", "palette2",       NULL,
  "SD",                                     "property",         "SD",                      "funkyrect", "palette2",       NULL,
  "CV",                                     "property",         "CV",                      "funkyrect", "palette2",       NULL,
  "FZG",                                    "property",         "FZG",                     "funkyrect", "palette2",       NULL,
  "FGO",                                    "property",         "FGO",                     "funkyrect", "palette2",       NULL,
  "RMS",                                    "property",         "RMS",                     "funkyrect", "palette2",       NULL,
  "RMZ",                                    "property",         "RMZ",                     "funkyrect", "palette2",       NULL,
  "RDM",                                    "property",         "RDM",                     "funkyrect", "palette2",       NULL,
  "acc_MARS-Seq",                           "platform",         "MARS-Seq",                "funkyrect", "palette2",       NULL,
  "acc_10X Genomics",                       "platform",         "10X Genomics",            "funkyrect", "palette2",       NULL,
  "acc_Smart-seq2",                         "platform",         "Smart-seq2",              "funkyrect", "palette2",       NULL,
  "acc_Mix sources1",                       "platform",         "Mix sources1",            "funkyrect", "palette2",       NULL,
  "acc_CEL-seq",                            "platform",         "CEL-seq",                 "funkyrect", "palette2",       NULL,
  "acc_Fluidigm C1",                        "platform",         "Fluidigm C1",             "funkyrect", "palette2",       NULL,
  "acc_CEL-seq2",                           "platform",         "CEL-seq2",                "funkyrect", "palette2",       NULL,
  "acc_Mix sources2",                       "platform",         "Mix sources2",            "funkyrect", "palette2",       NULL,
  "acc_Drop-seq",                           "platform",         "Drop-seq",                "funkyrect", "palette2",       NULL,
  "acc_inDrop",                             "platform",         "inDrop",                  "funkyrect", "palette2",       NULL,
  "acc_Microwell-seq",                      "platform",         "Microwell-seq",           "funkyrect", "palette2",       NULL,
  "acc_Smart-seq",                          "platform",         "Smart-seq",               "funkyrect", "palette2",       NULL,
  "acc_ST",                                 "platform",         "ST",                      "funkyrect", "palette2",       NULL,
  "acc_HDST",                               "platform",         "HDST",                    "funkyrect", "palette2",       NULL,
  "acc_10X Visium",                         "platform",         "10X Visium",              "funkyrect", "palette2",       NULL,
  "acc_Slide-Seq",                          "platform",         "Slide-Seq",               "funkyrect", "palette2",       NULL,
  "acc_Slide-SeqV2",                        "platform",         "Slide-SeqV2",             "funkyrect", "palette2",       NULL,
  "acc_seqFISH",                            "platform",         "seqFISH",                 "funkyrect", "palette2",       NULL,
  "acc_seqFISH+",                           "platform",         "seqFISH+",                "funkyrect", "palette2",       NULL,
  "acc_osmFISH",                            "platform",         "osmFISH",                 "funkyrect", "palette2",       NULL,
  "acc_sci-Space",                          "platform",         "sci-Space",               "funkyrect", "palette2",       NULL,
  "acc_MERFISH",                            "platform",         "MERFISH",                 "funkyrect", "palette2",       NULL,
  "acc_Stereo-Seq",                         "platform",         "Stereo-Seq",              "funkyrect", "palette2",       NULL,
  "acc_scRNA-seq data",                     "technique",        "scRNA-seq",               "funkyrect", "palette2",       NULL,
  "acc_spatial transcriptome data",         "technique",        "ST technology",           "funkyrect", "palette2",       NULL
)

### column grouping

column_groups <- tribble(
  ~ Experiment,     ~Category,                                       ~group,           ~palette,
  "Method",         "",                                              "method",         "palette1",
  "Accuracy",       "Per metric",                                    "metric",         "palette2",
  "Accuracy",       "Per property",                                  "property",       "palette2",
  "Accuracy",       "Per platform",                                  "platform",       "palette2",
  "Accuracy",       "Per technique",                                 "technique",      "palette2"
)


accuracy_sup_plot <- funky_heatmap(data = overall_data,
                                   column_info = column_info,
                                   column_groups = column_groups,
                                   row_info = row_info,
                                   row_groups = row_groups,
                                   palettes = palettes,
                                   scale_column = FALSE,
                                   expand = c(xmin = 1, xmax = 4, ymin = 0.5, ymax = 0.5),
                                   col_annot_offset = 4)

ggsave(accuracy_sup_plot,
       filename = "/Users/duohongrui/Desktop/sim-article/figures/Supp_Fig2.pdf",
       units = "in",
       width = 16,
       height = 19)


###--------------------------------------------------------------------------###
###                       Supplementary Fig.3 (Functionality)
###--------------------------------------------------------------------------###

### define column information
column_info <- column_info <- tribble(
  ~ id,                                     ~group,             ~name,                     ~geom,       ~palette,         ~options,
  "id",                                     "method",           "",                        "text",      NA,               list(hjust = 0, width = 7, size = 4, legend = FALSE),
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
  "Trajectory_score",                       "trajectory",       "Trajectory Score",        "funkyrect", "palette3",       NULL,
  "HIM",                                    "trajectory",       "HIM",                     "funkyrect", "palette3",       NULL,
  "F1_branches",                            "trajectory",       "F1_branches",             "funkyrect", "palette3",       NULL,
  "F1_milestones",                          "trajectory",       "F1_milestones",           "funkyrect", "palette3",       NULL,
  "Cor_dist",                               "trajectory",       "Cor_dist",                "funkyrect", "palette3",       NULL,
  "func_MARS-Seq",                          "platform",         "MARS-Seq",                "funkyrect", "palette3",       NULL,
  "func_10X Genomics",                      "platform",         "10X Genomics",            "funkyrect", "palette3",       NULL,
  "func_Smart-seq2",                        "platform",         "Smart-seq2",              "funkyrect", "palette3",       NULL,
  "func_Mix sources1",                      "platform",         "Mix sources1",            "funkyrect", "palette3",       NULL,
  "func_CEL-seq",                           "platform",         "CEL-seq",                 "funkyrect", "palette3",       NULL,
  "func_Fluidigm C1",                       "platform",         "Fluidigm C1",             "funkyrect", "palette3",       NULL,
  "func_CEL-seq2",                          "platform",         "CEL-seq2",                "funkyrect", "palette3",       NULL,
  "func_Mix sources2",                      "platform",         "Mix sources2",            "funkyrect", "palette3",       NULL,
  "func_Drop-seq",                          "platform",         "Drop-seq",                "funkyrect", "palette3",       NULL,
  "func_inDrop",                            "platform",         "inDrop",                  "funkyrect", "palette3",       NULL,
#  "func_Microwell-seq",                     "platform",         "Microwell-seq",           "funkyrect", "palette3",       NULL,
  "func_Smart-seq",                         "platform",         "Smart-seq",               "funkyrect", "palette3",       NULL,
  "func_ST",                                "platform",         "ST",                      "funkyrect", "palette3",       NULL,
  "func_HDST",                              "platform",         "HDST",                    "funkyrect", "palette3",       NULL,
  "func_10X Visium",                        "platform",         "10X Visium",              "funkyrect", "palette3",       NULL,
  "func_Slide-Seq",                         "platform",         "Slide-Seq",               "funkyrect", "palette3",       NULL,
  "func_Slide-SeqV2",                       "platform",         "Slide-SeqV2",             "funkyrect", "palette3",       NULL,
  "func_seqFISH",                           "platform",         "seqFISH",                 "funkyrect", "palette3",       NULL,
  "func_seqFISH+",                          "platform",         "seqFISH+",                "funkyrect", "palette3",       NULL,
  "func_osmFISH",                           "platform",         "osmFISH",                 "funkyrect", "palette3",       NULL,
  "func_sci-Space",                         "platform",         "sci-Space",               "funkyrect", "palette3",       NULL,
  "func_MERFISH",                           "platform",         "MERFISH",                 "funkyrect", "palette3",       NULL,
  "func_Stereo-Seq",                        "platform",         "Stereo-Seq",              "funkyrect", "palette3",       NULL,
  "func_scRNA-seq data",                    "technique",        "scRNA-seq",               "funkyrect", "palette3",       NULL,
  "func_spatial transcriptome data",        "technique",        "ST technology",           "funkyrect", "palette3",       NULL
)

### column grouping

column_groups <- tribble(
  ~ Experiment,     ~Category,            ~group,           ~palette,
  "Method",         "",                   "method",         "palette1",
  "Functionality",  "Group",              "group",          "palette3",
  "Functionality",  "DEGs",               "DEGs",           "palette3",
  "Functionality",  "Batch",              "batch",          "palette3",
  "Functionality",  "Trajectory",         "trajectory",     "palette3",
  "Functionality",  "Per platform",       "platform",       "palette3",
  "Functionality",  "Per technique",      "technique",      "palette3"
)


functionality_sup_plot <- funky_heatmap(data = overall_data,
                                        column_info = column_info,
                                        column_groups = column_groups,
                                        row_info = row_info,
                                        row_groups = row_groups,
                                        palettes = palettes,
                                        scale_column = FALSE,
                                        expand = c(xmin = 1, xmax = 4, ymin = 0.5, ymax = 0.5),
                                        col_annot_offset = 4)

ggsave(functionality_sup_plot,
       filename = "/Users/duohongrui/Desktop/sim-article/figures/Supp_Fig3.pdf",
       units = "in",
       width = 18,
       height = 19)
