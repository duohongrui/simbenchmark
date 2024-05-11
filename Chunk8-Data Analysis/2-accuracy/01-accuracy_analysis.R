library(stringr)
library(purrr)
library(dplyr)
library(tidyr)

### names of metrics and properties
a <- readRDS("../data_properties_result/BEARscc_ERCC_data1_GSE54006.rds")
property_name <- data.frame("abbr" = c("LS", "FZC", "CCC", "TMM", "ELS", "FCO", "RLZ", "ME", "SD", "CV", "FZG", "FGO", "DIS","RMS", "RMZ", "RDM"),
                            "name" = c("library",
                                       "zero fraction of cells",
                                       "cell correlation",
                                       "TMM",
                                       "effective library", "cell outlier",
                                       "library size vs zero fraction of cells",
                                       "mean expression",
                                       "sd",
                                       "cv",
                                       "zero fraction of genes",
                                       "gene outlier",
                                       "dispersion",
                                       "mean expression vs sd",
                                       "mean expression vs zero fraction of genes",
                                       "mean expression vs dispersion"))
metrics <- str_split(names(a), pattern = "_", simplify = T)[,1]
metrics[c(29:33, 64:68, 34:35, 69:74)] <- c(rep("bhattacharyya", 10), "multiKS", "KDE", rep("multiKS", 3), rep("KDE", 3))
property <- c("LS", "FZC", "CCC", "TMM", "ELS", "FCO",
              "LS", "FZC", "CCC", "TMM", "ELS",
              "LS", "FZC", "CCC", "TMM", "ELS", "FCO",
              "LS", "FZC", "CCC", "TMM", "ELS", "FCO",
              "LS", "FZC", "CCC", "TMM", "ELS",
              "LS", "FZC", "CCC", "TMM", "ELS",
              "RLZ", "RLZ",
              "ME", "SD", "CV", "FZG", "DIS", "FGO",
              "ME", "SD", "CV", "FZG", "DIS",
              "ME", "SD", "CV", "FZG", "DIS", "FGO",
              "ME", "SD", "CV", "FZG", "DIS", "FGO",
              "ME", "SD", "CV", "FZG", "DIS",
              "ME", "SD", "CV", "FZG", "DIS",
              "RMS", "RMZ", "RDM",
              "RMS", "RMZ", "RDM")
metrics_property <- paste0(metrics, "_", property)
rm(a)


### Data to Tibble
data_list <- list.files("../data_properties_result/")
all_result <- list()
for (i in data_list) {
  result <- readRDS(file.path("../data_properties_result", i))
  all_result[[i]] <- result
}

accuracy_data <- purrr::map_dfr(1:length(all_result), .f = function(index){
  result <- as.numeric(all_result[[index]])
  names(result) <- metrics_property
  data_name <- names(all_result)[index]
  split_names <- str_split(data_name, pattern = "_", simplify = TRUE)
  method_name <- split_names[, 1]
  names(method_name) <- "Method"
  if(str_starts(split_names[, 2], pattern = "data")){
    data_name <- split_names[, 2]
  }else{
    data_name <- split_names[, 3]
  }
  names(data_name) <- "Data"
  c(method_name, data_name, result)
}) %>% 
  mutate(
    across(all_of(metrics_property), as.numeric)
  )

### Delete DIS
DIS_col <- grep("DIS", colnames(accuracy_data))
accuracy_data <- accuracy_data[, -DIS_col]
metrics_property <- metrics_property[-(DIS_col-2)]
### normalize values which are not in [0, 1] for every dataset
check_range <- purrr::map_dfr(.x = metrics_property, .f = function(x){
  range_check <- as.numeric(range(accuracy_data %>% pull(x), na.rm = TRUE))
  names(range_check) <- c("low", "high")
  range_check
})
rownames(check_range) <- metrics_property
normalized_columes <- colnames(accuracy_data)[grep(colnames(accuracy_data), pattern = "MAD|MAE|RMSE|bhattacharyya|KDE")]

accuracy_data <- accuracy_data %>% 
  group_by(Data) %>% 
  mutate(
    across(all_of(normalized_columes), ~ pnorm((.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
  ) %>% 
  ungroup()

### Subtract values by 1
subtract_columes_prefix <- c("MAD", "KS", "MAE", "RMSE", "bhattacharyya", "multiKS", "KDE")
subtract_columes <- colnames(accuracy_data)[which(str_split(colnames(accuracy_data),
                                                            pattern = "_", simplify = T)[,1] %in% subtract_columes_prefix)]
accuracy_data <- accuracy_data %>% 
  mutate(
    across(all_of(subtract_columes), ~ 1 - .x)
  )

### scale for every metric [0, 1]
accuracy_data <- accuracy_data %>% 
  mutate(
    across(all_of(colnames(accuracy_data)[3:70]), ~ (.x - min(.x, na.rm = TRUE)) / (max(.x, na.rm = TRUE) - min(.x, na.rm = TRUE)))
  )
saveRDS(accuracy_data, file = "Chunk8-Data Analysis/2-accuracy/accuracy_data.rds")
### turn to long table
accuracy_long_data <- accuracy_data %>% 
  pivot_longer(., cols = 3:ncol(.), names_to = "metric_property", values_to = "value")

### split property to a single variable
accuracy_long_data <- accuracy_long_data %>% 
  separate(., col = metric_property, sep = "_", into = c("metric", "property"))
saveRDS(accuracy_long_data, file = "Chunk8-Data Analysis/2-accuracy/accuracy_long_data.rds")



###--------------------------------------------------------------------------###
###                      Accuracy scores for spatial data
###--------------------------------------------------------------------------###
data_list <- list.files("../spatial_metrics_results/")
all_result <- list()
for (i in data_list) {
  result <- readRDS(file.path("../spatial_metrics_results", i))
  all_result[[i]] <- result
}
spatial_accuracy <- map_dfr(all_result, function(x){x}) %>% 
  filter(data_property != "libraryvscellzero")
spatial_accuracy <- spatial_accuracy %>% 
  group_by(Data, data_property) %>% 
  mutate(
    across(all_of("spatial_KDE"), ~ pnorm((.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
  )
spatial_accuracy <- spatial_accuracy %>% 
  group_by(Data, data_property) %>% 
  mutate(
    across(all_of(colnames(spatial_accuracy)[4:5]), ~ 1 - .x)
  )
### scale for every metric [0, 1]
spatial_accuracy <- spatial_accuracy %>% 
  group_by(Data, data_property) %>% 
  mutate(
    across(all_of(colnames(spatial_accuracy)[4:5]), ~ (.x - min(.x, na.rm = TRUE)) / (max(.x, na.rm = TRUE) - min(.x, na.rm = TRUE)))
  )
spatial_accuracy <- spatial_accuracy %>% 
  pivot_longer(cols = c(4:5), names_to = "metric", values_to = "value") %>% 
  rename(., property = data_property)
saveRDS(spatial_accuracy, file = "Chunk8-Data Analysis/2-accuracy/spatial_accuracy_data.rds")
