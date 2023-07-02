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

### normalize values which are not in [0, 1] for every dataset
check_range <- purrr::map_dfr(.x = metrics_property, .f = function(x){
  range_check <- as.numeric(range(accuracy_data %>% pull(x), na.rm = TRUE))
  names(range_check) <- c("low", "high")
  range_check
})
rownames(check_range) <- metrics_property
normalized_columes <- colnames(accuracy_data)[grep(colnames(accuracy_data), pattern = "MAD|MAE|RMSE|bhattacharyya|KDE")]
# normalized_columes <- c("MAD_LS", "MAD_TMM", "MAD_ELS", "MAD_FCO", "MAD_CCC",
#                         "MAE_LS", "MAE_TMM", "MAE_ELS", "MAE_FCO",
#                         "RMSE_LS", "RMSE_TMM", "RMSE_ELS", "RMSE_FCO",
#                         "bhattacharyya_LS", "bhattacharyya_FZC", "bhattacharyya_CCC", "bhattacharyya_TMM", "bhattacharyya_ELS",
#                         "KDE_RLZ",
#                         "MAD_ME", "MAD_SD", "MAD_CV", "MAD_DIS",
#                         "MAE_ME", "MAE_SD", "MAE_CV", "MAE_DIS",
#                         "RMSE_ME", "RMSE_SD", "RMSE_CV", "RMSE_DIS",
#                         "bhattacharyya_ME", "bhattacharyya_SD", "bhattacharyya_CV", "bhattacharyya_FZG", "bhattacharyya_FZG", "bhattacharyya_DIS",
#                         "KDE_RMS", "KDE_RMZ", "KDE_RDM")

accuracy_data <- accuracy_data %>% 
  group_by(Data) %>% 
  mutate(
    across(all_of(normalized_columes), ~ pnorm((.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
  ) %>% 
  ungroup()

### Subtract values by 1
subtract_columes_prefix <- c("MAD", "KS", "MAE", "RMSE", "bhattacharyya", "multiKS")
subtract_columes <- colnames(accuracy_data)[which(str_split(colnames(accuracy_data),
                                                            pattern = "_", simplify = T)[,1] %in% subtract_columes_prefix)]
accuracy_data <- accuracy_data %>% 
  mutate(
    across(all_of(subtract_columes), ~ 1 - .x)
  )


# ### NA and NaN
# accuracy_data <- accuracy_data %>% 
#   mutate(
#     across(3:76, ~ replace_na(.x, 0))
#   )
saveRDS(accuracy_data, file = "Chunk8-Data Analysis/accuracy/accuracy_data.rds")
### turn to long table
accuracy_long_data <- accuracy_data %>% 
  pivot_longer(., cols = 3:ncol(.), names_to = "metric_property", values_to = "value")

### split property to a single variable
accuracy_long_data <- accuracy_long_data %>% 
  separate(., col = metric_property, sep = "_", into = c("metric", "property"))
saveRDS(accuracy_long_data, file = "Chunk8-Data Analysis/accuracy/accuracy_long_data.rds")
