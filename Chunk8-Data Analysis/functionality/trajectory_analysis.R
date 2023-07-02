library(stringr)
library(purrr)
library(dplyr)
library(tidyr)


data_list <- list.files("../trajectory_evaluation/")
a <- readRDS("../trajectory_evaluation/dyngen_data100_stimulated-dendritic-cells-PAM_shalek.rds.rds")
metric_name <- names(a)


### read data into a list
all_result <- list()
for (i in data_list) {
  result <- readRDS(file.path("../trajectory_evaluation", i))
  all_result[[i]] <- result
}

### turn to a tibble
trajectory_data <- purrr::map_dfr(1:length(all_result), .f = function(index){
  result <- as.numeric(all_result[[index]])
  names(result) <- metric_name
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
    across(all_of(metric_name), as.numeric)
  )

saveRDS(trajectory_data, file = "Chunk8-Data Analysis/functionality/trajectory_data.rds")
### turn to long table
trajectory_long_data <- trajectory_data %>% 
  pivot_longer(., cols = 3:ncol(.), names_to = "metric", values_to = "value")
saveRDS(trajectory_long_data, file = "Chunk8-Data Analysis/functionality/trajectory_long_data.rds")

