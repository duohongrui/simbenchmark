library(stringr)
library(purrr)
library(dplyr)
library(tidyr)


data_list <- list.files("../DEGs_evaluation/")
a <- readRDS("../DEGs_evaluation/powsimR_data22_GSE86469.rds")
metric_name <- names(a)


### read data into a list
all_result <- list()
for (i in data_list) {
  result <- readRDS(file.path("../DEGs_evaluation", i))
  if(is.null(result[["distribution_score"]])){
    result[["distribution_score"]] <- NA
  }
  all_result[[i]] <- result
}

### turn to a tibble
DEGs_data <- purrr::map_dfr(1:length(all_result), .f = function(index){
  result <- as.numeric(all_result[[index]])
  if(length(result) == 6){
    result <- c(result, NA)
  }
  if(length(result) > 7){
    result <- result[1:7]
  }
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

saveRDS(DEGs_data, file = "Chunk8-Data Analysis/3-functionality/DEGs_data.rds")
### turn to long table
DEGs_long_data <- DEGs_data %>% 
  pivot_longer(., cols = 3:ncol(.), names_to = "metric", values_to = "value")
saveRDS(DEGs_long_data, file = "Chunk8-Data Analysis/3-functionality/DEGs_long_data.rds")

