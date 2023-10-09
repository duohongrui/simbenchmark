library(stringr)
library(purrr)
library(dplyr)
library(tidyr)


data_list <- list.files("../batch_evaluation/")
a <- readRDS("../batch_evaluation/SCRIP-BGP-commonBCV_data11_Figshare_Aorta.rds")
metric_name <- names(a)[-c(8, 9)]


### read data into a list
all_result <- list()
for (i in data_list) {
  result <- readRDS(file.path("../batch_evaluation", i))
  name <- names(result)[-c(8, 9)]
  result <- map(c(1:7), .f = function(x){
    mean(result[[x]])
  }) %>% setNames(name)
  all_result[[i]] <- result
}

### turn to a tibble
batch_data <- purrr::map_dfr(1:length(all_result), .f = function(index){
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

### normalize some values
#### cms, MM, LISI, pcr
batch_data <- batch_data %>% 
  group_by(Data) %>% 
  mutate(
    across(all_of(c("cms", "mm", "LISI", "pcr")), ~ pnorm((.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
  ) %>% 
  ungroup()
#### AWS_batch
batch_data <- batch_data %>% 
  mutate(
    across(AWS_batch, abs)
  )

### Subtract values by 1
colume_name <- c("cms", "LISI", "shannon_entropy")
batch_data <- batch_data %>% 
  mutate(
    across(all_of(colume_name), ~ 1 - .x)
  )
### scale for every metric [0, 1]
batch_data <- batch_data %>% 
  mutate(
    across(all_of(colnames(batch_data)[3:9]), ~ (.x - min(.x, na.rm = TRUE)) / (max(.x, na.rm = TRUE) - min(.x, na.rm = TRUE)))
  )
saveRDS(batch_data, file = "Chunk8-Data Analysis/3-functionality/batch_data.rds")
### turn to long table
batch_long_data <- batch_data %>% 
  pivot_longer(., cols = 3:ncol(.), names_to = "metric", values_to = "value")
saveRDS(batch_long_data, file = "Chunk8-Data Analysis/3-functionality/batch_long_data.rds")
