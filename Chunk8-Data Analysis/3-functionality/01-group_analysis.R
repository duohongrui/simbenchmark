library(stringr)
library(purrr)
library(dplyr)
library(tidyr)


data_list <- list.files("../group_evaluation/")
a <- readRDS("../group_evaluation/ESCO_data102_spatial_B1.rds")
metric_name <- names(a)
rm(a)

### read data into a list
all_result <- list()
for (i in data_list) {
  result <- readRDS(file.path("../group_evaluation", i))
  all_result[[i]] <- result
}

### turn to a tibble
group_data <- purrr::map_dfr(1:length(all_result), .f = function(index){
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

### Filter and save Inf in CDI metric
inf_data <- group_data %>% 
  filter(is.infinite(CDI))
saveRDS(inf_data, file = "Chunk8-Data Analysis/3-functionality/inf_data_group.rds")

### normalize some values
group_data <- group_data %>% 
  mutate(
    across(all_of("CDI"), ~ replace(.x, .x == "Inf", values = NA))
  ) %>% 
  group_by(Data) %>% 
  mutate(
    across(all_of(metric_name[-2]), ~ pnorm((.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
  ) %>% 
  ungroup()

### Subtract values by 1
colume_name <- c("CDI", "connectivity", "DB_index")
group_data <- group_data %>% 
  mutate(
    across(all_of(colume_name), ~ 1 - .x)
  )
### scale for every metric [0, 1]
group_data <- group_data %>% 
  mutate(
    across(all_of(colnames(group_data)[3:8]), ~ (.x - min(.x, na.rm = TRUE)) / (max(.x, na.rm = TRUE) - min(.x, na.rm = TRUE)))
  )

saveRDS(group_data, file = "Chunk8-Data Analysis/3-functionality/group_data.rds")
### turn to long table
group_long_data <- group_data %>% 
  pivot_longer(., cols = 3:ncol(.), names_to = "metric", values_to = "value")
saveRDS(group_long_data, file = "Chunk8-Data Analysis/3-functionality/group_long_data.rds")

