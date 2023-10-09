library(stringr)
library(purrr)
library(dplyr)
library(tidyr)
group_data <- readRDS("Chunk8-Data Analysis/3-functionality/group_data.rds")
DEGs_data <- readRDS("Chunk8-Data Analysis/3-functionality/DEGs_data.rds")
batch_data <- readRDS("Chunk8-Data Analysis/3-functionality/batch_data.rds")
trajectory_data <- readRDS("Chunk8-Data Analysis/3-functionality/trajectory_data.rds")


functionality_data <- group_data %>% 
  full_join(., DEGs_data, by = c("Method", "Data")) %>% 
  full_join(., batch_data, by = c("Method", "Data")) %>% 
  full_join(., trajectory_data, by = c("Method", "Data"))

saveRDS(functionality_data, file = "Chunk8-Data Analysis/3-functionality/functionality_data.rds")

### turn to a long tibble
functionality_long_data <- functionality_data %>% 
  pivot_longer(., cols = 3:ncol(.), names_to = "metric", values_to = "value")
saveRDS(functionality_long_data, file = "Chunk8-Data Analysis/3-functionality/functionality_long_data.rds")
