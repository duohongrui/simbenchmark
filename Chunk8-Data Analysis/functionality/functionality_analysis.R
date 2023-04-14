library(stringr)
library(purrr)
library(dplyr)
library(tidyr)
group_data <- readRDS("Chunk8-Data Analysis/functionality/group_data.rds")
batch_data <- readRDS("Chunk8-Data Analysis/functionality/batch_data.rds")
DEGs_data <- readRDS("Chunk8-Data Analysis/functionality/DEGs_data.rds")


functionality_data <- group_data %>% 
  full_join(., batch_data, by = c("Method", "Data")) %>% 
  full_join(., DEGs_data, by = c("Method", "Data"))

# ### Fill NA and NaN with 0
# functionality_data <- functionality_data %>% 
#   mutate(
#     across(3:ncol(.), ~ replace_na(.x, 0))
#   )
saveRDS(functionality_data, file = "Chunk8-Data Analysis/functionality/functionality_data.rds")

### turn to a long tibble
functionality_long_data <- functionality_data %>% 
  pivot_longer(., cols = 3:ncol(.), names_to = "metric", values_to = "value")
saveRDS(functionality_long_data, file = "Chunk8-Data Analysis/functionality/functionality_long_data.rds")
