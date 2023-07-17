library(stringr)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

data_list <- list.files("../scalability2")

### read data into a list
all_result <- list()
for (i in data_list) {
  result <- readRDS(file.path("../scalability2", i))
  all_result[[i]] <- result
}

### turn to a tibble
scalability_data <- purrr::map_dfr(1:length(all_result), .f = function(index){
  all_result[[index]]
})
scalability_data2 <- purrr::map_dfr(1:length(all_result), .f = function(index){
  all_result[[index]]
})
scalability_data <- rbind(scalability_data, scalability_data2)
rm(scalability_data2)
saveRDS(scalability_data, file = "./Chunk8-Data Analysis/scalability/scalability_data.rds")
### scalability score

#### aggregate by repeat time
score_data <- scalability_data %>% 
  group_by(method, cell_num, gene_num) %>% 
  summarise(
    estimation_time = mean(estimation_time, na.rm = TRUE),
    estimation_memory = mean(estimation_memory, na.rm = TRUE),
    simulation_time = mean(simulation_time, na.rm = TRUE),
    simulation_memory = mean(simulation_memory, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  group_by(cell_num, gene_num) %>% 
  mutate(
    across(all_of(c("estimation_time",
                    "estimation_memory",
                    "simulation_time",
                    "simulation_memory")), ~ 1 - pnorm((.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
  ) %>% 
  ungroup()
  
### scale for every metric [0, 1]
score_data <- score_data %>% 
  mutate(
    across(all_of(colnames(score_data)[4:7]), ~ (.x - min(.x, na.rm = TRUE)) / (max(.x, na.rm = TRUE) - min(.x, na.rm = TRUE)))
  )

#### aggregate by method
score_data <- score_data %>% 
  group_by(method) %>% 
  summarise(
    estimation_time = mean(estimation_time, na.rm = TRUE),
    estimation_memory = mean(estimation_memory, na.rm = TRUE),
    simulation_time = mean(simulation_time, na.rm = TRUE),
    simulation_memory = mean(simulation_memory, na.rm = TRUE)
  ) %>% 
  ungroup()
score_data[42, c(3,5)] <- NaN
### time and memory score
score_data <- score_data %>% 
  mutate(
    time_score = apply(.[, -1], MARGIN = 1, function(x){mean(c(x[1], x[3]), na.rm = TRUE)}),
    memory_score = apply(.[, -1], MARGIN = 1, function(x){mean(c(x[2], x[4]), na.rm = TRUE)}),
    scalability_score = apply(.[, -1], MARGIN = 1, function(x){mean(c(x[1], x[2], x[3], x[4]), na.rm = TRUE)})
  )

### correlation values derived from "Shape constrained additive models and Random Forest" section
scam_model <- readRDS("./Chunk8-Data Analysis/scalability/scam_model.rds")
all_methods_cor <- map_dfc(scam_model, .f = function(x){
  id_name <- names(x)[1]
  step <- str_split(id_name, "_", simplify = TRUE)[1]
  feature <- str_split(id_name, "_", simplify = TRUE)[2]
  method_name <- names(x[[1]])
  
  method_cor <- map_dfr(method_name, function(y){
    cor_list <- x[[3]][[3]]
    per_method <- cor_list[[y]]
    if(is.null(per_method)){
      value <- c(NA)
    }else{
      value <- c(per_method$cor_value)
    }
    names(value) <- c("cor")
    value
  })
  ## numeric correlation values
  method_cor <- method_cor %>% 
    mutate(
      cor = as.numeric(cor)
    )
  colnames(method_cor) <- paste0(step, "_", feature, "_cor")
  method_cor
})
all_methods_cor <- all_methods_cor %>% 
  mutate(
    method = names(scam_model[[1]][[1]])
  ) %>% 
  relocate(method, .before = estimation_time_cor)
score_data <- score_data %>% 
  full_join(all_methods_cor, by = "method")

saveRDS(score_data, file = "Chunk8-Data Analysis/scalability/score_data.rds")



################################################################################
#######################            Plot Rect Data      #########################
################################################################################
scalability_data <- readRDS("./Chunk8-Data Analysis/scalability/scalability_data.rds")

### scalability score

#### aggregate by repeat time
score_data <- scalability_data %>% 
  group_by(method, cell_num, gene_num) %>% 
  summarise(
    estimation_time = mean(estimation_time, na.rm = TRUE),
    estimation_memory = mean(estimation_memory, na.rm = TRUE),
    simulation_time = mean(simulation_time, na.rm = TRUE),
    simulation_memory = mean(simulation_memory, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  filter(cell_num == 100 & gene_num == 1000 |
         cell_num == 1000 & gene_num == 1000 |
         cell_num == 1000 & gene_num == 10000 |
         cell_num == 10000 & gene_num == 1000) %>% 
  mutate(
    estimation_time_text = case_when(
      estimation_time < 1 ~ "<1s",
      estimation_time >= 1 & estimation_time < 60 ~ paste0(as.character(round(estimation_time)), "s"),
      estimation_time >= 60 & estimation_time < 3600 ~ paste0(estimation_time %/% 60, "m"),
      estimation_time > 3600 ~ paste0(">", estimation_time %/% 3600, 'h')
    ),
    simulation_time_text = case_when(
      simulation_time < 1 ~ "<1s",
      simulation_time >= 1 & simulation_time < 60 ~ paste0(as.character(round(simulation_time)), "s"),
      simulation_time >= 60 & simulation_time < 3600 ~ paste0(simulation_time %/% 60, "m"),
      simulation_time > 3600 ~ paste0(">", simulation_time %/% 3600, 'h')
    ),
    estimation_memory_text = case_when(
      estimation_memory < 10 ~ paste0(round(estimation_memory, digits = 1), 'M'),
      estimation_memory >= 10 & estimation_memory < 1024 ~ paste0(round(estimation_memory), 'M'),
      estimation_memory >= 1024 ~ paste0(">", estimation_memory %/%1024, 'G')
    ),
    simulation_memory_text = case_when(
      simulation_memory < 10 ~ paste0(round(simulation_memory, digits = 1), 'M'),
      simulation_memory < 1024 ~ paste0(round(simulation_memory), 'M'),
      simulation_memory >= 1024 ~ paste0(">", simulation_memory %/%1024, 'G')
    )
  ) %>% 
  group_by(cell_num, gene_num) %>% 
  mutate(
    across(all_of(c("estimation_time",
                    "estimation_memory",
                    "simulation_time",
                    "simulation_memory")), ~ 1 - (pnorm((.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE))))
  ) %>% 
  ungroup() %>% 
  mutate(
    estimation_score = (estimation_time + estimation_memory) / 2,
    simulation_score = (simulation_time + simulation_memory) / 2,
    estimation_text = paste0(estimation_time_text, "\n", estimation_memory_text),
    estimation_text = case_when(
      estimation_text == "NA\nNA" ~ "",
      TRUE ~ estimation_text
    ),
    simulation_text = paste0(simulation_time_text, "\n", simulation_memory_text),
    simulation_text = case_when(
      simulation_text == "NA\nNA" ~ "",
      TRUE ~ simulation_text
    ),
    data = paste0(as.character(cell_num), "x", as.character(gene_num))
  )

data_100_1000_est <- score_data %>% 
  filter(cell_num == 100, gene_num == 1000) %>%
  select(method, estimation_score, estimation_text)
colnames(data_100_1000_est) <- c("method", paste0("scalability_100_1000_", colnames(data_100_1000_est)[2:3]))

data_100_1000_sim <- score_data %>% 
  filter(cell_num == 100, gene_num == 1000) %>%
  select(method, simulation_score, simulation_text)
colnames(data_100_1000_sim) <- c("method", paste0("scalability_100_1000_", colnames(data_100_1000_sim)[2:3]))

data_1000_1000_est <- score_data %>% 
  filter(cell_num == 1000, gene_num == 1000) %>%
  select(method, estimation_score, estimation_text)
colnames(data_1000_1000_est) <- c("method", paste0("scalability_1000_1000_", colnames(data_1000_1000_est)[2:3]))

data_1000_1000_sim <- score_data %>% 
  filter(cell_num == 1000, gene_num == 1000) %>%
  select(method, simulation_score, simulation_text)
colnames(data_1000_1000_sim) <- c("method", paste0("scalability_1000_1000_", colnames(data_1000_1000_sim)[2:3]))

data_1000_10000_est <- score_data %>% 
  filter(cell_num == 1000, gene_num == 10000) %>%
  select(method, estimation_score, estimation_text)
colnames(data_1000_10000_est) <- c("method", paste0("scalability_1000_10000_", colnames(data_1000_10000_est)[2:3]))

data_1000_10000_sim <- score_data %>% 
  filter(cell_num == 1000, gene_num == 10000) %>%
  select(method, simulation_score, simulation_text)
colnames(data_1000_10000_sim) <- c("method", paste0("scalability_1000_10000_", colnames(data_1000_10000_sim)[2:3]))

data_10000_1000_est <- score_data %>% 
  filter(cell_num == 10000, gene_num == 1000) %>%
  select(method, estimation_score, estimation_text)
colnames(data_10000_1000_est) <- c("method", paste0("scalability_10000_1000_", colnames(data_10000_1000_est)[2:3]))

data_10000_1000_sim <- score_data %>% 
  filter(cell_num == 10000, gene_num == 1000) %>%
  select(method, simulation_score, simulation_text)
colnames(data_10000_1000_sim) <- c("method", paste0("scalability_10000_1000_", colnames(data_10000_1000_sim)[2:3]))

score_data <- data_100_1000_est %>% 
  full_join(data_100_1000_sim, by = "method")%>% 
  full_join(data_1000_1000_est, by = "method")%>% 
  full_join(data_1000_1000_sim, by = "method")%>% 
  full_join(data_1000_10000_est, by = "method")%>% 
  full_join(data_1000_10000_sim, by = "method")%>% 
  full_join(data_10000_1000_est, by = "method")%>% 
  full_join(data_10000_1000_sim, by = "method")

saveRDS(score_data, file = "./Chunk8-Data Analysis/scalability/scalability_plot_data.rds")
