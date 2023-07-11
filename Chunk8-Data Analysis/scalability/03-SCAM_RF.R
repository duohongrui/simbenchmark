################################################################################
##########    Shape constrained additive models and Random Forest   ############
################################################################################
data_list <- list.files("../scalability/")

### read data into a list
all_result <- list()
for (i in data_list) {
  result <- readRDS(file.path("../scalability", i))
  all_result[[i]] <- result
}

### turn to a tibble
scalability_data <- purrr::map_dfr(1:length(all_result), .f = function(index){
  all_result[[index]]
})

# #### aggregate by repeat time
# score_data <- scalability_data %>% 
#   group_by(method, cell_num, gene_num) %>% 
#   summarise(
#     estimation_time = mean(estimation_time, na.rm = TRUE),
#     estimation_memory = mean(estimation_memory, na.rm = TRUE),
#     simulation_time = mean(simulation_time, na.rm = TRUE),
#     simulation_memory = mean(simulation_memory, na.rm = TRUE)
#   ) %>% 
#   ungroup()

### execution data
scam_data_list <- list.files("/Volumes/Elements/sim_bench/simulation_data/")
scam_data <- map_dfr(1:length(scam_data_list), .f = function(x){
  print(x)
  data_tmp <- readRDS(paste0("/Volumes/Elements/sim_bench/simulation_data/", scam_data_list[x]))
  method_name_split <- str_split(scam_data_list[x], "_", simplify = TRUE)
  method_name <- method_name_split[1]
  if(method_name_split[2] == "ERCC"){
    data_id <- method_name_split[3]
  }else{
    data_id <- method_name_split[2]
  }
  data <- c(method_name,
            data_id,
            data_tmp$sim_data_info$cell_num,
            data_tmp$sim_data_info$gene_num,
            data_tmp$sim_data_info$estimate_time,
            data_tmp$sim_data_info$estimate_memory,
            data_tmp$sim_data_info$simulation_time,
            data_tmp$sim_data_info$simulation_memory)
  names(data) <- c("method",
                   "data_id",
                   "cell_num",
                   "gene_num",
                   "estimate_time",
                   "estimate_memory",
                   "simulation_time",
                   "simulation_memory")
  data
})
scam_data$cell_num <- as.numeric(scam_data$cell_num)
scam_data$gene_num <- as.numeric(scam_data$gene_num)
scam_data$estimate_time <- as.numeric(scam_data$estimate_time)
scam_data$estimate_memory <- as.numeric(scam_data$estimate_memory)
scam_data$simulation_time <- as.numeric(scam_data$simulation_time)
scam_data$simulation_memory <- as.numeric(scam_data$simulation_memory)
saveRDS(scam_data, file = "./Chunk8-Data Analysis/scalability/execution_detection_data.rds")

colnames(scam_data) <- c("method",
                         "data_id",
                         "cell_num",
                         "gene_num",
                         "estimation_time",
                         "estimation_memory",
                         "simulation_time",
                         "simulation_memory")
scam_data[c(1,2), c(6,8)] <- NA
scam_data <- scam_data %>% 
  select(-2) %>% 
  bind_rows(., scalability_data %>% select(-4)) %>% 
  mutate(
    row = row_number()
  )
scam_data <- scam_data %>% 
  filter(estimation_memory != 0 | is.na(estimation_memory))
scam_data[scam_data == 0] <- 0.01

# scam_data <- scam_data %>% 
#   select(-2) %>% 
#   bind_rows(., score_data) %>% 
#   mutate(
#     row = row_number()
#   )
### scGAN NA for memory
scam_data[scam_data$method == "scGAN", "estimation_memory"] <- NA
scam_data[scam_data$method == "scGAN", "simulation_memory"] <- NA


#### split data
set.seed(65)
train_data <- scam_data %>% 
  group_by(method) %>% 
  sample_frac(0.8) %>% 
  ungroup()
test_data <- scam_data[scam_data$row %in% setdiff(1:nrow(scam_data), train_data$row), ]
saveRDS(train_data, file = "./Chunk8-Data Analysis/scalability/scam_train_data.rds")
saveRDS(test_data, file = "./Chunk8-Data Analysis/scalability/scam_test_data.rds")

########################## scam
#### time of estimation
estimation_time_model <- scam_function(step = "estimation", feature = "time")
est_time_table <- bind_prediction_result(estimation_time_model)
est_time_correlation <- correlation_plot(est_time_table)
est_time_correlation$correlation_plot
#### memory of estimation
estimation_memory_model <- scam_function(step = "estimation", feature = "memory")
est_memory_table <- bind_prediction_result(estimation_memory_model)
est_memory_correlation <- correlation_plot(est_memory_table)
est_memory_correlation$correlation_plot
#### time of simulation
simulation_time_model <- scam_function(step = "simulation", feature = "time")
sim_time_table <- bind_prediction_result(simulation_time_model)
sim_time_correlation <- correlation_plot(sim_time_table)
sim_time_correlation$correlation_plot
#### memory of simulation
simulation_memory_model <- scam_function(step = "simulation", feature = "memory")
sim_memory_table <- bind_prediction_result(simulation_memory_model)
sim_memory_correlation <- correlation_plot(sim_memory_table)
sim_memory_correlation$correlation_plot
scam_model <- list("estimation_time_model" = list("estimation_time_model" = estimation_time_model,
                                                  "est_time_table" = est_time_table,
                                                  "est_time_correlation" = est_time_correlation),
                   "estimation_memory_model" = list("estimation_memory_model" = estimation_memory_model,
                                                    "est_memory_table" = est_memory_table,
                                                    "est_memory_correlation" = est_memory_correlation),
                   "simulation_time_model" = list("simulation_time_model" = simulation_time_model,
                                                  "sim_time_table" = sim_time_table,
                                                  "sim_time_correlation" = sim_time_correlation),
                   "simulation_memory_model" = list("simulation_memory_model" = simulation_memory_model,
                                                    "sim_memory_table" = sim_memory_table,
                                                    "sim_memory_correlation" = sim_memory_correlation))
saveRDS(scam_model, file = "./Chunk8-Data Analysis/scalability/scam_model.rds")

ggsave(plot = wrap_plots(map(scam_model, function(x){x[[3]][[2]]}), ncol = 2),
       filename = "./Chunk8-Data Analysis/scalability/scam_model.pdf",
       height = 8,
       width = 8,
       units = "in")

########################## RF
#### time of estimation
estimation_time_model <- scam_function(step = "estimation", feature = "time", model_selected = "RF")
est_time_table <- bind_prediction_result(estimation_time_model)
est_time_correlation <- correlation_plot(est_time_table)
est_time_correlation$correlation_plot
#### memory of estimation
estimation_memory_model <- scam_function(step = "estimation", feature = "memory", model_selected = "RF")
est_memory_table <- bind_prediction_result(estimation_memory_model)
est_memory_correlation <- correlation_plot(est_memory_table)
est_memory_correlation$correlation_plot
#### time of simulation
simulation_time_model <- scam_function(step = "simulation", feature = "time", model_selected = "RF")
sim_time_table <- bind_prediction_result(simulation_time_model)
sim_time_correlation <- correlation_plot(sim_time_table)
sim_time_correlation$correlation_plot
#### memory of simulation
simulation_memory_model <- scam_function(step = "simulation", feature = "memory", model = "RF")
sim_memory_table <- bind_prediction_result(simulation_memory_model)
sim_memory_correlation <- correlation_plot(sim_memory_table)
sim_memory_correlation$correlation_plot
RF_model <- list("estimation_time_model" = list("estimation_time_model" = estimation_time_model,
                                                "est_time_table" = est_time_table,
                                                "est_time_correlation" = est_time_correlation),
                 "estimation_memory_model" = list("estimation_memory_model" = estimation_memory_model,
                                                  "est_memory_table" = est_memory_table,
                                                  "est_memory_correlation" = est_memory_correlation),
                 "simulation_time_model" = list("simulation_time_model" = simulation_time_model,
                                                "sim_time_table" = sim_time_table,
                                                "sim_time_correlation" = sim_time_correlation),
                 "simulation_memory_model" = list("simulation_memory_model" = simulation_memory_model,
                                                  "sim_memory_table" = sim_memory_table,
                                                  "sim_memory_correlation" = sim_memory_correlation))
saveRDS(RF_model, file = "./Chunk8-Data Analysis/scalability/RF_model.rds")
ggsave(plot = wrap_plots(map(RF_model, function(x){x[[3]][[2]]}), ncol = 2),
       filename = "./Chunk8-Data Analysis/scalability/RF_model.pdf",
       height = 8,
       width = 8,
       units = "in")