#------------------------------------------------------------------------------#
# This file contains one simulation methods:
# 1. scGAN
#------------------------------------------------------------------------------#
library(simpipe)
library(simmethods)
library(dplyr)
library(stringr)

## data list
data_list <- list.files("../preprocessed_data/")

method <- "scGAN"

for(i in 1:length(data_list)){
  file_name <- data_list[i]
  ## read data
  message("-------------------------------------------------------------------")
  message(file_name)
  message("Read data...")
  data <- readRDS(file.path("../preprocessed_data", file_name))
  data_info <- data$data_info
  if(dynwrap::is_wrapper_with_expression(data$data)){
    counts <- t(data$data[["counts"]])
  }else{
    next
  }
  ## data meta info
  data_id <- data_info$id
  group <- data_info$group_condition
  treatment <- data_info$treatment
  cluster_info <- data_info$cluster_info
  
  ## 1) group
  if(is.null(group)){
    if(!is.null(cluster_info)){
      group <- as.numeric(as.factor(cluster_info))
    }
  }
  
  if(is.null(group)){
    next
  }
  
  if(dynwrap::is_wrapper_with_expression(data$data)){
    valid_cells <- round(nrow(data$data$counts) * 0.2)
  }else{
    valid_cells <- round(ncol(data$data) * 0.2)
  }
  
  ## data save file
  save_name <- paste0(method, "_", data_id)
  message(save_name)
  
  ## estimation (SCRIP only uses splatEstimate function to estimate parameters from real data)
  message("Estimating...")
  try_result <- try(
    estimation_result <- simpipe::estimate_parameters(
      ref_data = counts,
      method = method,
      other_prior = list(min_genes = 0,
                         min_cells = 0,
                         GPU = 1,
                         res = NULL,
                         valid_cells = valid_cells,
                         max_steps = 100000),
      seed = 1,
      verbose = TRUE,
      use_docker = FALSE),
    silent = FALSE,
    outFile = paste0("../error_text/", save_name, "_", "estimation_error.txt"))
  
  if("try-error" %in% class(try_result)){
    next
  }else{
    ### save
    saveRDS(estimation_result, paste0("../estimation_result/", save_name, "_estimation_result.rds"))
  }
  
  ## simulation
  message("Simulation...")
  try_result <- try(
    simulation_result <- simpipe::simulate_datasets(
      parameters = estimation_result,
      n = 1,
      seed = 1,
      return_format = "list",
      verbose = TRUE,
      use_docker = FALSE),
    silent = FALSE,
    outFile = paste0("../error_text/", save_name, "_", "simulation_error.txt"))
  
  if("try-error" %in% class(try_result)){
    next
  }
  
  ## estimation step monitor
  message("Save information of simulated data...")
  est_time <- estimation_result[[1]][["estimate_detection"]][1, 2]
  est_peak_memory <- estimation_result[[1]][["estimate_detection"]][1, 4]
  
  ## simulation step monitor
  sim_time <- simulation_result[[1]][["simulate_detection"]][1, 2]
  sim_peak_memory <- simulation_result[[1]][["simulate_detection"]][1, 4]
  
  ## simulated data info
  sim_counts <- simulation_result[[1]][["simulate_result"]][["count_data"]]
  ## simulated cell info
  ### group
  sim_col_data <- simulation_result[[1]][["simulate_result"]][["col_meta"]]
  if("group" %in% colnames(sim_col_data)){
    group_num <- length(unique(sim_col_data$group))
  }else{
    group_num <- "Not known"
  }

  sim_data_info <- list(sim_data_id = save_name,
                        method = method,
                        ref_data_platform = data_info$platform,
                        cell_num = ncol(sim_counts),
                        gene_num = nrow(sim_counts),
                        group = group_num,
                        batch = 0,
                        de_gene_num = 0,
                        estimate_time = est_time,
                        estimate_memory = est_peak_memory,
                        simulation_time = sim_time,
                        simulation_memory = sim_peak_memory)
  message("Save data...")
  save_result <- list(sim_data = simulation_result[[1]][["simulate_result"]],
                      sim_data_info = sim_data_info)
  ### save
  saveRDS(save_result, paste0("../simulation_data/", save_name, ".rds"))
  message("Done...")
  message("-------------------------------------------------------------------")
}
