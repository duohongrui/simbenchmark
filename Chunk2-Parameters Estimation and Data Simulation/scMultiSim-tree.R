#------------------------------------------------------------------------------#
# This file contains one simulation methods:
# 1. scMultiSim-tree
#------------------------------------------------------------------------------#
library(simpipe)
library(simmethods)
library(dplyr)
library(stringr)

## data list
data_list <- list.files("../preprocessed_data/")

method <- "scMultiSim-tree"

for(i in 1:length(data_list)){
  file_name <- data_list[i]
  ## read data
  message("-------------------------------------------------------------------")
  message(file_name)
  message("Read data...")
  data <- readRDS(file.path("../preprocessed_data", file_name))
  data_info <- data$data_info
  if(dynwrap::is_wrapper_with_expression(data$data)){
    counts <- t(as.matrix(data$data[["counts"]]))
  }else{
    next
  }
  ## data meta info
  data_id <- data_info$id
  ERCC <- data_info$ERCC
  dilution <- data_info$dilution_factor
  volumn <- data_info$volume
  group <- data_info$group_condition
  treatment <- data_info$treatment
  batch_info <- data_info$batch_info
  cluster_info <- data_info$cluster_info
  
  ## 1) group
  if(is.null(group)){
    if(!is.null(cluster_info)){
      group <- as.numeric(as.factor(cluster_info))
    }
  }
  group.prob <- as.numeric(table(group)/length(group))
  
  traj <- stringr::str_split(method, pattern = "-", simplify = TRUE)[2]
  method_id <- stringr::str_split(method, pattern = "-", simplify = TRUE)[1]
  
  ## data save file
  save_name <- paste0(method, "_", data_id)
  message(save_name)
  
  ## estimation (SCRIP only uses splatEstimate function to estimate parameters from real data)
  message("Estimating...")
  try_result <- try(
    estimation_result <- simpipe::estimate_parameters(
      ref_data = counts,
      method = method_id,
      other_prior = list(group.condition = group),
      seed = 1,
      verbose = TRUE),
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
      ref_data = counts,
      other_prior = list(traj = TRUE,
                         prob.group = group.prob),
      n = 1,
      seed = 1,
      return_format = "list",
      verbose = TRUE),
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
  ## simulated gene info
  sim_row_data <- simulation_result[[1]][["simulate_result"]][["row_meta"]]
  ### DEGs
  if("de_gene" %in% colnames(sim_row_data)){
    de_gene_num <- sum(sim_row_data$de_gene == "yes")
  }else{
    de_gene_num <- 0
  }
  sim_data_info <- list(sim_data_id = save_name,
                        method = method,
                        ref_data_platform = data_info$platform,
                        cell_num = ncol(sim_counts),
                        gene_num = nrow(sim_counts),
                        group = group_num,
                        batch = 0,
                        de_gene_num = de_gene_num,
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
