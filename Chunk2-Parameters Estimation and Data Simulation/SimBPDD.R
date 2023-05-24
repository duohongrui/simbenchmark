#------------------------------------------------------------------------------#
# This file contains one simulation method:
# 1. SimBPDD
#------------------------------------------------------------------------------#
library(simpipe)
library(dplyr)
library(stringr)

## data list
data_list <- list.files("../preprocessed_data/")

method <- "SimBPDD"

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
    counts <- data$data
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
  if(is.null(group)){
    next
  }
  if(length(unique(group)) != 2){
    next
  }
  prob.group <- as.numeric(table(group)[1])/length(group)
  prob.group <- append(prob.group, 1 - prob.group)
  
  ## data save file
  save_name <- paste0(method, "_", data_id)
  message(save_name)
  
  ## simulation
  message("Simulation...")
  try_result <- try(
    simulation_result <- simpipe::simulate_datasets(
      ref_data = counts,
      method = method,
      other_prior = list(prob.group = prob.group,
                         degree = 1/3),
      n = 1,
      seed = 1,
      return_format = "list",
      verbose = TRUE),
    silent = FALSE,
    outFile = paste0("../error_text/", save_name, "_", "simulation_error.txt"))
  
  if("try-error" %in% class(try_result)){
    next
  }
  
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
  } else group_num <- 1
  ### batch
  batch_num <- 1
  ## simulated gene info
  sim_row_data <- simulation_result[[1]][["simulate_result"]][["row_meta"]]
  ### DEGs
  de_gene_num <- "Not known"
  sim_data_info <- list(sim_data_id = save_name,
                        method = method,
                        ref_data_platform = data_info$platform,
                        cell_num = ncol(sim_counts),
                        gene_num = nrow(sim_counts),
                        group = group_num,
                        batch = batch_num,
                        de_gene_num = de_gene_num,
                        estimate_time = NA,
                        estimate_memory = NA,
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
