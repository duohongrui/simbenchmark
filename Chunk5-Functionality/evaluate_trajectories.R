library(tibble)
ref_data_list <- list.files("../preprocessed_data/")

for(i in ref_data_list){
  ref_data <- readRDS(file.path("../preprocessed_data", i))
  if(!dynwrap::is_data_wrapper(ref_data$data)){
    next
  }
  message(paste0(i, "..."))
  data_name <- stringr::str_split(i, pattern = "_", simplify = TRUE)[1]
  data_list <- list.files("../simulation_data/", pattern = data_name)
  
  for(w in data_list){
    message(paste0("--------", w, "", "--------"))
    sim_data <- readRDS(file.path("../simulation_data", w))
    
    if(data$sim_data_info$group >= 2 & "group" %in% colnames(data$sim_data$col_meta)){
      sim_data_grouping <- data$sim_data$col_meta$group
    }else{
      sim_data_grouping <- NULL
    }
    ### calculate trajectory metrics
    traj_result <- simutils::calculate_trajectory_properties(
      ref_data = ref_data,
      sim_data = sim_data$sim_data$count_data,
      sim_data_grouping = sim_data_grouping,
      seed = 666,
      verbose = TRUE
    )
    ### save results
    saveRDS(traj_result, file = paste0("../trajectory/", i))
    
    saveRDS(dplyr::lst(
      HIM = traj_result$HIM,
      F1_branches = traj_result$F1_branches,
      F1_milestones = traj_result$F1_milestones,
      Cor_dist = traj_result$Cor_dist$cor_dist$cor_dist
    ),
      file = paste0("../trajectory_evaluation/", i))
  }
}

# est <- simmethods::Splat_estimation(t(data$data$counts), seed = 111)
# sim <- simmethods::Splat_simulation(est$estimate_result,
#                                     return_format = "list",
#                                     seed = 111, other_prior = list(paths = TRUE))
# traj_result <- simutils::calculate_trajectory_properties(
#   ref_data = data$data,
#   sim_data = sim$simulate_result$count_data,
#   sim_data_grouping = NULL,
#   seed = 666,
#   verbose = TRUE
# )
