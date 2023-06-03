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
    if(file.exists(paste0("../trajectory_evaluation/", w))){
      next
    }
    if(!stringr::str_starts(w,
                            pattern = "PROSSTT|TedSim|dyntoy|dyngen|SymSim|VeloSim|MFA|Splat-paths|SplatPop-paths|phenopath|ESCO-traj|ESCO-tree|SCRIP-paths")){
      next
    }
    message(paste0("--------", w, "", "--------"))
    sim_data <- readRDS(file.path("../simulation_data", w))
    
    ### add fake cells for TedSim
    if(stringr::str_detect(w, pattern = "TedSim")){
      message("The number of cells is not the power of 2, and we will synthesize some extra cells base on your data...")
      ref_data <- simutils::synthesize_cells(t(as.matrix(ref_data$data$counts)),
                                             group = ref_data$data_info$cluster_info,
                                             seed = 1,
                                             verbose = TRUE)
      message("Performing trajectory inference by Slingshot for reference data...")
      ref_model <- dynwrap::infer_trajectory(dataset = ref_data,
                                             method = tislingshot::ti_slingshot(),
                                             parameters = NULL,
                                             give_priors = NULL,
                                             seed = 1,
                                             verbose = TRUE)
      ref_model <- dynwrap::add_expression(ref_model,
                                           counts = ref_data$counts,
                                           expression = ref_data$expression)
      ref_data <- list(data = ref_model)
    }
    
    if(sim_data$sim_data_info$group >= 2 & "group" %in% colnames(sim_data$sim_data$col_meta)){
      sim_data_grouping <- sim_data$sim_data$col_meta$group
    }else{
      sim_data_grouping <- NULL
    }
    ### calculate trajectory metrics
    error <- try(
      traj_result <- simutils::calculate_trajectory_properties(
        ref_data = ref_data$data,
        sim_data = as.matrix(sim_data$sim_data$count_data),
        sim_data_grouping = sim_data_grouping,
        seed = 666,
        verbose = TRUE
      ),
      silent = TRUE,
      outFile = paste0("../trajectory_error/", w, "_traj_error.txt")
    )
    if(is(error, "try-error")){
      next
    }
    ### save results
    saveRDS(traj_result, file = paste0("../trajectory/", w))
    
    saveRDS(dplyr::lst(
      HIM = traj_result$HIM,
      F1_branches = traj_result$F1_branches,
      F1_milestones = traj_result$F1_milestones,
      Cor_dist = traj_result$Cor_dist$cor_dist$cor_dist
    ),
      file = paste0("../trajectory_evaluation/", w))
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
