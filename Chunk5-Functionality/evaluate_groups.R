library(tibble)
data_list <- list.files("../simulation_data/")

for(i in data_list){
  data <- readRDS(file.path("../simulation_data", i))
  
  if(data$sim_data_info$group >= 2 & "group" %in% colnames(data$sim_data$col_meta) |
     data$sim_data_info$group >= 2 & "plate" %in% colnames(data$sim_data$col_meta)){
  
    message(i)
    
    if("plate" %in% colnames(data$sim_data$col_meta)){
      data$sim_data$col_meta$group <- data$sim_data$col_meta$plate
    }
    
    message("Calculating distance matrix...")
    dist <- parallelDist::parDist(t(as.matrix(data$sim_data$count_data)), threads = 8)
    
    message("1-Calculating CDI...")
    error <- try(
      CDI <- simutils::calculate_CDI(as.matrix(data$sim_data$count_data), cluster_info = data$sim_data$col_meta$group),
      silent = TRUE
    )
    if("try-error" %in% class(error)){
      CDI <- NA
    }else{
      CDI <- min(CDI[1, 1], CDI[1, 2])
    }
    
    message("2-Calculating ROUGE...")
    error <- try(
      ROUGE <- simutils::calculate_ROGUE(as.matrix(data$sim_data$count_data), cluster_info = data$sim_data$col_meta$group),
      silent = TRUE
    )
    if("try-error" %in% class(error)){
      ROUGE <- NA
    }
    
    message("3-Calculating silhouette...")
    error <- try(
      silhouette <- simutils::calculate_silhouette(dist, cluster_info = data$sim_data$col_meta$group),
      silent = TRUE
    )
    if("try-error" %in% class(error)){
      silhouette <- NA
    }
    
    message("4-Calculating dunn...")
    error <- try(
      dunn <- simutils::calculate_dunn(dist, cluster_info = data$sim_data$col_meta$group),
      silent = TRUE
    )
    if("try-error" %in% class(error)){
      dunn <- NA
    }
    
    message("5-Calculating connectivity...")
    error <- try(
      connectivity <- simutils::calculate_connectivity(dist, cluster_info = data$sim_data$col_meta$group),
      silent = TRUE
    )
    if("try-error" %in% class(error)){
      connectivity <- NA
    }
    
    message("6-Calculating DB index...")
    error <- try(
      DB_index <- simutils::calculate_DB_index(as.matrix(data$sim_data$count_data), cluster_info = data$sim_data$col_meta$group),
      silent = TRUE
    )
    if("try-error" %in% class(error)){
      DB_index <- NA
    }
    
    group_metrics <- dplyr::lst(CDI,
                                ROUGE,
                                silhouette,
                                dunn,
                                connectivity,
                                DB_index)
    
    saveRDS(group_metrics, file.path("../group_evaluation", i))
    
  }else{
    next
  }
}


                      