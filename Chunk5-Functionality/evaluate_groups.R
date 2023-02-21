library(tibble)
data_list <- list.files()

for(i in data_list){
  data <- readRDS()
  
  if(data$sim_data_info$group >= 2 & "group" %in% colnames(data$sim_data$col_meta) |
     data$sim_data_info$group >= 2 & "plate" %in% colnames(data$sim_data$col_meta)){
  
    message(i)
    
    if("plate" %in% colnames(data$sim_data$col_meta)){
      data$sim_data$col_meta$group <- data$sim_data$col_meta$plate
    }
    
    message("Calculating distance matrix...")
    dist <- parallelDist::parDist(t(data$sim_data$count_data), threads = 5)
    
    message("1-Calculating CDI...")
    CDI <- simutils::calculate_CDI(data$sim_data$count_data, cluster_info = data$sim_data$col_meta$group)
    CDI <- min(CDI[1, 1], CDI[1, 2])
    
    message("2-Calculating ROUGE...")
    ROUGE <- simutils::calculate_ROGUE(data$sim_data$count_data, cluster_info = data$sim_data$col_meta$group)
    
    message("3-Calculating silhouette...")
    silhouette <- simutils::calculate_silhouette(dist, cluster_info = data$sim_data$col_meta$group)
    
    message("4-Calculating dunn...")
    dunn <- simutils::calculate_dunn(dist, cluster_info = data$sim_data$col_meta$group)
    
    message("5-Calculating connectivity...")
    connectivity <- simutils::calculate_connectivity(dist, cluster_info = data$sim_data$col_meta$group)
    
    message("6-Calculating DB index...")
    DB_index <- simutils::calculate_DB_index(data$sim_data$count_data, cluster_info = data$sim_data$col_meta$group)
    
    group_metrics <- dplyr::lst(dist,
                                CDI,
                                ROUGE,
                                silhouette,
                                dunn,
                                connectivity,
                                DB_index)
    
  }else{
    next
  }
}


                      