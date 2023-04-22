library(tibble)
data_list <- list.files("../simulation_data/")

for(i in data_list){
  data <- readRDS(file.path("../simulation_data", i))
  if(data$sim_data_info$batch >= 2 & "batch" %in% colnames(data$sim_data$col_meta)){
    message(i)
    batch_info <- as.character(data$sim_data$col_meta$batch)
    k <- round(sqrt(ncol(data$sim_data$count)))
    message(k)
    batch_matrics <- simutils::calculate_batch_properties(data = as.matrix(data$sim_data$count),
                                                          batch_info = batch_info,
                                                          verbose = T,
                                                          k = k)
    saveRDS(batch_matrics, file.path("../batch_evaluation", i))
    
  }else{
    next
  }
}


