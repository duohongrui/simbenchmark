library(tibble)
data_list <- list.files("F:/sim_bench/simulation_data/", pattern = "^SPARSim_")

for(i in data_list){
  data <- readRDS(file.path("F:/sim_bench/simulation_data", i))
  if(data$sim_data_info$batch >= 2 & "batch" %in% colnames(data$sim_data$col_meta)){
    message(i)
    batch_info <- as.character(data$sim_data$col_meta$batch)
    batch_info <- as.numeric(as.factor(batch_info))
    batch_matrics <- simutils::calculate_batch_properties(data$sim_data$count,
                                                          batch_info = batch_info,
                                                          verbose = T)
    saveRDS(batch_matrics, file.path("F:/sim_bench/batch_evaluation", i))
    
  }else{
    next
  }
}


