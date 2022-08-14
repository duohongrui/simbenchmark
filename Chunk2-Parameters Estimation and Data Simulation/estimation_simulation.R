## data list
data_list <- list.files("/Users/duohongrui/Desktop/preprocessed_data/")

for(i in 1:length(data_list)){
  file_name <- data_list[i]
  ## read data
  data <- readRDS(file.path("/Users/duohongrui/Desktop/preprocessed_data", file_name))
  data_info <- data$data_info
  if(dynwrap::is_wrapper_with_expression(data$data)){
    counts <- data$data[["counts"]]
  }else{
    counts <- data$data
  }
  ## data meta info
  ERCC <- data_info$ERCC
  dilution <- data_info$dilution_factor
  volumn <- data_info$volume
  group <- data_info$group_condition
  treatment <- data_info$treatment
  batch_info <- data_info$batch_info
  cluster_info <- data_info$cluster_info
  
  ## estimation
  estimation_result <- simpipe::estimate_parameters(
    ref_data = counts,
    method = "Splat",
    other_prior = list(group_condition = group),
    seed = 111,
    verbose = TRUE,
    use_docker = FALSE
  )
  
}


