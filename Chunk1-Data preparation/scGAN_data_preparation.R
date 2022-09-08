data_list <- list.files("/Users/duohongrui/Desktop/preprocessed_data")
save_path <- "./scGAN_preprocessed_data"

for(i in 1:101){
  data <- readRDS(file.path("/Users/duohongrui/Desktop/preprocessed_data", data_list[i]))
  message("-------------------------")
  message(data_list[i])
  ## group
  data_info <- data$data_info
  group <- data_info$group_condition
  cluster_info <- data_info$cluster_info
  if(is.null(group)){
    if(!is.null(cluster_info)){
      group <- as.numeric(as.factor(cluster_info))
    }else{
      next
    }
  }
  group <- group-1
  data_id <- data_info$id
  if(dynwrap::is_wrapper_with_expression(data$data)){
    data <- t(data$data[["counts"]])
  }else{
    data <- data$data
  }
  result <- simutils::scgan_data_conversion(data = data,
                                            data_i = data_id,
                                            group = group,
                                            save_to_path = save_path,
                                            verbose = FALSE)
  message("Done")
  message("-------------------------")
  Sys.sleep(2)
}
