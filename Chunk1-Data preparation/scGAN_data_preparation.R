data_list <- list.files("/Users/duohongrui/Desktop/preprocessed_data")
save_path <- "./scGAN_preprocessed_data"

for(i in 1:101){
  data <- readRDS(file.path("/Users/duohongrui/Desktop/preprocessed_data", data_list[i]))
  message("-------------------------")
  message(data_list[i])
  data_info <- data$data_info
  data_id <- data_info$id
  if(dynwrap::is_wrapper_with_expression(data$data)){
    data <- t(data$data[["counts"]])
  }else{
    data <- data$data
  }
  result <- simutils::scgan_data_conversion(data = data,
                                            data_id = data_id,
                                            save_to_path = save_path,
                                            verbose = FALSE)
  message("Done")
  message("-------------------------")
  Sys.sleep(2)
}
