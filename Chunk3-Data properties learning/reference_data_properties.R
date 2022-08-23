data_list <- list.files("/Users/duohongrui/Desktop/preprocessed_data/")

for(i in 1:50){
  ## Read data
  data_name <- data_list[i]
  data <- readRDS(paste0("/Users/duohongrui/Desktop/preprocessed_data/", data_name))
  data_id <- data$data_info$id
  ## count
  if(dynwrap::is_wrapper_with_expression(data$data)){
    counts <- t(data$data[["counts"]])
  }else{
    counts <- data$data
  }
  ## cell properties
  cell_properties <- simutils::cell_properties(data = counts, verbose = TRUE)
  ## gene properties
  gene_properties <- simutils::gene_properties(data = counts, verbose = TRUE)
  ## result
  result <- dplyr::lst(cell_properties,
                       gene_properties)
  ## save name
  save_name <- paste0(data_id, "_properties.rds")
  ## save path
  save_path <- file.path("/Volumes/Elements/ref_data_properties", save_name)
  ## save
  saveRDS(result, file = save_path)
  message(i)
  message(glue::glue("{data_id} is done..."))
  message("----------------------------------")
}

