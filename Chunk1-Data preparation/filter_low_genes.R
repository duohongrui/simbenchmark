data_list <- list.files("/Users/duohongrui/Desktop/preprocessed_data/")

for(i in 1:length(data_list)){
  ### read data
  message(data_list[i])
  ref_data <- readRDS(paste0("/Users/duohongrui/Desktop/preprocessed_data/", data_list[i]))
  ### count data
  if(dynwrap::is_wrapper_with_expression(ref_data$data)){
    next
  }else{
    counts <- as.matrix(ref_data$data)
  }
  cat(paste0("Raw gene number: ", nrow(counts), "\n"))
  
  ### gene filter
  index <- rowSums(counts) == 0
  counts <- counts[!index, ]
  cat(paste0("Retained gene number: ", nrow(counts), "\n"))
  
  ### save
  ref_data$data <- counts
  ref_data$data_info$gene_num <- nrow(counts)
  if(!dir.exists("/Users/duohongrui/Desktop/filter_preprocessed_data/")){
    dir.create("/Users/duohongrui/Desktop/filter_preprocessed_data/")
  }
  message("Save...")
  saveRDS(ref_data, file = paste0("/Users/duohongrui/Desktop/filter_preprocessed_data/", data_list[i]))
}

