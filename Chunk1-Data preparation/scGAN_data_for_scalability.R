data <- readRDS("/Users/duohongrui/Desktop/preprocessed_data/data1_GSE54006.rds")
group_condition <- data[["data_info"]][["group_condition"]]
data <- data$data
save_path <- "./scGAN_preprocessed_data"

gradient_num <- data.frame("cell" = c(100, 200, 500, 800, 1000, 2000, 3000, 5000, 8000, 10000, rep(1000, 10)),
                           "gene" = c(rep(1000, 10), 500, 800, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 10000))
gradient_num <- data.frame("cell" = c(150, 300, 400, 600, 700, 900, 1500, 4000, 6000, 7000, 9000, rep(1000, 9)),
                           "gene" = c(rep(1000, 11), 1000, 1500, 3500, 4500, 5500, 6500, 7500, 8500, 9500))
cluster_info <- list()
for(i in 1:20){
  
  cell_num <- gradient_num[i, 1]
  print(cell_num)
  gene_num <- gradient_num[i, 2]
  print(gene_num)
  
  set.seed(i)
  sample_index <- sample(ncol(data), size = cell_num, replace = TRUE)
  set.seed(i)
  gene_index <- sample(nrow(data), size = gene_num, replace = TRUE)
  
  sub_data <- data[gene_index, sample_index]
  group <- group_condition[sample_index]
  rownames(sub_data) <- paste0("Gene", 1:nrow(sub_data))
  colnames(sub_data) <- paste0("Cell", 1:ncol(sub_data))
  
  message("-------------------------")
  group <- group - 1
  data_id <- paste0("scalability_", cell_num, "_", gene_num)
  result <- simutils::scgan_data_conversion(data = sub_data,
                                            data_i = data_id,
                                            group = group,
                                            save_to_path = save_path,
                                            verbose = FALSE)
  cluster_info[[data_id]] <- group + 1
  message("Done")
  message("-------------------------")
  Sys.sleep(2)
}

saveRDS(cluster_info, file = "/Users/duohongrui/Desktop/scGAN_preprocessed_data/cluster_info.rds")
