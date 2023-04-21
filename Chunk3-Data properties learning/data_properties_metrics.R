library(simpipe)
sim_data_list <- list.files("../simulation_data/", pattern = "^SPsimSeq")

for(i in sim_data_list){
  print(i)
  message("Read simulated data...")
  sim_result <- readRDS(file.path("../simulation_data", i))
  sim_data <- as.data.frame(sim_result$sim_data$count_data)
  
  if(any(colSums(sim_data)==0)){
    next
  }
  
  message("Calculate cell properties...")
  sim_data_cell_properties <- simutils::cell_properties(sim_data, verbose = TRUE)
  message("Calculate gene properties...")
  sim_data_gene_properties <- simutils::gene_properties(sim_data, verbose = TRUE)
  
  ### save data property
  sim_data_properties <- dplyr::lst(cell_properties = sim_data_cell_properties,
                                    gene_properties = sim_data_gene_properties)
  
  saveRDS(sim_data_properties, file.path("../sim_data_properties",
                                         paste0(sim_result[["sim_data_info"]][["sim_data_id"]], ".rds")))
  
  ### Read real data property
  message("Read real data properties...")
  data_name <- stringr::str_extract(i, pattern = "data[0-9]+[_]")
  ref_data_properties <- readRDS(file.path("../ref_data_properties",
                                           list.files("../ref_data_properties", pattern = data_name)))
  ref_data_cell_properties <- ref_data_properties$cell_properties
  ref_data_gene_properties <- ref_data_properties$gene_properties
  
  ### Calculate similarity
  
  cell_num <- length(ref_data_cell_properties[[1]]) != length(sim_data_cell_properties[[1]])
  gene_num <- length(ref_data_gene_properties[[1]]) != length(sim_data_gene_properties[[1]])
  
  if(cell_num & gene_num){
    message("The numbers of cell(&gene) are not identical...")
    next
  }
  
  message("Calculating metrics...")
  result <- simpipe::data_properties_summary(ref_data_cell_properties = ref_data_cell_properties,
                                             sim_data_cell_properties = sim_data_cell_properties,
                                             ref_data_gene_properties = ref_data_gene_properties,
                                             sim_data_gene_properties = sim_data_gene_properties,
                                             ncore = 3)
  ### save result
  message("Save results...")
  saveRDS(result, file.path("../data_properties_result", i))
}
