sim_data_list <- list.files("../../splatter_package_methods/simulation_data/")
ref_data_list <- list.files("/Users/duohongrui/Desktop/preprocessed_data/")

all_result <- list()
for(ref_data_name in ref_data_list){
  ### read reference data
  print(ref_data_name)
  ref_data <- readRDS(paste0("/Users/duohongrui/Desktop/preprocessed_data/", ref_data_name))
  if(dynwrap::is_wrapper_with_expression(ref_data$data)){
    ref_data <- t(ref_data$data$counts)
  }else{
    ref_data <- ref_data$data
  }
  ref_data_cell_properties <- simutils::cell_properties(ref_data, verbose = TRUE)
  ref_data_gene_properties <- simutils::gene_properties(ref_data, verbose = TRUE)
  
  data_num <- stringr::str_split(ref_data_name, pattern = "_", simplify = TRUE)[1]
  
  for(sim_data_name in sim_data_list){
    ## read simulated data
    sim_data <- readRDS(paste0("./", sim_data_name))
    sim_data <- sim_data$sim_data$count_data
    sim_data_num <- stringr::str_split(sim_data_name, pattern = "_", simplify = TRUE)[2]
    if(data_num != sim_data_num){
      next
    }
    print(sim_data_name)
    sim_data_cell_properties <- simutils::cell_properties(sim_data, verbose = TRUE)
    sim_data_gene_properties <- simutils::gene_properties(sim_data, verbose = TRUE)

    data_properties <- simpipe::data_properties_summary(
      ref_data_cell_properties = ref_data_cell_properties,
      ref_data_gene_properties = ref_data_gene_properties,
      sim_data_cell_properties = sim_data_cell_properties,
      sim_data_gene_properties = sim_data_gene_properties
    )
    
    all_result[[sim_data_name]] <- data_properties
    
    saveRDS(all_result, file = "/Users/duohongrui/Desktop/accuracy.rds")
  }
}

#### Calculate spatial-level metrics
sim_data_list <- list.files("../../splatter_package_methods/simulation_data/")
ref_data_list <- list.files("/Users/duohongrui/Desktop/preprocessed_data/", pattern = "spatial")

all_result <- list()
for(ref_data_name in ref_data_list){
  ### read reference data
  print(ref_data_name)
  ref_data <- readRDS(paste0("/Users/duohongrui/Desktop/preprocessed_data/", ref_data_name))
  if(dynwrap::is_wrapper_with_expression(ref_data$data)){
    ref_data <- t(ref_data$data$counts)
  }else{
    ref_data <- ref_data$data
  }
  ref_data_cell_properties <- simutils::cell_properties(ref_data, verbose = TRUE)
  ref_data_gene_properties <- simutils::gene_properties(ref_data, verbose = TRUE)
  data_num <- stringr::str_split(ref_data_name, pattern = "_", simplify = TRUE)[1]
  for(sim_data_name in sim_data_list){
    ## read simulated data
    sim_data <- readRDS(paste0("./", sim_data_name))
    sim_data <- sim_data$sim_data$count_data
    sim_data_num <- stringr::str_split(sim_data_name, pattern = "_", simplify = TRUE)[2]
    if(data_num != sim_data_num){
      next
    }
    match_result <- simutils::match_cells(ref_data = ref_data, sim_data, t = TRUE)
    print(sim_data_name)
    sim_data_cell_properties <- simutils::cell_properties(sim_data, verbose = TRUE)
    sim_data_gene_properties <- simutils::gene_properties(sim_data, verbose = TRUE)
    
    spatial_data_propety <- purrr::map_dfr(names(sim_data_cell_properties), .f = function(x){
      sim_property <- cbind(spatial_coordinate,
                            data.frame("x" = sim_data_cell_properties[[x]]))
      real_property <- cbind(spatial_coordinate,
                             data.frame("x" = ref_data_cell_properties[[x]]))
      spatial-multiKS <- fasano.franceschini.test::fasano.franceschini.test(real_property, sim_property)
      spatial-KDE <- ks::kde.test(real_property, sim_property)
      tibble::tibble(
        "Method" = method,
        "Data" = data_name,
        "spatial-multiKS" = spatial-multiKS,
        "spatial-KDE" = spatial-KDE
      )
    })
    saveRDS(spatial_data_propety, file = "/Users/duohongrui/Desktop/accuracy.rds")
  }
}






