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
library(dplyr)
sim_data_list <- list.files("../simulation_data/", pattern = "spatial")
ref_data_list <- list.files("../preprocessed_data/", pattern = "spatial")

for(ref_data_name in ref_data_list){
  ### read reference data
  print(which(ref_data_list %in% ref_data_name))
  print(ref_data_name)
  ref_data <- readRDS(paste0("../preprocessed_data/", ref_data_name))
  if(dynwrap::is_wrapper_with_expression(ref_data$data)){
    ref_count_data <- t(as.matrix(ref_data$data$counts))
  }else{
    ref_count_data <- ref_data$data
  }
  ref_spatial_data <- ref_data$data_info$spatial_coordinate %>% 
    tibble::rownames_to_column("reference")
  ref_data_cell_properties <- simutils::cell_properties(ref_count_data, verbose = TRUE)
  data_num <- stringr::str_split(ref_data_name, pattern = "_", simplify = TRUE)[1]
  
  sim_data_list <- list.files("../sim_data_properties/", pattern = data_num)
  for(sim_data_name in sim_data_list){
    ## read simulated data
    cat(sim_data_name)
    sim_data_cell_properties <- readRDS(paste0("../sim_data_properties/", sim_data_name))
    sim_data <- readRDS(paste0("../simulation_data/", sim_data_name))
    sim_count_data <- sim_data$sim_data$count_data %>% as.matrix()
    method <- stringr::str_split(sim_data_name, "_", simplify = TRUE)[1]
    if(!identical(dim(sim_count_data), dim(ref_count_data))){
      next
    }
    if(method == "SRTsim" |
       method == "scDesign3"){
      match_result <- data.frame("reference" = colnames(ref_count_data),
                                 "simulation" = colnames(sim_count_data),
                                 "x" = ref_spatial_data$x,
                                 "y" = ref_spatial_data$y)
      arrange_index <- 1:ncol(ref_count_data)
    }else{
      match_result <- simutils::match_cells(ref_data = t(ref_count_data), t(sim_count_data), algorithm = "Hungarian")
      cell_match <- match_result$cell_pair
      match_result <- dplyr::left_join(cell_match, ref_spatial_data, "reference")
      match_result <- match_result[, -3]
      arrange_index <- match(match_result$reference, colnames(ref_count_data))
    }
    save_name <- paste0("../spatial_metrics_results/", method, "_", data_num, "_cell", ".rds")
    if(file.exists(save_name)){
      next
    }
    error <- try(
      spatial_data_cell_propety <- purrr::map_dfr(names(sim_data_cell_properties$cell_properties), .f = function(x){
        print(x)
        sim_cell_property <- sim_data_cell_properties$cell_properties[[x]]
        ref_cell_property <- ref_data_cell_properties[[x]][arrange_index]
        if(length(sim_cell_property) != ncol(sim_count_data)){
          tibble::tibble()
        }else{
          sim_property <- match_result %>% 
            select(3,4) %>% 
            mutate("property" = sim_cell_property)
          real_property <- match_result %>% 
            select(3,4) %>% 
            mutate("property" = ref_cell_property)
          spatial_multiKS <- fasano.franceschini.test::fasano.franceschini.test(real_property, sim_property)
          spatial_KDE <- ks::kde.test(real_property, sim_property)
          tibble::tibble(
            "Method" = method,
            "Data" = data_num,
            "data_property" = x,
            "spatial_multiKS" = mean(spatial_multiKS$estimate),
            "spatial_KDE" = spatial_KDE$zstat
          )
        }
      }), silent = TRUE
    )
    if(is(error, "try-error")){
      next
    }
    saveRDS(spatial_data_cell_propety, file = paste0("../spatial_metrics_results/", method, "_", data_num, "_cell", ".rds"))
  }
}


