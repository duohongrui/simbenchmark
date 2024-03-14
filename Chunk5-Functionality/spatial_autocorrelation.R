library(dplyr)
data_list <- list.files("../simulation_data/", pattern = "spatial")

for(file_name in data_list){
  data_name <- stringr::str_split(file_name, pattern = "_", simplify = TRUE)[, 2]
  method <- stringr::str_split(file_name, pattern = "_", simplify = TRUE)[, 1]
  
  simulated_data <- readRDS(paste0("../simulation_data/", file_name))
  
  if(simulated_data$sim_data_info$group >= 2 & "de_gene" %in% colnames(simulated_data$sim_data$row_meta) |
     simulated_data$sim_data_info$group >= 2 & "DEstatus" %in% colnames(simulated_data$sim_data$row_meta)){
    message(file_name)
    data <- readRDS(paste0("../preprocessed_data/", list.files("../preprocessed_data/", pattern = data_name)))
    if(dynwrap::is_wrapper_with_expression(data$data)){
      ref_data <- t(as.matrix(data$data$counts))
    }else{
      ref_data <- data$data
    }
    locations <- data$data_info$spatial_coordinate
    count_data <- simulated_data$sim_data$count_data %>% as.matrix()
    if(ncol(count_data) != ncol(ref_data)){
      next
    }
    
    cat("Performing PCA...\n")
    error3 <- try(
      pca_ref <- prcomp(t(ref_data), rank. = 50),
      silent = TRUE
    )
    if(is(error3, "try-error")){
      next
    }
    error1 <- try(
      pca_sim <- prcomp(t(count_data), rank. = 50),
      silent = TRUE
    )
    if(is(error1, "try-error")){
      next
    }
    cat("Performing Correlation...\n")
    cor_result <- WGCNA::cor(t(rbind(pca_ref$x, pca_sim$x)), method = 'spearman')
    index <- dim(cor_result)[1]/2
    cor_result <- cor_result[(index+1):(index*2), 1:index]
    error2 <- try(
      Hungarian_result <- RcppHungarian::HungarianSolver(-cor_result),
      silent = TRUE
    )
    if(is(error2, "try-error")){
      next
    }
    match_value <- purrr::map(1:ncol(cor_result), .f = function(x){
      cor_result[Hungarian_result[["pairs"]][x, 1], Hungarian_result[["pairs"]][x, 2]]
    }) %>% unlist()
    cell_pair <- data.frame("reference" = colnames(cor_result)[Hungarian_result[["pairs"]][, 2]],
                            "simulation" = rownames(cor_result)[Hungarian_result[["pairs"]][, 1]],
                            "match_value" = match_value)
    print(utils::head(cell_pair))
    
    if(stringr::str_starts(file_name, pattern = "^scDD_")){
      all_svgs_index <- stringr::str_starts(simulated_data[["sim_data"]][["row_meta"]][["DEstatus"]], pattern = "^D")
      all_svgs <- rownames(count_data)[all_svgs_index]
    }else{
      all_svgs <- rownames(count_data)[which(simulated_data[["sim_data"]][["row_meta"]][["de_gene"]] == "yes")]
    }
    
    #### All SVGs
    index <- match(cell_pair$reference, colnames(ref_data))
    count_data <- count_data[all_svgs, ]
    moran_result <- tibble()
    print(length(all_svgs))
    for(i in 1:length(all_svgs)){
      candicate_data <- locations[index, ] %>%
        mutate(
          "gene" = count_data[i, ]
        )
      error <- try(
        result <- moranfast::moranfast(candicate_data$gene, candicate_data$x, candicate_data$y),
        silent = TRUE
      )
      if(is(error, "try-error")){
        next
      }
      moran_result <- moran_result %>% 
        rbind(tibble(
          "method" = method,
          "data" = data_name,
          "gene" = rownames(count_data)[i],
          "Moran" = result$observed
        ))
    }
    saveRDS(moran_result, file = paste0("../spatial_autocorrelation/", file_name))
  }else{
    next
  }
}
