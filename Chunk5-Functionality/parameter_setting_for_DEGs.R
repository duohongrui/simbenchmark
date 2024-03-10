data <- readRDS("../preprocessed_data/data1_GSE54006.rds")
data_info <- data$data_info
counts <- data$data
## data meta info
group <- data_info$group_condition

#### DEGs
result <- simutils::perform_DEA(data = counts,
                                group = group,
                                method = "edgeRQLFDetRate")
# saveRDS(result, "../DEA_result/data1_GSE54006.rds")
de_genes_per_group <- lapply(result, function(df){
  rownames(df)[df$PValue < 0.05]
})
de_genes <- unique(BiocGenerics::Reduce(x = de_genes_per_group, f = union))
prob.group <- as.numeric(table(group))/length(group)
de.prob <- length(de_genes)/nrow(counts)


estimation_result <- simpipe::estimate_parameters(
  ref_data = counts,
  method = "Splat",
  other_prior = NULL,
  seed = 666,
  verbose = TRUE)

de.facLoc <- c(0.01, 0.05, 0.1, 0.2, 0.5, 1, 2)
de.facScale <- c(0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.5)
for(i in de.facLoc){
  for(j in de.facScale){
    ### other information
    other_prior <- list(prob.group = prob.group,
                        de.prob = de.prob,
                        batchCells = ncol(counts),
                        de.facLoc = i,
                        de.facScale = j,
                        nGenes = nrow(counts))
    
    simulation_result <- simpipe::simulate_datasets(
      parameters = estimation_result,
      other_prior = other_prior,
      n = 3,
      seed = c(111, 666, 888),
      return_format = "list",
      verbose = TRUE)
    saveRDS(simulation_result,
            paste0("../parameters_tuning/DEGs/simulation_result_", i, "_", j, ".rds"))
  }
}


########### Evaluation of Simulated Cell Groups
library(tibble)
data_list <- list.files("../parameters_tuning/DEGs/", pattern = "^simulation")
for(i in data_list[34:49]){
  data <- readRDS(file.path("../parameters_tuning/DEGs", i))
  message(i)
  
  evl_group <- purrr::map_dfr(1:length(data), .f = function(x){
    simulated <- data[[x]]
    message("Calculating distance matrix...")
    dist <- parallelDist::parDist(t(as.matrix(simulated$simulate_result$count_data)), threads = 6)
    
    message("1-Calculating CDI...")
    error <- try(
      CDI <- simutils::calculate_CDI(as.matrix(simulated$simulate_result$count_data),
                                     cluster_info = simulated$simulate_result$col_meta$group),
      silent = TRUE
    )
    if("try-error" %in% class(error)){
      CDI <- NA
    }else{
      CDI <- min(CDI[1, 1], CDI[1, 2])
    }
    
    message("2-Calculating ROUGE...")
    error <- try(
      ROUGE <- simutils::calculate_ROGUE(as.matrix(simulated$simulate_result$count_data),
                                         cluster_info = simulated$simulate_result$col_meta$group),
      silent = TRUE
    )
    if("try-error" %in% class(error)){
      ROUGE <- NA
    }
    
    message("3-Calculating silhouette...")
    error <- try(
      silhouette <- simutils::calculate_silhouette(dist, cluster_info = simulated$simulate_result$col_meta$group),
      silent = TRUE
    )
    if("try-error" %in% class(error)){
      silhouette <- NA
    }
    
    message("4-Calculating dunn...")
    error <- try(
      dunn <- simutils::calculate_dunn(dist, cluster_info = simulated$simulate_result$col_meta$group),
      silent = TRUE
    )
    if("try-error" %in% class(error)){
      dunn <- NA
    }
    
    message("5-Calculating connectivity...")
    error <- try(
      connectivity <- simutils::calculate_connectivity(dist, cluster_info = simulated$simulate_result$col_meta$group),
      silent = TRUE
    )
    if("try-error" %in% class(error)){
      connectivity <- NA
    }
    
    message("6-Calculating DB index...")
    error <- try(
      DB_index <- simutils::calculate_DB_index(as.matrix(simulated$simulate_result$count_data),
                                               cluster_info = simulated$simulate_result$col_meta$group),
      silent = TRUE
    )
    if("try-error" %in% class(error)){
      DB_index <- NA
    }
    tibble::tibble(
      "meanlog" = stringr::str_split(i, pattern = "_", simplify = TRUE)[3],
      "sdlog" = gsub(".rds","",stringr::str_split(i, pattern = "_", simplify = TRUE)[4]),
      "time" = x,
      "CDI" = CDI,
      "ROUGE" = ROUGE,
      "silhouette" = silhouette,
      "dunn" = dunn,
      "connectivity" = connectivity,
      "DB_index" = DB_index
    )
  })
  saveRDS(evl_group, paste0("../parameters_tuning/DEGs/", "group_", i))
}


