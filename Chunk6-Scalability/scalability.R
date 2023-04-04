example_data <- readRDS("../preprocessed_data/data42_GSE65525_subset3.rds")
data <- example_data$data
group_condition <- example_data$data_info$group_condition


gradient_num <- data.frame("cell" = c(100, 200, 500, 800, 1000, 2000, 3000, 5000, 8000, 10000, rep(1000, 10)),
                           "gene" = c(rep(1000, 10), 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 10000))

## First class of methods which users can custom cell and gene number

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
  
  for(n in 1:3){
    ### estimation
    est <- simpipe::estimate_parameters(ref_data = sub_data,
                                        method = c("Simple",
                                                   "Kersplat",
                                                   "Splat",
                                                   "SplatPop",
                                                   "Lun"),
                                        seed = 111,
                                        verbose = TRUE,
                                        use_docker = FALSE)
    time <- lapply(est, function(x){x[["estimate_detection"]][1,2]})
    method_name <- stringr::str_split(names(time), "_", simplify = TRUE)[, 2]
    est_time <- as.numeric(time)
    est_memory <- as.numeric(lapply(est, function(x){x[["estimate_detection"]][1,4]}))
    
    ### simulation
    sim <- simpipe::simulate_datasets(parameters = est,
                                      seed = 111,
                                      return_format = "list",
                                      verbose = TRUE)
    sim_time <- as.numeric(lapply(sim, function(x){x[["simulate_detection"]][1,2]}))
    sim_memory <- as.numeric(lapply(sim, function(x){x[["simulate_detection"]][1,4]}))
    
    scala <- tibble::tibble("method" = method_name,
                            "cell_num" = cell_num,
                            "gene_num" = gene_num,
                            "repeat_time" = n,
                            "estimation_time" = est_time,
                            "estimation_memory" = est_memory,
                            "simulation_time" = sim_time,
                            "simulation_memory" = sim_memory)
    Sys.sleep(1)
    saveRDS(scala, file = paste0("../scalability/",
                                 "class01_",
                                 cell_num,
                                 "_",
                                 gene_num,
                                 "_", n,
                                 ".rds"))
  }
}


## Second class of different SCRIP modes

modes <- c("GP-trendedBCV", "GP-commonBCV", "BGP-commonBCV", "BP", "BGP-trendedBCV")

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
  
  
  scala_result <- purrr::map(
    .x = modes,
    .f = function(mode){
      ### estimation
      est <- simpipe::estimate_parameters(ref_data = sub_data,
                                          method = c("SCRIP"),
                                          seed = 111,
                                          verbose = TRUE,
                                          use_docker = FALSE)
      time <- lapply(est, function(x){x[["estimate_detection"]][1,2]})
      method_name <- stringr::str_split(names(time), "_", simplify = TRUE)[, 2]
      est_time <- as.numeric(time)
      est_memory <- as.numeric(lapply(est, function(x){x[["estimate_detection"]][1,4]}))
      
      ### simulation
      sim <- simpipe::simulate_datasets(ref_data = sub_data,
                                        parameters = est,
                                        other_prior = list(mode = mode),
                                        seed = 111,
                                        return_format = "list",
                                        verbose = TRUE)
      sim_time <- as.numeric(lapply(sim, function(x){x[["simulate_detection"]][1,2]}))
      sim_memory <- as.numeric(lapply(sim, function(x){x[["simulate_detection"]][1,4]}))
      
      tibble::tibble("method" = paste0(method_name, "-", mode),
                     "cell_num" = cell_num,
                     "gene_num" = gene_num,
                     "repeat_time" = n,
                     "estimation_time" = est_time,
                     "estimation_memory" = est_memory,
                     "simulation_time" = sim_time,
                     "simulation_memory" = sim_memory)
    })
  scala <- purrr::map_dfr(scala_result, .f = rbind)
  Sys.sleep(1)
  saveRDS(scala, file = paste0("../scalability/",
                               "class02_",
                               cell_num,
                               "_",
                               gene_num,
                               "_", n,
                               ".rds"))
}


## Third class of different SCRIP modes

third_class <- c("scDesign", "SPsimSeq")

for(i in 1:20){
  cell_num <- gradient_num[i, 1]
  print(cell_num)
  gene_num <- gradient_num[i, 2]
  print(gene_num)
  
  set.seed(i)
  sub_data <- scater::mockSCE(ncells = cell_num,
                              ngenes = gene_num)
  sub_data <- SingleCellExperiment::counts(sub_data)
  
  scala_result <- purrr::map(
    .x = third_class,
    .f = function(method){
      
      ### simulation
      sim <- simpipe::simulate_datasets(ref_data = sub_data,
                                        method = method,
                                        seed = 111,
                                        return_format = "list",
                                        verbose = TRUE)
      sim_time <- as.numeric(lapply(sim, function(x){x[["simulate_detection"]][1,2]}))
      sim_memory <- as.numeric(lapply(sim, function(x){x[["simulate_detection"]][1,4]}))
      
      tibble::tibble("method" = method,
                     "cell_num" = cell_num,
                     "gene_num" = gene_num,
                     "repeat_time" = n,
                     "estimation_time" = NA,
                     "estimation_memory" = NA,
                     "simulation_time" = sim_time,
                     "simulation_memory" = sim_memory)
    })
  scala <- purrr::map_dfr(scala_result, .f = rbind)
  Sys.sleep(1)
  saveRDS(scala, file = paste0("../scalability/",
                               "class03_",
                               cell_num,
                               "_",
                               gene_num,
                               "_", n,
                               ".rds"))
}





