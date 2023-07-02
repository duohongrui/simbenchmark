example_data <- readRDS("../preprocessed_data/data42_GSE65525_subset3.rds")
data <- example_data$data
group_condition <- example_data$data_info$group_condition


gradient_num <- data.frame("cell" = c(100, 200, 500, 800, 1000, 2000, 3000, 5000, 8000, 10000, rep(1000, 10)),
                           "gene" = c(rep(1000, 10), 500, 800, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 10000))
gradient_num <- data.frame("cell" = c(150, 300, 400, 600, 700, 900, 1500, 4000, 6000, 7000, 9000, rep(1000, 9)),
                           "gene" = c(rep(1000, 11), 1000, 1500, 3500, 4500, 5500, 6500, 7500, 8500, 9500))
## First class of methods which users can custom cell and gene number
method <- c("Simple",
            "Kersplat",
            "Splat",
            "SplatPop",
            "Lun")
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
                                        method = method,
                                        seed = 111,
                                        verbose = TRUE)
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
  
  for(n in 1:3){
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
}


## Third class of scDesign and SPsimSeq

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
  
  for(n in 1:3){
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
}

## Forth class of methods which need group information

example_data <- readRDS("../preprocessed_data/data42_GSE65525_subset3.rds")
data <- example_data$data
group_condition <- example_data$data_info$group_condition

fouth_class <- c("SPARSim", "powsimR", "POWSC", "scDesign2", "muscat")

for(i in 1:20){
  
  cell_num <- gradient_num[i, 1]
  print(cell_num)
  gene_num <- gradient_num[i, 2]
  print(gene_num)
  
  set.seed(i*10)
  sample_index <- sample(ncol(data), size = cell_num, replace = TRUE)
  set.seed(i*10)
  gene_index <- sample(nrow(data), size = gene_num, replace = TRUE)
  
  sub_data <- data[gene_index, sample_index]
  group <- group_condition[sample_index]
  rownames(sub_data) <- paste0("Gene", 1:nrow(sub_data))
  colnames(sub_data) <- paste0("Cell", 1:ncol(sub_data))
  
  for(n in 1:3){
    scala_result <- purrr::map(
      .x = fouth_class,
      .f = function(method){
        
        ### estimation
        est <- simpipe::estimate_parameters(ref_data = sub_data,
                                            method = method,
                                            other_prior = list(group.condition = group),
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
        
        tibble::tibble("method" = method,
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
                                 "class04_",
                                 cell_num,
                                 "_",
                                 gene_num,
                                 "_", n,
                                 ".rds"))
  }
}




## Fifth class of methods
library(ESCO)
fifth_class <- c("ESCO", "zinbwave", "hierarchicell", "dropsim")

for(i in 1:20){
  
  cell_num <- gradient_num[i, 1]
  print(cell_num)
  gene_num <- gradient_num[i, 2]
  print(gene_num)
  
  set.seed(i)
  sub_data <- SingleCellExperiment::counts(scater::mockSCE(ncells = cell_num, ngenes = gene_num, nspikes = 0))
  
  rownames(sub_data) <- paste0("Gene", 1:nrow(sub_data))
  colnames(sub_data) <- paste0("Cell", 1:ncol(sub_data))
  
  for(n in 1:3){
    scala_result <- purrr::map(
      .x = fifth_class,
      .f = function(method){
        
        ### estimation
        est <- simpipe::estimate_parameters(ref_data = sub_data,
                                            method = method,
                                            seed = 111,
                                            verbose = TRUE,
                                            use_docker = FALSE)
        time <- lapply(est, function(x){x[["estimate_detection"]][1,2]})
        method_name <- stringr::str_split(names(time), "_", simplify = TRUE)[, 2]
        est_time <- as.numeric(time)
        est_memory <- as.numeric(lapply(est, function(x){x[["estimate_detection"]][1,4]}))
        
        ### simulation
        if(method == "ESCO"){
          sim <- simpipe::simulate_datasets(parameters = NULL,
                                            method = "ESCO",
                                            other_prior = list(nCells = cell_num,
                                                               nGenes = gene_num),
                                            seed = 111,
                                            return_format = "list",
                                            verbose = TRUE)
        }else{
          sim <- simpipe::simulate_datasets(parameters = est,
                                            other_prior = list(nCells = cell_num,
                                                               nGenes = gene_num),
                                            seed = 111,
                                            return_format = "list",
                                            verbose = TRUE)
        }
        
        sim_time <- as.numeric(lapply(sim, function(x){x[["simulate_detection"]][1,2]}))
        sim_memory <- as.numeric(lapply(sim, function(x){x[["simulate_detection"]][1,4]}))
        
        tibble::tibble("method" = method,
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
                                 "class05_",
                                 cell_num,
                                 "_",
                                 gene_num,
                                 "_", n,
                                 ".rds"))
  }
}


## Sixth class of methods
sixth_class <- c("SparseDC", "zinbwaveZinger")
library(simmethods)
for(i in 1:20){
  
  cell_num <- gradient_num[i, 1]
  print(cell_num)
  gene_num <- gradient_num[i, 2]
  print(gene_num)
  
  set.seed(i)
  sub_data <- SingleCellExperiment::counts(scater::mockSCE(ncells = cell_num, ngenes = gene_num, nspikes = 0))
  set.seed(i)
  group <- sample(1:2, ncol(sub_data), replace = TRUE)
  rownames(sub_data) <- paste0("Gene", 1:nrow(sub_data))
  colnames(sub_data) <- paste0("Cell", 1:ncol(sub_data))
  
  for(n in 1:3){
    scala_result <- purrr::map(
      .x = sixth_class,
      .f = function(method){
        
        ### estimation
        est <- simpipe::estimate_parameters(ref_data = sub_data,
                                            method = method,
                                            other_prior = list(group.condition = group,
                                                               nclusters = 2),
                                            seed = 111,
                                            verbose = TRUE,
                                            use_docker = FALSE)
        time <- lapply(est, function(x){x[["estimate_detection"]][1,2]})
        method_name <- stringr::str_split(names(time), "_", simplify = TRUE)[, 2]
        est_time <- as.numeric(time)
        est_memory <- as.numeric(lapply(est, function(x){x[["estimate_detection"]][1,4]}))
        
        ### simulation
        sim <- simpipe::simulate_datasets(parameters = est,
                                          ref_data = sub_data,
                                          other_prior = list(group.condition = group),
                                          seed = 111,
                                          return_format = "list",
                                          verbose = TRUE)
        
        sim_time <- as.numeric(lapply(sim, function(x){x[["simulate_detection"]][1,2]}))
        sim_memory <- as.numeric(lapply(sim, function(x){x[["simulate_detection"]][1,4]}))
        
        tibble::tibble("method" = method,
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
                                 "class06_",
                                 cell_num,
                                 "_",
                                 gene_num,
                                 "_", n,
                                 ".rds"))
  }
}

## Seventh class of methods (BEARscc)
example_data <- readRDS("../preprocessed_data/data23_GSE62270.rds")
data <- example_data$data
ERCC_count <- data[grep(rownames(data), pattern = "^ERCC"), ]
data <- data[-grep(rownames(data), pattern = "^ERCC"), ]
library(SingleCellExperiment)

ninth_class <- c("BEARscc")
for(i in 1:20){
  
  cell_num <- gradient_num[i, 1]
  print(cell_num)
  gene_num <- gradient_num[i, 2]
  print(gene_num)
  
  set.seed(i)
  sample_index <- sample(ncol(data), size = cell_num, replace = TRUE)
  set.seed(i)
  gene_index <- sample(nrow(data), size = gene_num-92, replace = TRUE)
  
  sub_data <- rbind(data[gene_index, sample_index],
                    ERCC_count[, sample_index])
  
  for(n in 1:3){
    scala_result <- purrr::map(
      .x = ninth_class,
      .f = function(method){
        
        ### estimation
        est <- simpipe::estimate_parameters(ref_data = sub_data,
                                            method = method,
                                            other_prior = list(dilution.factor = 50000,
                                                               volume = 0.03),
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
                                          verbose = FALSE)
        
        sim_time <- as.numeric(lapply(sim, function(x){x[["simulate_detection"]][1,2]}))
        sim_memory <- as.numeric(lapply(sim, function(x){x[["simulate_detection"]][1,4]}))
        
        tibble::tibble("method" = method,
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
                                 "class07_",
                                 cell_num,
                                 "_",
                                 gene_num,
                                 "_", n,
                                 ".rds"))
  }
}


## Eighth class of methods (scDD)
library(scDD)
data(scDatExSim)
data(scDatEx)
data <- SingleCellExperiment::normcounts(scDatExSim)
group_condition <- SingleCellExperiment::colData(scDatExSim)$condition

eighth_class <- c("scDD")
for(i in 1:20){
  
  cell_num <- gradient_num[i, 1]
  print(cell_num)
  gene_num <- gradient_num[i, 2]
  print(gene_num)
  
  set.seed(i+666)
  sample_index <- sample(ncol(data), size = cell_num, replace = TRUE)
  set.seed(i+666)
  gene_index <- sample(nrow(data), size = gene_num, replace = TRUE)

  sub_data <- data[gene_index, sample_index]
  group <- group_condition[sample_index]
  rownames(sub_data) <- paste0("Gene", 1:nrow(sub_data))
  colnames(sub_data) <- paste0("Cell", 1:ncol(sub_data))
  
  for(n in 1:3){
    scala_result <- purrr::map(
      .x = eighth_class,
      .f = function(method){
        
        ### estimation
        est <- simpipe::estimate_parameters(ref_data = sub_data,
                                            method = method,
                                            other_prior = list(group.condition = group),
                                            seed = 111,
                                            verbose = TRUE,
                                            use_docker = FALSE)
        time <- lapply(est, function(x){x[["estimate_detection"]][1,2]})
        method_name <- stringr::str_split(names(time), "_", simplify = TRUE)[, 2]
        est_time <- as.numeric(time)
        est_memory <- as.numeric(lapply(est, function(x){x[["estimate_detection"]][1,4]}))
        
        ### simulation
        parameters <- splatter::newSCDDParams()
        parameters <- splatter::setParams(parameters,
                                          list(nDE = gene_num/5,
                                               nDP = gene_num/5,
                                               nDM = gene_num/5,
                                               nDB = gene_num/5,
                                               nEE = gene_num/5,
                                               nEP = 0,
                                               SCdat = scDatEx,
                                               seed = i))
        simulate_detection <- peakRAM::peakRAM(
          simulate_result <- splatter::scDDSimulate(parameters,
                                                    verbose = TRUE)
        )
        sim_time <- simulate_detection[1,2]
        sim_memory <- simulate_detection[1,4]
        tibble::tibble("method" = method,
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
                                 "class08_",
                                 cell_num,
                                 "_",
                                 gene_num,
                                 "_", n,
                                 ".rds"))
  }
}
## BPPARAM = BiocParallel::MulticoreParam(workers = 2)

## Ninth class of methods (BACiCS)
seventh_class <- c("BASiCS")
for(i in 1:20){
  
  cell_num <- gradient_num[i, 1]
  print(cell_num)
  gene_num <- gradient_num[i, 2]
  print(gene_num)
  
  set.seed(i*10)
  sub_data <- cbind(matrix(rpois(cell_num * gene_num/2, 2), nrow = gene_num, ncol = cell_num/2),
                    matrix(rpois(cell_num * gene_num/2, 6), nrow = gene_num, ncol = cell_num/2))
  batch <- c(rep(1, cell_num/2), rep(2, cell_num/2))
  rownames(sub_data) <- paste0("Gene", 1:nrow(sub_data))
  colnames(sub_data) <- paste0("Cell", 1:ncol(sub_data))
  
  for(n in 1:3){
    scala_result <- purrr::map(
      .x = seventh_class,
      .f = function(method){
        
        ### estimation
        est <- simpipe::estimate_parameters(ref_data = sub_data,
                                            method = method,
                                            other_prior = list(batch.condition = batch,
                                                               n = 6000),
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
        
        tibble::tibble("method" = method,
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
                                 "class09_",
                                 cell_num,
                                 "_",
                                 gene_num,
                                 "_", n,
                                 ".rds"))
  }
}


## Tenth class of methods (Lun2)
tenth_class <- c("Lun2")
for(i in 1:20){
  
  cell_num <- gradient_num[i, 1]
  print(cell_num)
  gene_num <- gradient_num[i, 2]
  print(gene_num)
  
  set.seed(i*10)
  sub_data <- cbind(matrix(rpois(cell_num * gene_num/2, 2), nrow = gene_num, ncol = cell_num/2),
                    matrix(rpois(cell_num * gene_num/2, 6), nrow = gene_num, ncol = cell_num/2))
  group <- c(rep(1, cell_num/2), rep(2, cell_num/2))
  rownames(sub_data) <- paste0("Gene", 1:nrow(sub_data))
  colnames(sub_data) <- paste0("Cell", 1:ncol(sub_data))
  
  for(n in 1:3){
    scala_result <- purrr::map(
      .x = tenth_class,
      .f = function(method){
        
        ### estimation
        est <- simpipe::estimate_parameters(ref_data = sub_data,
                                            method = method,
                                            other_prior = list(group.condition = group),
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
        
        tibble::tibble("method" = method,
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
                                 "class10_",
                                 cell_num,
                                 "_",
                                 gene_num,
                                 "_", n,
                                 ".rds"))
  }
}


## Eleventh class of methods (zingeR)
library(simmethods)
data <- simmethods::data
group_condition <- as.numeric(simmethods::group_condition)

eleventh_class <- c("zingeR")
for(i in 10:20){
  
  cell_num <- gradient_num[i, 1]
  print(cell_num)
  gene_num <- gradient_num[i, 2]
  print(gene_num)
  
  set.seed(i)
  sample_index <- sample(ncol(data), size = cell_num, replace = TRUE)
  set.seed(i)
  gene_index <- sample(nrow(data), size = gene_num, replace = TRUE)
  
  sub_data <- data[gene_index, sample_index]
  while(any(colSums(sub_data) == 0)){
    message("Resampling...")
    zero_cell_index <- which(colSums(sub_data) == 0)
    new_index <- sample(ncol(data), size = length(zero_cell_index), replace = TRUE)
    for(j in 1:length(new_index)){
      sub_data[, zero_cell_index[j]] <- data[gene_index, new_index[j]]
    }
  }
  group <- group_condition[sample_index]
  rownames(sub_data) <- paste0("Gene", 1:nrow(sub_data))
  colnames(sub_data) <- paste0("Cell", 1:ncol(sub_data))
  
  for(n in 1:3){
    scala_result <- purrr::map(
      .x = eleventh_class,
      .f = function(method){
        
        ### estimation
        est <- simpipe::estimate_parameters(ref_data = sub_data,
                                            method = method,
                                            other_prior = list(group.condition = group),
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
                                          other_prior = list(group.condition = group),
                                          seed = 111,
                                          return_format = "list",
                                          verbose = TRUE)
        
        sim_time <- as.numeric(lapply(sim, function(x){x[["simulate_detection"]][1,2]}))
        sim_memory <- as.numeric(lapply(sim, function(x){x[["simulate_detection"]][1,4]}))
        
        tibble::tibble("method" = method,
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
                                 "class11_",
                                 cell_num,
                                 "_",
                                 gene_num,
                                 "_", n,
                                 ".rds"))
  }
}



## Twelfth class of methods
library(simmethods)
data <- readRDS("../preprocessed_data/data85_cellbench-SC4_luyitian.rds")
group_condition <- data$data_info$cluster_info
data <- t(data$data$counts)

twelfth_class <- c("SCRIP-paths", "PROSSTT", "TedSim", "dyntoy", "SymSim", "VeloSim", "MFA", "phenopath")
for(i in 1:20){
  
  cell_num <- gradient_num[i, 1]
  print(cell_num)
  gene_num <- gradient_num[i, 2]
  print(gene_num)
  
  set.seed(i * 10)
  sample_index <- sample(ncol(data), size = cell_num, replace = TRUE)
  set.seed(i * 10)
  gene_index <- sample(nrow(data), size = gene_num, replace = TRUE)
  
  sub_data <- data[gene_index, sample_index]
  group <- as.numeric(factor(group_condition[sample_index]))
  rownames(sub_data) <- paste0("Gene", 1:nrow(sub_data))
  colnames(sub_data) <- paste0("Cell", 1:ncol(sub_data))
  
  for(n in 1:3){
    scala_result <- purrr::map(
      .x = twelfth_class,
      .f = function(method){
        
        ### SCRIP method
        if(method == "SCRIP-paths"){
          sub_method <- "SCRIP"
        }else{
          sub_method <- method
        }
        
        ### TedSim
        if(sub_method == "TedSim"){
          if(cell_num == 800 | cell_num == 3000 | cell_num == 5000 | cell_num == 10000){
            tibble::tibble("method" = "TedSim",
                           "cell_num" = NA,
                           "gene_num" = NA,
                           "repeat_time" = NA,
                           "estimation_time" = NA,
                           "estimation_memory" = NA,
                           "simulation_time" = NA,
                           "simulation_memory" = NA)
          }else{
            ### estimation
            est <- simpipe::estimate_parameters(ref_data = sub_data,
                                                method = sub_method,
                                                other_prior = list(group.condition = group),
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
                                              other_prior = list(paths = TRUE),
                                              seed = 111,
                                              return_format = "list",
                                              verbose = TRUE)
            
            sim_time <- as.numeric(lapply(sim, function(x){x[["simulate_detection"]][1,2]}))
            sim_memory <- as.numeric(lapply(sim, function(x){x[["simulate_detection"]][1,4]}))
            
            tibble::tibble("method" = method,
                           "cell_num" = cell_num,
                           "gene_num" = gene_num,
                           "repeat_time" = n,
                           "estimation_time" = est_time,
                           "estimation_memory" = est_memory,
                           "simulation_time" = sim_time,
                           "simulation_memory" = sim_memory)
          }
        }else{
          ### estimation
          est <- simpipe::estimate_parameters(ref_data = sub_data,
                                              method = sub_method,
                                              other_prior = list(group.condition = group),
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
                                            other_prior = list(paths = TRUE),
                                            seed = 111,
                                            return_format = "list",
                                            verbose = TRUE)
          
          sim_time <- as.numeric(lapply(sim, function(x){x[["simulate_detection"]][1,2]}))
          sim_memory <- as.numeric(lapply(sim, function(x){x[["simulate_detection"]][1,4]}))
          
          tibble::tibble("method" = method,
                         "cell_num" = cell_num,
                         "gene_num" = gene_num,
                         "repeat_time" = n,
                         "estimation_time" = est_time,
                         "estimation_memory" = est_memory,
                         "simulation_time" = sim_time,
                         "simulation_memory" = sim_memory)
        }
      })
    scala <- purrr::map_dfr(scala_result, .f = rbind)
    Sys.sleep(1)
    saveRDS(scala, file = paste0("../scalability/",
                                 "class12_",
                                 cell_num,
                                 "_",
                                 gene_num,
                                 "_", n,
                                 ".rds"))
  }
}



## Thirteenth class of methods 
data <- readRDS("../preprocessed_data/data85_cellbench-SC4_luyitian.rds")
group_condition <- data$data_info$cluster_info
data <- t(data$data$counts)

thirteenth_class <- c("Splat-paths", "SplatPop-paths", "ESCO-traj", "ESCO-tree")
for(i in 1:20){
  
  cell_num <- gradient_num[i, 1]
  print(cell_num)
  gene_num <- gradient_num[i, 2]
  print(gene_num)
  
  set.seed(i * 10)
  sample_index <- sample(ncol(data), size = cell_num, replace = TRUE)
  set.seed(i * 10)
  gene_index <- sample(nrow(data), size = gene_num, replace = TRUE)
  
  sub_data <- data[gene_index, sample_index]
  group <- as.numeric(factor(group_condition[sample_index]))
  rownames(sub_data) <- paste0("Gene", 1:nrow(sub_data))
  colnames(sub_data) <- paste0("Cell", 1:ncol(sub_data))
  
  for(n in 1:3){
    scala_result <- purrr::map(
      .x = thirteenth_class,
      .f = function(method){
        
        ### SCRIP method
        if(method == "Splat-paths"){
          sub_method <- "Splat"
          est_prior <- NULL
          sim_prior <- list(paths = TRUE)
        }
        if(method == "SplatPop-paths"){
          sub_method <- "SplatPop"
          est_prior <- NULL
          sim_prior <- list(paths = TRUE)
        }
        if(method == "ESCO-traj"){
          sub_method <- "ESCO"
          est_prior <- list(tree = TRUE,
                            group.condition = group)
          sim_prior <- list(type = "traj")
        }
        if(method == "ESCO-tree"){
          sub_method <- "ESCO"
          est_prior <- list(tree = TRUE,
                            group.condition = group)
          sim_prior <- list(type = "tree")
        }
        
        ### estimation
        est <- simpipe::estimate_parameters(ref_data = sub_data,
                                            method = sub_method,
                                            other_prior = est_prior,
                                            seed = 111,
                                            verbose = TRUE,
                                            use_docker = FALSE)
        time <- lapply(est, function(x){x[["estimate_detection"]][1,2]})
        method_name <- stringr::str_split(names(time), "_", simplify = TRUE)[, 2]
        est_time <- as.numeric(time)
        est_memory <- as.numeric(lapply(est, function(x){x[["estimate_detection"]][1,4]}))
        
        ### simulation
        sim <- simpipe::simulate_datasets(parameters = est,
                                          other_prior = sim_prior,
                                          seed = 111,
                                          return_format = "list",
                                          verbose = TRUE)
        
        sim_time <- as.numeric(lapply(sim, function(x){x[["simulate_detection"]][1,2]}))
        sim_memory <- as.numeric(lapply(sim, function(x){x[["simulate_detection"]][1,4]}))
        
        tibble::tibble("method" = method,
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
                                 "class13_",
                                 cell_num,
                                 "_",
                                 gene_num,
                                 "_", n,
                                 ".rds"))
  }
}



## Fourteenth class of methods 
data <- readRDS("../preprocessed_data/data85_cellbench-SC4_luyitian.rds")
group_condition <- data$data_info$cluster_info
data <- t(data$data$counts)

fourteenth_class <- c("dyngen")
for(i in 1:20){
  
  cell_num <- gradient_num[i, 1]
  print(cell_num)
  gene_num <- gradient_num[i, 2]
  print(gene_num)
  
  set.seed(10)
  sample_index <- sample(ncol(data), size = cell_num, replace = TRUE)
  set.seed(10)
  gene_index <- sample(nrow(data), size = gene_num, replace = TRUE)
  
  sub_data <- data[gene_index, sample_index]
  group <- as.numeric(factor(group_condition[sample_index]))
  rownames(sub_data) <- paste0("Gene", 1:nrow(sub_data))
  colnames(sub_data) <- paste0("Cell", 1:ncol(sub_data))
  
  for(n in 1:3){
    scala_result <- purrr::map(
      .x = fourteenth_class,
      .f = function(method){
        
        ### estimation
        est <- simpipe::estimate_parameters(ref_data = sub_data,
                                            method = method,
                                            other_prior = list(group.condition = group),
                                            seed = 111,
                                            verbose = TRUE,
                                            use_docker = FALSE)
        time <- lapply(est, function(x){x[["estimate_detection"]][1,2]})
        method_name <- stringr::str_split(names(time), "_", simplify = TRUE)[, 2]
        est_time <- as.numeric(time)
        est_memory <- as.numeric(lapply(est, function(x){x[["estimate_detection"]][1,4]}))
        
        ### simulation
        sim <- simpipe::simulate_datasets(parameters = est,
                                          other_prior = list(nCells = cell_num,
                                                             nGenes = gene_num),
                                          seed = 111,
                                          return_format = "list",
                                          verbose = TRUE)
        
        sim_time <- as.numeric(lapply(sim, function(x){x[["simulate_detection"]][1,2]}))
        sim_memory <- as.numeric(lapply(sim, function(x){x[["simulate_detection"]][1,4]}))
        
        tibble::tibble("method" = method,
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
                                 "class14_",
                                 cell_num,
                                 "_",
                                 gene_num,
                                 "_", n,
                                 ".rds"))
  }
}


## Fifteenth class of methods (CancerInSilico)
data <- readRDS("../preprocessed_data/data85_cellbench-SC4_luyitian.rds")
data <- t(data$data$counts)

fifteenth_class <- c("CancerInSilico")
for(i in 1:20){
  
  cell_num <- gradient_num[i, 1]
  print(cell_num)
  gene_num <- gradient_num[i, 2]
  print(gene_num)
  
  set.seed(i * 10)
  sample_index <- sample(ncol(data), size = cell_num, replace = TRUE)
  set.seed(i * 10)
  gene_index <- sample(nrow(data), size = gene_num, replace = TRUE)
  
  sub_data <- data[gene_index, sample_index]
  rownames(sub_data) <- paste0("Gene", 1:nrow(sub_data))
  colnames(sub_data) <- paste0("Cell", 1:ncol(sub_data))
  
  for(n in 1:3){
    scala_result <- purrr::map(
      .x = fifteenth_class,
      .f = function(method){
        
        ### estimation
        est <- simpipe::estimate_parameters(ref_data = sub_data,
                                            method = method,
                                            seed = 111,
                                            verbose = TRUE,
                                            use_docker = FALSE)
        time <- lapply(est, function(x){x[["estimate_detection"]][1,2]})
        method_name <- stringr::str_split(names(time), "_", simplify = TRUE)[, 2]
        est_time <- as.numeric(time)
        est_memory <- as.numeric(lapply(est, function(x){x[["estimate_detection"]][1,4]}))
        
        ### simulation
        sim <- simpipe::simulate_datasets(parameters = est,
                                          other_prior = list(nCells = cell_num),
                                          seed = 111,
                                          return_format = "list",
                                          verbose = TRUE)
        
        sim_time <- as.numeric(lapply(sim, function(x){x[["simulate_detection"]][1,2]}))
        sim_memory <- as.numeric(lapply(sim, function(x){x[["simulate_detection"]][1,4]}))
        
        tibble::tibble("method" = method,
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
                                 "class15_",
                                 cell_num,
                                 "_",
                                 gene_num,
                                 "_", n,
                                 ".rds"))
  }
}


## Sixtheenth class of SimBPDD
sixteenth_class <- c("SimBPDD")

for(i in 1:20){
  cell_num <- gradient_num[i, 1]
  print(cell_num)
  gene_num <- gradient_num[i, 2]
  print(gene_num)
  
  set.seed(i)
  sub_data <- scater::mockSCE(ncells = cell_num,
                              ngenes = gene_num)
  sub_data <- SingleCellExperiment::counts(sub_data)
  
  for(n in 1:3){
    scala_result <- purrr::map(
      .x = sixteenth_class,
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
                                 "class16_",
                                 cell_num,
                                 "_",
                                 gene_num,
                                 "_", n,
                                 ".rds"))
  }
}


## seventeenth class of methods which users can custom cell and gene number
method <- c("scDesign3")
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
  rownames(sub_data) <- paste0("Gene", 1:nrow(sub_data))
  colnames(sub_data) <- paste0("Cell", 1:ncol(sub_data))
  
  for(n in 1:3){
    ### estimation
    est <- simpipe::estimate_parameters(ref_data = sub_data,
                                        method = method,
                                        seed = 111,
                                        verbose = TRUE)
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
                                 "class17_",
                                 cell_num,
                                 "_",
                                 gene_num,
                                 "_", n,
                                 ".rds"))
  }
}