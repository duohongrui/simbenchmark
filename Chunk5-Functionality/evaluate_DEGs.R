library(tibble)
data_list <- list.files("G:/sim_bench/simulation_data/", pattern = "^Lun")

for(i in data_list[44:90]){
  data <- readRDS(file.path("G:/sim_bench/simulation_data", i))
  if(data$sim_data_info$group >= 2 & "de_gene" %in% colnames(data$sim_data$row_meta)){
    message(i)
    
    sim_data <- as.matrix(data$sim_data$count_data)
    group <- as.character(data$sim_data$col_meta$group)
    if(rlang::is_empty(group)){
      group <- paste0("Group", as.numeric(data$sim_data$col_meta$plates))
    }
    group_combn <- utils::combn(unique(group), 2)
    de_genes <- rownames(sim_data)[which(data[["sim_data"]][["row_meta"]][["de_gene"]] == "yes")]
    sim_DEGs <- list()
    valid_DEGs_distribution <- list()
    distribution_score <- c()
    
    if(stringr::str_starts(i, pattern = "Splat") | stringr::str_starts(i, pattern = "SCRIP")){
      for(conb in 1:ncol(group_combn)){
        conb1 <- group_combn[1, conb]
        conb2 <- group_combn[2, conb]
        conb_name <- paste0(group_combn[, conb], collapse = "vs")
        message(conb_name)
        fac1 <- data[["sim_data"]][["row_meta"]][, stringr::str_ends(colnames(data[["sim_data"]][["row_meta"]]), pattern = conb1)]
        fac2 <- data[["sim_data"]][["row_meta"]][, stringr::str_ends(colnames(data[["sim_data"]][["row_meta"]]), pattern = conb2)]
        index <- fac1 != fac2
        DEGs <- rownames(data$sim_data$count_data)[index]
        sim_DEGs[[conb_name]] <- DEGs
        
        ### Distribution
        message("Distribution of null data...")
        col1 <- which(data[["sim_data"]][["col_meta"]][["group"]] %in% conb1)
        col2 <- which(data[["sim_data"]][["col_meta"]][["group"]] %in% conb2)
        sub_data <- sim_data[!index, c(col1, col2)]
        sub_group <- c(rep(conb1, length(col1)), rep(conb2, length(col2)))
        error <- try(sub_DEA_result <- simutils::perform_DEA(data = sub_data,
                                                             group = sub_group,
                                                             method = "edgeRQLFDetRate",
                                                             verbose = TRUE))
        if(class(error) == "try-error"){
          distribution_score <- append(distribution_score, NA)
          valid_DEGs_distribution[[conb_name]] <- NA
        }else{
          valid_DEGs_distribution[[conb_name]] <- sub_DEA_result[[1]]
          p_values <- sub_DEA_result[[1]][["PValue"]]
          uniform_result <- simutils::test_uni_distribution(p_values)
          distribution_score <- append(distribution_score, uniform_result[["score"]])
        }
      }
      distribution_score <- mean(distribution_score)
      saveRDS(valid_DEGs_distribution, file.path("G:/sim_bench/valid_DEGs_distribution", i))
    }
    
    ### True proportions of DEGs
    message("True proportions of DEGs...")
    error2 <- try(DEGs_result <- simutils::true_DEGs_proportion(sim_data = sim_data,
                                                                group = group,
                                                                group_combn = group_combn,
                                                                sim_DEGs = sim_DEGs,
                                                                DEA_method = "edgeRQLFDetRate",
                                                                verbose = TRUE))
    if(class(error2) == "try-error"){
      true_proportion <- NA
    }else{
      saveRDS(DEGs_result, file.path("G:/sim_bench/DEGs_result", i))
      true_proportion <- DEGs_result[["weighted_true_prop"]]
    }
    
    ### SVM
    message("SVM...")
    error3 <- try(SVM_result <- simutils::model_predict(data = sim_data,
                                                        group = group,
                                                        de_genes = de_genes,
                                                        method = "SVM",
                                                        verbose = TRUE))
    if(class(error3) == "try-error"){
      AUC <- NA
      Accuracy <- NA
      Precision <- NA
      Recall <- NA
      F1 <- NA
    }else{
      saveRDS(SVM_result, file.path("G:/sim_bench/SVM_model", i))
      AUC <- as.numeric(SVM_result$roc$auc)
      Accuracy <- unname(SVM_result[["conf_matrix"]][["overall"]][1])
      if(length(unique(group)) == 2){
        Precision <- unname(SVM_result[["conf_matrix"]][["byClass"]]["Precision"])
        Recall <- unname(SVM_result[["conf_matrix"]][["byClass"]]["Recall"])
        F1 <- unname(SVM_result[["conf_matrix"]][["byClass"]]["F1"])
      }else{
        Precision <- mean(SVM_result[["conf_matrix"]][["byClass"]][, "Precision"], na.rm = TRUE)
        Recall <- mean(SVM_result[["conf_matrix"]][["byClass"]][, "Recall"], na.rm = TRUE)
        F1 <- mean(SVM_result[["conf_matrix"]][["byClass"]][, "F1"], na.rm = TRUE)
      }
    }
    
    saveRDS(dplyr::lst(distribution_score,
                       true_proportion,
                       Accuracy,
                       Precision,
                       Recall,
                       F1),
            file.path("G:/sim_bench/DEGs_evaluation", i))
    
  }else{
    next
  }
}
