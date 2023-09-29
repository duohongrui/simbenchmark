library(tidyverse)
library(tibble)
library(ggrepel)
library(ggpubr)
library(patchwork)

### plot trend
plot_trend <- function(scalability_long_data,
                       step2,
                       feature2,
                       dim,
                       class_palette){
  a <- map(unique(scalability_long_data$method), function(y){
    tmp <- scalability_long_data %>% 
      filter(`step` == step2, `feature` == feature2, get(dim) == 1000, `method` == y)
    if(dim == "cell_num"){
      tmp <- tmp %>% 
        mutate(
          coef = gene_coef,
          trend_class = gene_trend_class,
          x = gene_num
        )
    }else{
      tmp <- tmp %>% 
        mutate(
          coef = cell_coef,
          trend_class = cell_trend_class,
          x = cell_num
        )
    }
    if(y == "SCRIP-GP-commonBCV" |
       y == "SCRIP-GP-trendedBCV" |
       y == "SCRIP-BGP-commonBCV" |
       y == "SCRIP-BGP-trendedBCV"){
      y <- paste0(paste0(str_split(y, pattern = "-", simplify = TRUE)[1, ][1:2], collapse = "-"),
                  "-\n",
                  paste0(str_split(y, pattern = "-", simplify = TRUE)[1, ][3], collapse = "-"))
    }
    p <- tmp %>% 
      ggplot(aes(x = x, y = value))+
      geom_point(color = "white", size = 0.5)+
      stat_function(fun = function(x){
        tmp2 <- tmp %>% 
          pull(coef)
        tmp2 <- tmp2[[1]]
        tmp2[1] + tmp2[2]*log(x) + tmp2[3]*sqrt(x) + tmp2[4]*x + tmp2[5]*x^2 + tmp2[6]*x^3
      }, aes(color = trend_class)) +
      theme(panel.background = element_rect(fill = "black", linewidth = 2),
            panel.grid = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            legend.position = "none",
            title = element_text(size = 8))+
      scale_color_manual("", values = class_palette)+
      ggtitle(label = y, subtitle = tmp$trend_class[1])
    p
  }) %>% setNames(unique(scalability_long_data$method))
  return(a)
}


### model function for every method
scam_function <- function(step,
                          feature,
                          train_data,
                          test_data,
                          model_selected = "scam",
                          return_model = TRUE,
                          ntree = 500,
                          verbose = TRUE){
  methods <- unique(train_data$method)
  id <- paste0(step, "_", feature)
  predict_result <- map(methods, .f = function(x){
    if(x == "SPsimSeq" & step == "estimation" |
       x == "scDesign" & step == "estimation" | 
       x == "SimBPDD" & step == "estimation" | 
       x == "scGAN" & feature == "memory"){
      if(verbose){
        message(x)
      }
      cat(paste0("Passing modeling the ", feature, " in ", step, " for ", x, "\n"))
      lst(model = NULL,
          predict_result = NULL)
    }else{
      if(verbose){
        message(x)
      }
      ### training
      tmp <- train_data %>% 
        filter(method == x)
      ### BASiCS and TedSim NA
      if(x == "BASiCS" & feature == "memory" | x == "TedSim"){
        tmp <- tmp %>% 
          drop_na()
      }
      if(model_selected == "scam"){
        model <- scam::scam(get(id) ~ s(cell_num, gene_num, bs = "tedmi"), data = tmp) %>%
          strip::strip(keep = "predict")
      }
      if(model_selected == "RF"){
        model <- randomForest::randomForest(get(id) ~ cell_num + gene_num, data = tmp, ntree = ntree)
      }
      
      ### prediction
      tmp2 <- test_data %>% 
        filter(method == x)
      if(x == "BASiCS" & feature == "memory" | x == "TedSim"){
        tmp2 <- test_data %>% 
          filter(method == x) %>% 
          drop_na()
      }
      predicted_result <- predict(model, tmp2 %>% select(2,3))
      
      predict_result <- tibble("method" = x,
                               "real" = tmp2 %>% pull(id),
                               "predicted" = predicted_result)
      colnames(predict_result) <- c("method",
                                    paste0("real_", id),
                                    paste0("predicted_", id))
      lst(model,
          predict_result)
    }
  }) %>% setNames(unique(train_data$method))
  
  return(predict_result)
}



bind_prediction_result <- function(model){
  tmp_data <- map_dfr(1:length(model), .f = function(x){
    model[[x]][["predict_result"]]
  })
  col_name <- colnames(tmp_data)
  colnames(tmp_data) <- c("method", "real_data", "predict_data")
  tmp_data <- tmp_data %>% 
    filter(predict_data > 0)
  colnames(tmp_data) <- col_name
  tmp_data
}


correlation_plot <- function(table){
  #################### all results
  ### step and feature
  step <- str_split(colnames(table)[2], "_", simplify = TRUE)[2]
  feature <- str_split(colnames(table)[2], "_", simplify = TRUE)[3]
  
  colnames(table) <- c("method", "x", "y")
  table$x <- log10(table$x)
  table$y <- log10(table$y)
  
  cor_value <- cor(table$x, table$y, method = "pearson")
  
  if(feature == "time"){
    feature_id <- "execution time"
  }else{
    feature_id <- "memory usage"
  }
  
  p <- ggplot(table, aes(x = x, y = y)) +
    geom_point(size = 1) +
    stat_smooth(method="lm", se = TRUE) +
    stat_cor(method = "pearson") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour = "black"),
          panel.border = element_rect(colour = "black"),
          axis.ticks = element_line(colour = "black")) +
    xlab(paste0("Actual log10 ", feature_id, " in ", step)) +
    ylab(paste0("Predicted log10 ", feature_id, " in ", step))
  
  #################### for every method
  
  method_prediction_result <- map(unique(table$method), .f = function(z){
    ### method data
    tmp <- table %>% 
      filter(method == z)
    
    ### correlation
    cor_value_method <- round(cor(tmp$x, tmp$y), digits = 2)
    
    ### Plot
    p_method <- ggplot(tmp, aes(x = x, y = y)) +
      geom_point() +
      stat_smooth(method="lm", se = TRUE) +
      stat_cor(method = "pearson") +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text = element_text(colour = "black"),
            panel.border = element_rect(colour = "black"),
            axis.ticks = element_line(colour = "black")) +
      xlab(paste0("Actual log10 ", feature_id, " in ", step)) +
      ylab(paste0("Predicted log10 ", feature_id, " in ", step))
    
    list(cor_value = cor_value_method,
         plot = p_method)
    
  }) %>% setNames(unique(table$method))
  
  ### return all results
  list(correlation = cor_value,
       correlation_plot = p,
       method_prediction_result = method_prediction_result)
}


get_all_correlations <- function(correlation, method_names){
  values <- correlation[["correlation"]]
  for(x in method_names){
    values <- append(values, ifelse(is.null(correlation[["method_prediction_result"]][[x]][[1]]), NA, correlation[["method_prediction_result"]][[x]][[1]]))
  }
  names(values) <- c("correlation", method_names)
  return(values)
}


RF_gradient_function <- function(ntree, methods){
  RF_gradient_result <- map_dfr(ntree, function(tree){
    message("-----------------------------------")
    message(tree)
    message("-----------------------------------")
    ## iterate 10 times
    data_tmp <- list()
    for(n in 1:10){
      cat("-----------------------------------\n")
      cat(paste0(n, "\n"))
      cat("-----------------------------------\n")
      set.seed(n)
      train_data <- data %>% 
        group_by("method") %>% 
        sample_frac(0.8) %>% 
        ungroup()
      test_data <- data[data$row %in% setdiff(1:nrow(data), train_data$row), ]

      #### time of estimation
      message("time of estimation")
      est_time_correlation <- scam_function(step = "estimation",
                                            feature = "time",
                                            train_data = train_data,
                                            test_data = test_data,
                                            model_selected = "RF",
                                            ntree = tree,
                                            verbose = FALSE) %>% 
        bind_prediction_result() %>% 
        correlation_plot() %>% 
        get_all_correlations(method_names = methods)
      #est_time_correlation <- append(anno, est_time_correlation)
      
      #### memory of estimation
      message("memory of estimation")
      est_memory_correlation <- scam_function(step = "estimation",
                                              feature = "memory",
                                              train_data = train_data,
                                              test_data = test_data,
                                              model_selected = "RF",
                                              ntree = tree,
                                              verbose = FALSE) %>% 
        bind_prediction_result() %>% 
        correlation_plot() %>% 
        get_all_correlations(method_names = methods)
      #est_memory_correlation <- append(anno, est_memory_correlation)
      
      #### time of simulation
      message("time of simulation")
      sim_time_correlation <- scam_function(step = "simulation",
                                            feature = "time",
                                            train_data = train_data,
                                            test_data = test_data,
                                            model_selected = "RF",
                                            ntree = tree,
                                            verbose = FALSE) %>% 
        bind_prediction_result() %>% 
        correlation_plot() %>% 
        get_all_correlations(method_names = methods)
      #sim_time_correlation <- append(anno, sim_time_correlation)
      
      #### memory of simulation
      message("memory of simulation")
      sim_memory_correlation <- scam_function(step = "estimation",
                                              feature = "memory",
                                              train_data = train_data,
                                              test_data = test_data,
                                              model_selected = "RF",
                                              ntree = tree,
                                              verbose = FALSE) %>% 
        bind_prediction_result() %>% 
        correlation_plot() %>% 
        get_all_correlations(method_names = methods)
      #sim_memory_correlation <- append(anno, sim_memory_correlation)
      tmp <- tibble(
        step = c("estimation", "estimation", "simulation", "simulation"),
        feature = c("time", "memory", "time", "memory"),
        tree = tree,
        n = n
      )
      tmp2 <- rbind(est_time_correlation,
                    est_memory_correlation,
                    sim_time_correlation,
                    sim_memory_correlation)
      data_tmp[[n]] <- cbind(tmp, tmp2) %>% as_tibble()
    }
    data_tmp <- map_dfr(data_tmp, .f = function(x){x})
    data_tmp
  })
  return(RF_gradient_result)
}


RF_tree_plot_function <- function(Step, Feature, Methods, line_data, RF_gradient_result){
  
  boxplot_data <- RF_gradient_result %>% filter(step == Step,
                                                feature == Feature)
  plot_list <- map(methods, .f = function(x){
    ggplot() +
      geom_boxplot(data = boxplot_data,
                   aes(x = tree, y = get(x), fill = tree),
                   alpha = 0.6,
                   size = 0.25,
                   color = "#176691",
                   outlier.alpha = 0,
                   width = 0.5) +
      geom_jitter(data = boxplot_data,
                  aes(x = tree, y = get(x)),
                  inherit.aes = TRUE,
                  size = 0.02,
                  color = "gray30") +
      geom_line(data = line_data %>% filter(step == Step,
                                            feature == Feature,
                                            method == x),
                mapping = aes(x = tree,
                              y = mean_correlation),
                color = "#E41A1C",
                lwd = 0.25) +
      theme_pubr() +
      scale_fill_manual(values = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "YlGnBu"))(20)) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, size = 3),
        axis.text.y = element_text(size = 3),
        axis.title = element_text(size = 4),
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.length = unit(0.05, units = "cm")
      )+
      ylab(paste0("correlation (", x, ")"))
  }) %>% setNames(methods)
  return(plot_list)
}
