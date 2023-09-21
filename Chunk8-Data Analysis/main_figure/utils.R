accuracy_process_function <- function(accuracy){
  
  ### summarize accuracy score for quantification strategy
  accuracy_summary_quantification <- accuracy %>% 
    group_by(Method, metric, Platform, Type, Data, `Quantification Strategy`) %>% 
    summarise(
      value = mean(value, na.rm = TRUE),
    ) %>% 
    group_by(Method, metric, `Quantification Strategy`, Type, Data) %>% 
    summarise(
      value = mean(value, na.rm = TRUE),
    ) %>% 
    filter(Type == "scRNA-seq data") %>% 
    group_by(Method, `Quantification Strategy`, metric) %>% 
    summarise(
      value = mean(value, na.rm = TRUE),
    ) %>% 
    group_by(Method, `Quantification Strategy`) %>% 
    summarise(
      value = mean(value, na.rm = TRUE),
    ) %>% 
    drop_na() %>% 
    mutate(
      `Quantification Strategy` = paste0("acc_", `Quantification Strategy`)
    )
  
  
  ### summarize accuracy score for metrics
  accuracy_summary_per_metric <- accuracy %>% 
    group_by(Method, metric, Platform, Type, Data) %>% 
    summarise(
      value = mean(value, na.rm = TRUE),
    ) %>% 
    ungroup()
  
  ### summarize accuracy score for platform
  accuracy_summary_per_platform <- accuracy_summary_per_metric %>% 
    group_by(Method, metric, Platform) %>% 
    summarise(
      value = mean(value, na.rm = TRUE)
    ) %>% 
    group_by(Method, Platform, metric) %>% 
    summarise(
      value = mean(value, na.rm = TRUE)
    ) %>% 
    group_by(Method, Platform) %>% 
    summarise(
      value = mean(value, na.rm = TRUE)
    ) %>% 
    mutate(
      Platform = paste0("acc_", Platform)
    )
  
  ### summarize accuracy score for techniques
  accuracy_summary_per_technique <- accuracy_summary_per_metric %>% 
    group_by(Method, metric, Type, Platform) %>% 
    summarise(
      value = mean(value, na.rm = TRUE)
    ) %>% 
    group_by(Method, Type, metric) %>% 
    summarise(
      value = mean(value, na.rm = TRUE)
    ) %>% 
    group_by(Method, Type) %>% 
    summarise(
      value = mean(value, na.rm = TRUE)
    ) %>% 
    mutate(
      Type = paste0("acc_", Type)
    )
  
  ### summarize accuracy score for metrics
  accuracy_summary_per_metric <- accuracy_summary_per_metric %>% 
    group_by(Method, metric, Data) %>% 
    summarise(
      value = mean(value, na.rm = TRUE)
    ) %>% 
    group_by(Method, metric) %>% 
    summarise(
      value = mean(value, na.rm = TRUE)
    )
  
  ### summarize accuracy score
  accuracy_score <- accuracy_summary_per_metric %>% 
    group_by(Method) %>% 
    summarise(
      accuracy = mean(value, na.rm = TRUE)
    ) %>% 
    ungroup()
  
  ### accuracy_summary_per_property
  accuracy_summary_per_property <- accuracy %>% 
    group_by(Method, property, Data) %>% 
    summarise(
      value = mean(value, na.rm = TRUE)
    ) %>% 
    group_by(Method, property) %>% 
    summarise(
      value = mean(value, na.rm = TRUE)
    )
  
  ### add accuracy score
  accuracy_plot <- accuracy_score %>% 
    full_join(., accuracy_summary_per_metric %>% pivot_wider(names_from = "metric", values_from = "value"), by = "Method") %>% 
    full_join(., accuracy_summary_per_property %>% pivot_wider(names_from = "property", values_from = "value"), by = "Method") %>% 
    full_join(., accuracy_summary_per_platform %>% pivot_wider(names_from = "Platform", values_from = "value"), by = "Method") %>% 
    full_join(., accuracy_summary_per_technique %>% pivot_wider(names_from = "Type", values_from = "value"), by = "Method") %>% 
    full_join(., accuracy_summary_quantification %>% pivot_wider(names_from = "Quantification Strategy", values_from = "value"), by = "Method") %>% 
    arrange(desc(accuracy))
  
  return(accuracy_plot)
}



functionality_process_function <- function(functionality){
  
  ### summarize functionality score for quantification strategy
  functionality_summary_quantification <- functionality %>% 
    group_by(Method, metric, Platform, Type, Data, `Quantification Strategy`) %>% 
    summarise(
      value = mean(value, na.rm = TRUE),
    ) %>% 
    group_by(Method, metric, `Quantification Strategy`, Type, Data) %>% 
    summarise(
      value = mean(value, na.rm = TRUE),
    ) %>% 
    filter(Type == "scRNA-seq data") %>% 
    group_by(Method, `Quantification Strategy`, metric) %>% 
    summarise(
      value = mean(value, na.rm = TRUE),
    ) %>% 
    group_by(Method, `Quantification Strategy`) %>% 
    summarise(
      value = mean(value, na.rm = TRUE),
    ) %>% 
    drop_na() %>% 
    mutate(
      `Quantification Strategy` = paste0("func_", `Quantification Strategy`)
    )
  
  ### summarize functionality score for metrics
  functionality_summary_per_metric <- functionality %>% 
    group_by(Method, Platform, Type, metric, Data) %>% 
    summarise(
      value = mean(value, na.rm = TRUE)
    ) %>% 
    ungroup()
  
  ### summarize functionality score for platforms
  functionality_summary_per_technique <- functionality_summary_per_metric %>% 
    group_by(Method, metric, Platform) %>% 
    summarise(
      value = mean(value, na.rm = TRUE)
    ) %>% 
    ungroup() %>% 
    group_by(Method, Platform) %>% 
    summarise(
      value = mean(value, na.rm = TRUE)
    ) %>% 
    ungroup() %>% 
    mutate(
      Platform = paste0("func_", Platform)
    )
  
  ### summarize functionality score for techniques
  functionality_summary_per_platform <- functionality_summary_per_metric %>% 
    group_by(Method, metric, Type) %>% 
    summarise(
      value = mean(value, na.rm = TRUE)
    ) %>% 
    ungroup() %>% 
    group_by(Method, Type) %>% 
    summarise(
      value = mean(value, na.rm = TRUE)
    ) %>% 
    ungroup() %>% 
    mutate(
      Type = paste0("func_", Type)
    )
  
  
  ### summarize functionality score
  functionality_score <- functionality_summary_per_metric %>% 
    group_by(Method, metric, Data) %>% 
    summarise(
      value = mean(value, na.rm = TRUE)
    ) %>% 
    group_by(Method, metric) %>% 
    summarise(
      value = mean(value, na.rm = TRUE)
    ) %>% 
    group_by(Method) %>% 
    summarise(
      functionality = mean(value, na.rm = TRUE)
    )
  
  ### add functionality score to the metric result
  functionality_summary_per_metric <- functionality_summary_per_metric %>% 
    group_by(Method, metric, Data) %>% 
    summarise(
      value = mean(value, na.rm = TRUE)
    ) %>% 
    group_by(Method, metric) %>% 
    summarise(
      value = mean(value, na.rm = TRUE)
    ) %>% 
    pivot_wider(names_from = "metric", values_from = "value")
  col_names <- colnames(readRDS("./Chunk8-Data Analysis/functionality/functionality_data.rds"))[c(-1, -2)]
  functionality_summary_per_metric <- functionality_summary_per_metric[, c("Method", col_names)]
  
  functionality_plot <- functionality_score %>% 
    full_join(., functionality_summary_per_metric, by = "Method") %>% 
    full_join(., functionality_summary_per_platform %>% pivot_wider(names_from = "Type", values_from = "value"), by = "Method") %>% 
    full_join(., functionality_summary_per_technique %>% pivot_wider(names_from = "Platform", values_from = "value"), by = "Method") %>% 
    full_join(., functionality_summary_quantification %>% pivot_wider(names_from = "Quantification Strategy", values_from = "value"), by = "Method") %>% 
    arrange(desc(functionality))
  
  functionality_plot <- functionality_plot %>% 
    mutate(
      Group_score = apply(X = functionality_plot %>% select(3:8), MARGIN = 1, FUN = function(x){mean(x, na.rm = TRUE)}),
      DEGs_score = apply(X = functionality_plot %>% select(9:15), MARGIN = 1, FUN = function(x){mean(x, na.rm = TRUE)}),
      Batch_score = apply(X = functionality_plot %>% select(16:22), MARGIN = 1, FUN = function(x){mean(x, na.rm = TRUE)}),
      Trajectory_score = apply(X = functionality_plot %>% select(23:26), MARGIN = 1, FUN = function(x){mean(x, na.rm = TRUE)})
    ) %>% 
    relocate(Group_score, .before = "CDI") %>% 
    relocate(DEGs_score, .before = "distribution_score") %>% 
    relocate(Batch_score, .before = "cms") %>% 
    relocate(Trajectory_score, .before = "HIM")
  return(functionality_plot)
}