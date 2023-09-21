library(tidyverse)
library(ggpubr)
library(patchwork)
library(ggrepel)

technique_colors <- c("#7DAEF4", "#E0709E")
method_class_colors <- c(RColorBrewer::brewer.pal(12, "Set3")[4],
                         "#94D3C7","#AAAAE0","#ADCE65","#EAC5E3")


data_process <- function(function_data,
                         functionality,
                         info){
  data_info <- openxlsx::read.xlsx("./Chunk1-Data preparation/evaluation_datasets.xlsx", rowNames = TRUE) %>% 
    mutate(
      Data = str_split(Dataset.ID, "_", simplify = TRUE)[, 1]
    ) %>% 
    select(Data, Platform)
  methods_info <- openxlsx::read.xlsx("Chunk1-Data preparation/methods.xlsx") %>% 
    select(1, 3, 7:10) %>% 
    mutate(
      color = case_when(
        Category == "Class 1" ~ method_class_colors[1],
        Category == "Class 2" ~ method_class_colors[2],
        Category == "Class 3" ~ method_class_colors[3],
        Category == "Class 4" ~ method_class_colors[4],
        Category == "Class 5" ~ method_class_colors[5]
      )
    )
  
  platforms <- c("MARS-seq",
                 "10X Genomics",
                 "Smart-seq2",
                 "Mix sources1",
                 "CEL-seq",
                 "Fluidigm C1",
                 "CEL-seq2",
                 "Mix sources2",
                 "Drop-seq",
                 "inDrop",
                 "STRT-seq",
                 "Microwell-seq",
                 "Smart-seq",
                 "ST",
                 "HDST",
                 "10X Visium",
                 "Slide-seq",
                 "Slide-seqV2",
                 "seqFISH",
                 "seqFISH+",
                 "osmFISH",
                 "sci-Space",
                 "MERFISH",
                 "Stereo-seq")
  
  if(info == "platform"){
    platforms_data <- function_data %>% 
      group_by(Method, Data) %>% 
      summarise(
        value = mean(value, na.rm = TRUE)
      ) %>% 
      ungroup() %>% 
      left_join(data_info, by = "Data") %>% 
      relocate(Platform, .after = "Data") %>% 
      mutate(
        Platform = case_when(
          Platform == "Smart-seq2\r\n10X Genomics" ~ "Mix sources1",
          Platform == "CEL-seq\r\nCEL-seq2" ~ "Mix sources2",
          TRUE ~ Platform
        )
      ) %>% 
      relocate(Platform, .after = "Platform") %>% 
      group_by(Method, Platform) %>% 
      summarise(
        value = mean(value, na.rm = TRUE)
      ) %>% 
      ungroup() %>% 
      left_join(methods_info, by = "Method")
    
    platforms_data$Platform <- factor(platforms_data$Platform, levels = platforms)
    platforms_data$Method <- factor(platforms_data$Method,
                                    levels = methods_info$Method[which(methods_info$Method %in% unique(platforms_data$Method))])
    platforms_data$Category <- factor(platforms_data$Category, levels = paste0("Class ", 1:5))
    
    if(functionality == "Group"){
      platforms_data <- platforms_data %>% 
        filter(Platform != "Mix sources2",
               Platform != "Microwell-seq",
               Platform != "MERFISH",
               Platform != "Stereo-Seq")
    }
    return(platforms_data)
  }
  
  if(info == "technique"){
    techniques_data <- function_data %>% 
      group_by(Method, Data) %>% 
      summarise(
        value = mean(value, na.rm = TRUE)
      ) %>% 
      ungroup() %>% 
      left_join(data_info, by = "Data") %>% 
      relocate(Platform, .after = "Data") %>% 
      mutate(
        Type = case_when(
          Data %in% paste0("data", 1:101) ~ "scRNA-seq",
          Data %in% paste0("data", 101:152) ~ "ST technology"
        )
      ) %>% 
      relocate(Type, .after = "Platform") %>% 
      group_by(Method, Type) %>% 
      summarise(
        value = mean(value, na.rm = TRUE)
      ) %>% 
      ungroup() %>% 
      left_join(methods_info, by = "Method")
    ### Method factor
    techniques_data$Method <- factor(techniques_data$Method,
                                     levels = methods_info$Method[which(methods_info$Method %in% unique(techniques_data$Method))])
    ### Sorting
    techniques_data <- techniques_data %>% 
      arrange(desc(value))
    return(techniques_data)
  }
}




boxplot_platform <- function(function_data,
                             functionality){
  platforms_data <- data_process(function_data = function_data,
                                 functionality = functionality,
                                 info = "platform")
  plot <- ggboxplot(data = platforms_data,
                    x = "Category",
                    y = "value",
                    color = "Category",
                    ggtheme = theme_pubr(),
                    size = 0.2,
                    alpha = 0.6,
                    width = 0.8,
                    bxp.errorbar = FALSE,
                    outlier.size = 0,
                    add = "jitter",
                    add.params = list(size = 0.15))+
    ylab("Metric values") +
    scale_color_manual(values = method_class_colors) +
    scale_fill_manual(values = method_class_colors) +
    facet_wrap(.~ Platform, ncol = 4, strip.position = "top") +
    theme(axis.text = element_text(size = 4),
          axis.text.x = element_text(angle = 30, hjust = 1),
          axis.title = element_text(size = 5),
          axis.title.x = element_blank(),
          axis.ticks = element_line(linewidth = 0.1),
          axis.ticks.length = unit(0.05, 'cm'),
          axis.line = element_line(linewidth = 0.2),
          legend.text = element_text(size = 4),
          legend.position = c(0.5, 1.05),
          legend.direction = "horizontal",
          legend.background = element_blank(),
          legend.title = element_text(size = 4),
          legend.key.size = unit(0.5, 'cm'),
          strip.background = element_rect(fill = "white", color = "white"),
          strip.text = element_text(size = 4),
          strip.text.x = element_text(margin = margin(-0.01,0,-0.01,0, "cm")),
          plot.margin = unit(c(0.6,0.1,0.1,0.1), 'cm')) +
    ylab("Platform scores")
  return(plot)
  
}



function_metric_for_datasets <- function(function_data,
                                         functionality){
  ################################ Techniques ###################################
  techniques_data <- data_process(function_data = function_data,
                                  functionality = functionality,
                                  info = "technique")
  
  tmp <- techniques_data %>% 
    pivot_wider(names_from = "Type", values_from = "value") %>% 
    tibble::column_to_rownames("Method")
  
  method_color <- tmp[levels(techniques_data$Method), ]$color
  
  ### Barplot
  technique_bar <- ggplot(techniques_data)+
    geom_bar(aes(x = Method, y = value, fill = Type),
             stat = "identity",
             position = "dodge",
             width = 0.6) +
    theme_bw() +
    scale_fill_manual(values = technique_colors) +
    theme(
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 45,
                                 vjust = 1,
                                 hjust = 1,
                                 size = 5.5,
                                 color = "black"),
      axis.text.y = element_text(size = 4),
      axis.title.y = element_text(size = 5),
      legend.title = element_text(size = 5),
      legend.text = element_text(size = 4),
      legend.position = "bottom",
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.box = "horizontal",
      legend.key = element_rect(fill = NA),
      legend.key.height = unit(0.4, "cm"),
      legend.key.width = unit(0.4, "cm"),
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks = element_line(linewidth = 0.2, lineend = 0.3)
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    ylab(paste0("Scores of ", tolower(functionality), " function"))

  
  
  if(functionality != "Batch"){
    ### linear regression plot
    cor_data <- techniques_data %>% 
      pivot_wider(names_from = "Type", values_from = "value")
    cor_data <- cor_data %>% 
      mutate(
        across(all_of(c("scRNA-seq", "ST technology")), ~ replace_na(.x, 0))
      ) %>% 
      filter(`scRNA-seq` != 0, `ST technology` != 0)
    

    if(functionality == "Trajectory"){
      color <- method_class_colors[2]
    }else{
      color <- method_class_colors
    }
    
    fig <- ggplot(cor_data, aes(x = `ST technology`, y = `scRNA-seq`)) +
      geom_smooth(method = 'lm', linewidth = 0.7, color = "red") +
      geom_point(size = 0.7) +
      stat_cor(method = "pearson", size = 2) +
      geom_label_repel(aes(label = Method,
                           fill = Category),
                       color = "white",
                       size = 1.3,
                       box.padding = unit(0.4, "lines"),
                       label.padding = unit(0.15, "lines"),
                       point.padding = unit(0.3, "lines"),
                       segment.colour = "black",
                       max.overlaps = 50) +
      theme_bw() +
      theme(legend.position = "top",
            axis.text = element_text(size = 6, color = "black"),
            axis.title = element_text(size = 6),
            legend.title = element_text(size = 6),
            legend.text = element_text(size = 5),
            panel.grid = element_blank(),
            legend.key.width = unit(0.3, units = "cm"),
            legend.key.height = unit(0.3, units = "cm")) +
      xlab(paste0(functionality, " scores on the spatial datasets")) +
      ylab(paste0(functionality, " scores on the scRNA-seq datasets")) +
      scale_fill_manual(values = color)
    
  }
  
  
  
  ############################### Platforms ###################################
  platforms_data <- data_process(function_data = function_data,
                                 functionality = functionality,
                                 info = "platform")
  ### top 5 method be labelled
  top_5_methods <- platforms_data %>% 
    group_by(Platform) %>% 
    top_n(., 5, value) %>% 
    arrange(desc(value)) %>% 
    mutate(
      n = 1:n()
    ) %>% 
    ungroup() %>% 
    mutate(
      n = case_when(
        n == 6 ~ 5,
        TRUE ~ n
      )
    )
  
  ### bottom 5 method be labelled
  bottom_5_methods <- platforms_data %>% 
    mutate(
      value = 1 - value
    ) %>% 
    group_by(Platform) %>% 
    top_n(., 5, value) %>% 
    arrange(desc(value)) %>% 
    mutate(
      n = 1:n()
    ) %>% 
    ungroup() %>% 
    mutate(
      n = case_when(
        n == 6 ~ 5,
        TRUE ~ n
      ),
      n = -n
    )
  platform_bar <- ggbarplot(data = platforms_data,
                             x = "Method",
                             y = "value",
                             color = "Category",
                             fill = "Category",
                             ggtheme = theme_pubr(),
                             size = 1,
                             alpha = 0.6) +
    geom_text(data = top_5_methods,
              mapping = aes(x = Method, y = value, label = n), nudge_y = -0.2, size = 2) +
    geom_text(data = bottom_5_methods,
              mapping = aes(x = Method, y = value, label = n), nudge_y = -0.01, size = 2) +
    ylab("Metric values") +
    scale_color_manual(values = method_class_colors) +
    scale_fill_manual(values = method_class_colors) +
    facet_grid(Platform ~ .) +
    ylim(0, 1) +
    scale_y_continuous(labels = c("0", "0.8"), breaks = c(0, 0.8)) +
    theme(axis.text = element_text(size = 9),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9),
          legend.position = "top",
          legend.background = element_blank(),
          axis.title.x = element_blank(),
          strip.background = element_rect(fill = NA, color = "black"),
          strip.text.y.right = element_text(size = 8, angle = 0)) +
    ylab(paste0("Scores of ", tolower(functionality), " function"))
  
  if(functionality != "Batch"){
    return(list(technique_bar = technique_bar,
                fig = fig,
                platform_bar = platform_bar))
  }else{
    return(list(technique_bar = technique_bar,
                platform_bar = platform_bar))
  }
  
}


