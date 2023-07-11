library(tidyverse)
library(ggpubr)
library(patchwork)
library(ggrepel)

color <- c("#F1947C", "#6CB9B0")
method_class_colors <- c("#F97369", "#69ACD1", "#F9A24B", "#ADDD52", "#F9B4DA")

function_metric_for_datasets <- function(function_data,
                                         functionality,
                                         bar_plot_name,
                                         cor_plot_name,
                                         technique_bar_name){
  data_info <- openxlsx::read.xlsx("./Chunk1-Data preparation/evaluation_datasets.xlsx", rowNames = TRUE) %>% 
    mutate(
      Data = str_split(Dataset.ID, "_", simplify = TRUE)[, 1]
    ) %>% 
    select(Data, Platform)
  methods_info <- openxlsx::read.xlsx("Chunk1-Data preparation/methods.xlsx") %>% 
    select(1, 3, 6:9) %>% 
    mutate(
      color = case_when(
        Category == "Class 1" ~ method_class_colors[1],
        Category == "Class 2" ~ method_class_colors[2],
        Category == "Class 3" ~ method_class_colors[3],
        Category == "Class 4" ~ method_class_colors[4],
        Category == "Class 5" ~ method_class_colors[5]
      )
    )
  
  platforms <- c("MARS-Seq",
                 "10X Genomics",
                 "Smart-seq2",
                 "Mix sources1",
                 "CEL-seq",
                 "Fluidigm C1",
                 "CEL-seq2",
                 "Mix sources2",
                 "Drop-seq",
                 "inDrop",
                 "Microwell-seq",
                 "Smart-seq",
                 "ST",
                 "HDST",
                 "10X Visium",
                 "Slide-Seq",
                 "Slide-SeqV2",
                 "seqFISH",
                 "seqFISH+",
                 "osmFISH",
                 "sci-Space",
                 "MERFISH",
                 "Stereo-Seq")
  
  ################################ Platforms ###################################
  platforms_data <- function_data %>% 
    group_by(Method, Data) %>% 
    summarise(
      value = mean(value, na.rm = TRUE)
    ) %>% 
    ungroup() %>% 
    left_join(data_info, by = "Data") %>% 
    relocate(Platform, .after = "Data") %>% 
    mutate(
      Type = case_when(
        Data %in% paste0("data", 1:101) ~ "scRNA-seq data",
        Data %in% paste0("data", 101:152) ~ "spatial transcriptome data"
      )
    ) %>% 
    relocate(Type, .after = "Platform") %>% 
    group_by(Method, Type) %>% 
    summarise(
      value = mean(value, na.rm = TRUE)
    ) %>% 
    ungroup() %>% 
    left_join(methods_info, by = "Method")
  
  platforms_data$Method <- factor(platforms_data$Method, levels = methods_info$Method[which(methods_info$Method %in% unique(platforms_data$Method))])
  platforms_data <- platforms_data %>% 
    arrange(desc(value))
  
  platform_bar <- ggplot(platforms_data,
                         aes(x = Method, y = value, fill = Category))+
    geom_bar(stat = "identity",
             width = 0.7) +
    theme_bw() +
    facet_grid(Type ~ .) +
    scale_fill_manual(values = method_class_colors) +
    theme(
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6.5),
      axis.text.y = element_text(size = 6.5),
      axis.title.y = element_text(size = 8),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6),
      legend.position = "bottom",
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.box = "horizontal",
      legend.key = element_rect(fill = NA),
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      strip.text = element_text(size = 5)
    )
  
  ggsave(plot = platform_bar,
         filename = paste0("/Users/duohongrui/Desktop/sim-article/figures/", bar_plot_name),
         width = 12,
         height = 10,
         units = "cm")
  
  if(functionality == "Group"){
    cor_data <- platforms_data %>% 
      pivot_wider(names_from = "Type", values_from = "value")
    cor_data <- cor_data %>% 
      mutate(
        across(all_of(c("scRNA-seq data", "spatial transcriptome data")), ~ replace_na(.x, 0))
      ) %>% 
      filter(`scRNA-seq data` != 0, `spatial transcriptome data` != 0)
    cor_value <- sprintf("%0.2f", cor(cor_data$`scRNA-seq data`, cor_data$`spatial transcriptome data`))
    
    fig5b <- ggplot(cor_data, aes(x = `spatial transcriptome data`, y = `scRNA-seq data`)) +
      geom_abline(intercept = 0, slope = 1, color = "red") +
      geom_point(size = 1) +
      geom_label_repel(aes(label = Method,
                           fill = Category),
                       color = "white",
                       size = 1.5,
                       box.padding = unit(0.5, "lines"),
                       point.padding = unit(0.5, "lines"),
                       segment.colour = "black",
                       max.overlaps = 50) +
      annotate(geom = "text",
               x = 0.25,
               y = 0.8,
               label = paste0("Cor: ", cor_value),
               size = 3) +
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
      scale_fill_manual(values = method_class_colors)
    
    ggsave(plot = fig5b,
           filename = paste0("/Users/duohongrui/Desktop/sim-article/figures/", cor_plot_name),
           width = 9,
           height = 9,
           units = "cm")
  }
  
  
  
  ############################### Techniques ###################################
  techniques_data <- function_data %>% 
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
  
  techniques_data$Method <- factor(techniques_data$Method, levels = methods_info$Method[which(methods_info$Method %in% unique(techniques_data$Method))])
  techniques_data$Platform <- factor(techniques_data$Platform, levels = platforms)
  technique_bar <- ggbarplot(data = techniques_data,
                             x = "Method",
                             y = "value",
                             color = "Category",
                             fill = "Category",
                             ggtheme = theme_pubr(),
                             size = 1,
                             alpha = 0.6)+
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
          strip.text.y.right = element_text(size = 8, angle = 0))
  
  ggsave(plot = technique_bar,
         filename = paste0("/Users/duohongrui/Desktop/sim-article/figures/", technique_bar_name),
         width = 18,
         height = 25,
         units = "cm")
  
}


