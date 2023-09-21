library(ggplot2)
library(tidyverse)
library(ggdendro)
library(dendextend)
library(patchwork)

cor_plot_dendro <- function(data, 
                            type = "platform",
                            criteria = "accuracy",
                            column_names = NULL,
                            k = 3,
                            branch_lwd = 0.6,
                            class_title = FALSE,
                            title = NULL){
  
  if(type == "platform"){
    items <- c("MARS-Seq","10X Genomics","Smart-seq2",
               "Mix sources1","CEL-seq","Fluidigm C1",
               "CEL-seq2","Mix sources2","Drop-seq",
               "inDrop","Microwell-seq","Smart-seq",
               "ST","HDST","10X Visium","Slide-Seq",
               "Slide-SeqV2","seqFISH","seqFISH+",
               "osmFISH","sci-Space","MERFISH",
               "Stereo-Seq")
    data <- data %>% 
      select(all_of(c("id", "Category", items)))
  }else{
    if(is.null(column_names)){
      if(any(str_detect("id", colnames(data)))){
        stop("Please input method id")
      }
      if(any(str_detect("Category", colnames(data)))){
        stop("Please input category")
      }
    }else{
      data <- data %>% 
        select(all_of(c("id", "Category", column_names)))
    }
  }
  
  ### correlation
  value_dataframe <- data %>% select(c(-1,-2))
  value_dataframe <- value_dataframe[!rowSums(value_dataframe, na.rm = TRUE)==0,
                                     !colSums(value_dataframe, na.rm = TRUE)==0]
  cor_res <- Hmisc::rcorr(value_dataframe %>% as.matrix(), type = "pearson")
  cor_mat <- cor_res$r
  p_mat <- cor_res$P
  
  if(criteria == "accuracy" & !class_title){
    ### hclust
    clust <- hclust(dist(t(value_dataframe %>% as.matrix())), method = "ward.D")
    
    order <- rev(clust$labels[clust$order])
    
    ### correlation
    cor_mat <- cor_res$r[order, order]
    p_mat <- cor_res$P[order, order]
    
    ### dendrogram
    ggd <- as.ggdend(clust %>%
                       as.dendrogram() %>%
                       set("branches_k_color", k = k) %>% 
                       set("branches_lwd", branch_lwd))
    
    p1 <- ggplot(ggd, horiz = TRUE, theme = NULL, labels = FALSE) +
      theme_dendro()
  }else{
    order <- colnames(value_dataframe)
  }
  
  ### correlation data
  cor_mat_up <- cor_mat
  cor_mat_up[lower.tri(cor_mat_up, diag = T)] <- NA
  
  cor_mat_up_long <- cor_mat_up %>%
    as.data.frame() %>%
    mutate(x = factor(rownames(cor_mat_up), levels = order)) %>%
    pivot_longer(cols = !x, names_to = "y", values_to = "cor") %>%
    mutate(y = factor(y, levels = rev(order)))
  
  # P values
  p_mat_lower <- p_mat
  p_mat_lower[!lower.tri(p_mat_lower, diag = F)] <- 1
  
  p_mat_lower_long <- p_mat_lower %>%
    as.data.frame() %>%
    mutate(x = factor(rownames(p_mat_lower), levels = order)) %>%
    pivot_longer(cols = !x, names_to = "y", values_to = "p") %>%
    mutate(y = factor(y, levels = rev(order)))
  
  p_mat_lower_long$p_sig <- as.character(symnum(p_mat_lower_long$p,
                                                cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                                symbols = c("***", "**", "*", "")))
  p_mat_lower_long$p[which(p_mat_lower_long$p == 0)] <- sort(p_mat_lower_long$p)[which(sort(p_mat_lower_long$p) != 0)][1]
  p_mat_lower_long$p_sig[which(p_mat_lower_long$p_sig == "?")] <- ""
  
  p2 <- ggplot() +
    geom_point(data = p_mat_lower_long,
               aes(x, y, color = -log10(p)),
               size = 9,
               shape = 15)+
    geom_text(data = p_mat_lower_long,
              aes(x, y, label = p_sig))+
    scale_color_gradientn(colours = c("#ffffff", "#C8D39B", "#2D9042", "#FFE41F"))+
    geom_point(data = cor_mat_up_long,
               shape = 21,
               aes(x, y, size = cor, fill = cor),
               color = "gray50")+
    scale_fill_gradient2(high = 'orange', mid = 'gray92',low = 'navyblue', midpoint = mean(c(max(as.numeric(cor_mat), na.rm = TRUE),
                                                                                             min(as.numeric(cor_mat), na.rm = TRUE))))+
   scale_x_discrete(expand = c(0.025, 0.025))+
   scale_y_discrete(expand = c(0.025, 0.025))+
    geom_vline(aes(xintercept = seq(0.5, (0.5 + ncol(value_dataframe)), 1)), color = "#bbbbbb")+
    geom_hline(aes(yintercept = seq(0.5, (0.5 + ncol(value_dataframe)), 1)), color = "#bbbbbb")+
    xlab("")+
    ylab("")+
    theme_minimal()+
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 12, color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top",
          legend.direction = "horizontal")+
    guides(fill = guide_colorbar("correlation", title.position = "top"),
           color = guide_colorbar("-log10(P value)", title.position = "top"),
           size = "none") +
    ggtitle(label = ifelse(class_title, title, ""))
  
  ### barplot
  barplot_data <- data.frame("value" = apply(value_dataframe, 2,mean, na.rm = TRUE),
                             "platform" = names(apply(value_dataframe, 2,mean, na.rm = TRUE)),
                             row.names = names(apply(value_dataframe, 2,mean, na.rm = TRUE)))
  barplot_data <- barplot_data[order, ]
  barplot_data$platform <- factor(barplot_data$platform, levels = rev(order))
  barplot <- ggplot(data = barplot_data, aes(x = platform, y = value)) +
    geom_col(aes(fill = value), color = "gray50") +
    geom_text(aes(label = sprintf("%0.2f", value)),
              nudge_y = 0.1,
              size = 4) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.title = element_blank(),
      legend.position = "top",
      axis.text.y = element_blank()
    ) +
    coord_flip() +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, (max(barplot_data$value, na.rm = TRUE) + 0.3)),
                       breaks = seq(0,1,0.2),
                       labels = as.character(seq(0,1,0.2))) +
    scale_fill_gradient(low = "#6da7a6", high = "#c6e4f4") +
    guides(fill = guide_colorbar(paste0(criteria, " score"), title.position = "top"))
  
  
  if(criteria == "accuracy" & !class_title){
    return(p1 + p2 + barplot + plot_layout(widths = c(1,5,2)))
  }else{
    return(p2 + barplot + plot_layout(widths = c(5,2)))
  }
}


cor_plot_dendro_class <- function(data, 
                                  type = "platform",
                                  criteria = "accuracy",
                                  column_names = NULL,
                                  k = 3,
                                  branch_lwd = 0.6){
  if(!any(str_detect("id", colnames(data)))){
    stop("Please input method id")
  }
  if(!any(str_detect("Category", colnames(data)))){
    stop("Please input category")
  }
  
  method_class <- unique(data$Category)
  
  method_plot <- map(method_class, function(i){
    print(i)
    sub_data <- data %>% 
      filter(Category == i)
    if(!is.null(column_names)){
      sub_data <- sub_data %>% 
        select(all_of(c("id", "Category", column_names)))
      column_names <- column_names[!colSums(sub_data[, -c(1,2)], na.rm = T)==0]
    }
    
    plot <- cor_plot_dendro(data = sub_data, 
                            type = type,
                            criteria = criteria,
                            column_names = column_names,
                            k = k,
                            branch_lwd = branch_lwd,
                            class_title = TRUE,
                            title = i)
    plot
  }) %>% setNames(method_class)
  return(method_plot)
}



cor_plot <- function(data,
                     criteria,
                     iterate_items,
                     xlab_prefix,
                     ylab,
                     ncol){
  
  p <- platform_cor_plot <- map(iterate_items, .f = function(i){
    
    tmp <- data %>% 
      select(all_of(c(i, criteria))) %>% 
      drop_na()
    
    ggplot(data, aes(x = get(i), y = get(criteria))) +
      geom_point(size = 0.8) +
      geom_smooth(method = 'lm', linewidth = 0.8) +
      stat_cor(method = "pearson", size = 2) +
      theme_bw() +
      theme(legend.position = c(0.1, 0.8),
            axis.text = element_text(size = 4, color = "black"),
            axis.title = element_text(size = 6),
            legend.title = element_text(size = 6),
            legend.text = element_text(size = 4),
            panel.grid = element_blank(),
            title = element_text(size = 7)) +
      xlab(paste0(xlab_prefix, i, " data")) +
      ylab(ylab) +
      ggtitle(i)
  })
  return(patchwork::wrap_plots(p, ncol = ncol))
}