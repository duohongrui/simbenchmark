library(ggplot2)
library(tidyverse)
library(ggdendro)
library(dendextend)
library(patchwork)

cor_plot_dendro <- function(data, 
                            type = "platform",
                            criteria = "accuracy",
                            k = 3,
                            branch_lwd = 0.6){
  
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
      select(all_of(items))
  }else{
    data <- data
  }
  
  ### correlation
  cor_res <- Hmisc::rcorr(data %>% as.matrix(), type = "pearson")
  cor_mat <- cor_res$r
  p_mat <- cor_res$P
  
  if(criteria == "accuracy"){
    ### hclust
    clust <- hclust(dist(t(data)), method = "ward.D", )
    
    order <- rev(clust$labels[clust$order])
    
    ### correlation
    cor_mat <- cor_res$r[order, order]
    p_mat <- cor_res$P[order, order]
    
    ggd <- as.ggdend(clust %>%
                       as.dendrogram() %>%
                       set("branches_k_color", k = k) %>% 
                       set("branches_lwd", branch_lwd))
    
    p1 <- ggplot(ggd, horiz = TRUE, theme = NULL, labels = FALSE) +
      theme_dendro()
  }else{
    order <- colnames(data)
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
    geom_vline(aes(xintercept = seq(0.5, (0.5 + ncol(data)), 1)), color = "#bbbbbb")+
    geom_hline(aes(yintercept = seq(0.5, (0.5 + ncol(data)), 1)), color = "#bbbbbb")+
    xlab("")+
    ylab("")+
    theme_minimal()+
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right")+
    guides(fill = guide_colorbar("correlation", title.position = "top"),
           color = guide_colorbar("-log10(P value)", title.position = "top"),
           size = "none")
  
  if(criteria == "accuracy"){
    return(p1 + p2 + plot_layout(widths = c(1,5)))
  }else{
    return(p2)
  }
}



cor_plot <- function(data,
                     criteria,
                     iterate_items,
                     anno_y,
                     xlab_prefix,
                     ylab,
                     ncol){
  
  p <- platform_cor_plot <- map(iterate_items, .f = function(i){
    
    tmp <- data %>% 
      select(all_of(c(i, criteria))) %>% 
      drop_na()
    
    cor_value <- sprintf("%0.2f", cor(tmp %>% pull(i),
                                      tmp %>% pull(get(criteria))))
    
    ggplot(data, aes(x = get(i), y = get(criteria))) +
      geom_point(size = 0.8) +
      geom_smooth(method = 'lm', linewidth = 0.8) +
      stat_poly_eq(use_label(c("eq", "adj.R2", "p.value.label"), sep = "*\", \"*"),
                   formula = y ~ x,
                   parse = TRUE,
                   size = 1.5,
                   label.x = 0.98,
                   label.y = 0.02) +
      annotate(geom = "text",
               x = mean(tmp %>% pull(i), na.rm = TRUE),
               y = anno_y,
               label = paste0("Cor: ", cor_value),
               size = 2) +
      theme_bw() +
      theme(legend.position = c(0.1, 0.8),
            axis.text = element_text(size = 4),
            axis.title = element_text(size = 6),
            legend.title = element_text(size = 6),
            legend.text = element_text(size = 4),
            panel.grid = element_blank(),
            title = element_text(size = 7)) +
      xlab(paste0(xlab_prefix, i)) +
      ylab(ylab) +
      ggtitle(i)
  })
  return(patchwork::wrap_plots(p, ncol = ncol))
}