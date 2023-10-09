library(ggplot2)
library(tidyverse)
library(patchwork)

cor_plot <- function(data,
                     criteria,
                     iterate_items,
                     xlab_prefix,
                     ylab,
                     ncol){
  
  p <- map(iterate_items, .f = function(i){
    
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