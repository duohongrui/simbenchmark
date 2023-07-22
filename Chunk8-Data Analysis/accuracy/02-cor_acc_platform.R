library(tidyverse)
library(ggdendro)
library(dendextend)
library(patchwork)

overall_data <- readRDS("./Chunk8-Data Analysis/overall_data.rds") %>% 
  select(1:62)
colnames(overall_data)[38:62] <- str_split(colnames(overall_data)[38:62], "_", simplify = TRUE)[, 2]
platform <- colnames(overall_data)[38:60]

source("./Chunk8-Data Analysis/utils.plot.R")
###--------------------------------------------------------------------------###
###             Supplementary Fig.4 (platform linear regression)
###--------------------------------------------------------------------------###

platform_cor_plot <- cor_plot(data = overall_data,
                              criteria = "accuracy",
                              iterate_items = platform,
                              anno_y = 0.7,
                              xlab_prefix = "Accuracy scores of ",
                              ylab = "Overall accuracy scores",
                              ncol = 4)

ggsave(plot = platform_cor_plot,
       filename = "/Users/duohongrui/Desktop/sim-article/figures/Supp_Fig4.pdf",
       units = "cm",
       width = 18,
       height = 19)

###--------------------------------------------------------------------------###
###             Supplementary Fig.5 (03-property.R)
###--------------------------------------------------------------------------###

###--------------------------------------------------------------------------###
###             Supplementary Fig.6 (03-property.R)
###--------------------------------------------------------------------------###

###--------------------------------------------------------------------------###
###             Supplementary Fig.7 (Cor plot for platforms)
###--------------------------------------------------------------------------###

cor_plot_platform <- cor_plot_dendro(data = overall_data)

ggsave(plot = cor_plot_platform,
       filename = "/Users/duohongrui/Desktop/sim-article/figures/Supp_Fig7.pdf",
       units = "in",
       width = 10,
       height = 7.5)

###--------------------------------------------------------------------------###
###             Supplementary Fig.8 (Cor plot for metrics and properties)
###--------------------------------------------------------------------------###
cor_plot_metric <- cor_plot_dendro(data = overall_data %>% select(15:22), type = "metric")
cor_plot_property <- cor_plot_dendro(data = overall_data %>% select(23:37), type = "property")

ggsave(plot = cor_plot_metric / cor_plot_property + plot_annotation(tag_levels = "a"),
       filename = "/Users/duohongrui/Desktop/sim-article/figures/Supp_Fig8.pdf",
       units = "in",
       width = 8,
       height = 10)


# p2 <- quickcor(cor_mat,
#                type = "lower",
#                is.cor = TRUE,
#                show.diag = TRUE,
#                method = "pearson",
#                na.rm = TRUE,
#                use = "pairwise.complete.obs") +
#   geom_square()+
#   scale_fill_gradient2( high = 'orange', mid = 'white',low = 'navyblue', midpoint = 0.6) +
#   theme(
#     legend.position = "bottom"
#   )