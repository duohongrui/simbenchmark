library(tidyverse)
library(ggdendro)
library(dendextend)
library(patchwork)
library(ggpmisc)

overall_data <- readRDS("./Chunk8-Data Analysis/overall_data.rds") %>% 
  select(1:62)
colnames(overall_data)[39:62] <- str_split(colnames(overall_data)[39:62], "_", simplify = TRUE)[, 2]
platform <- colnames(overall_data)[39:62]

source("./Chunk8-Data Analysis/utils.plot.R")
###--------------------------------------------------------------------------###
###             Supplementary Fig.4 (platform linear regression)
###--------------------------------------------------------------------------###

platform_cor_plot <- cor_plot(data = overall_data,
                              criteria = "accuracy",
                              iterate_items = platform,
                              xlab_prefix = "Accuracy scores on ",
                              ylab = "Overall accuracy scores",
                              ncol = 4)

ggsave(plot = platform_cor_plot,
       filename = "/Users/duohongrui/Desktop/sim-article/figures/Supp_Fig_4.pdf",
       units = "cm",
       width = 18,
       height = 19)

###--------------------------------------------------------------------------###
###             Supplementary Fig.7 (03-property.R)
###--------------------------------------------------------------------------###

###--------------------------------------------------------------------------###
###             Supplementary Fig.8 (03-property.R)
###--------------------------------------------------------------------------###

###--------------------------------------------------------------------------###
###     Supplementary Fig.9 (Cor plot for platforms, metrics and properties)
###--------------------------------------------------------------------------###
cor_plot_platform <- cor_plot_dendro(data = overall_data)
cor_plot_metric <- cor_plot_dendro(data = overall_data, type = "metric", column_names = colnames(overall_data)[15:22])
cor_plot_property <- cor_plot_dendro(data = overall_data, type = "property", column_names = colnames(overall_data)[23:37])
ggsave(plot = cor_plot_platform / (cor_plot_property | cor_plot_metric) + plot_annotation(tag_levels = "a") + plot_layout(heights = c(2,1), widths = c(1.5,1)),
       filename = "/Users/duohongrui/Desktop/sim-article/figures/Supp_Fig9.pdf",
       units = "in",
       width = 14,
       height = 15)

###--------------------------------------------------------------------------###
###  Supplementary Fig.10-12 (Cor plot for platforms, metrics and properties of each class)
###--------------------------------------------------------------------------###
cor_plot_platform <- cor_plot_dendro_class(data = overall_data)
cor_plot_metric <- cor_plot_dendro_class(data = overall_data, type = "metric", column_names = colnames(overall_data)[15:22])
cor_plot_property <- cor_plot_dendro_class(data = overall_data, type = "property", column_names = colnames(overall_data)[23:37])


ggsave(plot = wrap_plots(cor_plot_platform, ncol = 2),
       filename = "/Users/duohongrui/Desktop/sim-article/figures/Supp_Fig10.pdf",
       units = "in",
       width = 18,
       height = 25.3)
ggsave(plot = wrap_plots(cor_plot_metric, ncol = 2),
       filename = "/Users/duohongrui/Desktop/sim-article/figures/Supp_Fig11.pdf",
       units = "in",
       width = 10,
       height = 14)
ggsave(plot = wrap_plots(cor_plot_property, ncol = 2),
       filename = "/Users/duohongrui/Desktop/sim-article/figures/Supp_Fig12.pdf",
       units = "in",
       width = 12,
       height = 16.9)


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