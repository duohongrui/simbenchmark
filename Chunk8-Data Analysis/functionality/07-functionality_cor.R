library(tidyverse)
library(ggpmisc)

overall_data <- readRDS("./Chunk8-Data Analysis/overall_data.rds")

function_data <- overall_data %>% 
  select(c(1:3, 12,69:121)) %>% 
  slice(1:43)
colnames(function_data)[33:57] <- str_split(colnames(function_data)[33:57], "_", simplify = TRUE)[, 2]
platform <- colnames(function_data)[35:57]

source("./Chunk8-Data Analysis/utils.plot.R")
###--------------------------------------------------------------------------###
###             Supplementary Fig.19 (platform linear regression)
###--------------------------------------------------------------------------###

technique_cor_plot <- cor_plot(data = function_data,
                               criteria = "functionality",
                               iterate_items = platform,
                               anno_y = 0.79,
                               xlab_prefix = "Functionality scores of ",
                               ylab = "Overall functionality scores",
                               ncol = 4)

ggsave(plot = technique_cor_plot,
       filename = "/Users/duohongrui/Desktop/sim-article/figures/Supp_Fig19.pdf",
       units = "cm",
       width = 18,
       height = 19)


###--------------------------------------------------------------------------###
###             Supplementary Fig.20 (metric linear regression)
###--------------------------------------------------------------------------###
metrics <- colnames(function_data)[c(6:11, 13:19, 21:27, 29:32)]
metric_cor_plot <- cor_plot(data = function_data,
                            criteria = "functionality",
                            iterate_items = metrics,
                            anno_y = 0.79,
                            xlab_prefix = "Functionality scores of ",
                            ylab = "Overall functionality scores",
                            ncol = 4)
ggsave(plot = metric_cor_plot,
       filename = "/Users/duohongrui/Desktop/sim-article/figures/Supp_Fig20.pdf",
       units = "cm",
       width = 18,
       height = 19)

###--------------------------------------------------------------------------###
###             Supplementary Fig.21 (Cor plot for platforms)
###--------------------------------------------------------------------------###
source("./Chunk8-Data Analysis/utils.plot.R")

cor_plot_platform <- cor_plot_dendro(data = function_data, criteria = "functionality")

ggsave(plot = cor_plot_platform,
       filename = "/Users/duohongrui/Desktop/sim-article/figures/Supp_Fig21.pdf",
       units = "in",
       width = 10,
       height = 7.5)

###--------------------------------------------------------------------------###
###             Supplementary Fig.22 (Cor plot for metrics)
###--------------------------------------------------------------------------###
source("./Chunk8-Data Analysis/utils.plot.R")

metrics <- colnames(function_data)[c(6:11, 13:19, 21:27, 29:32)]
cor_plot_metric <- cor_plot_dendro(data = function_data,
                                   type = "metric",
                                   column_names = metrics,
                                   criteria = "functionality")

ggsave(plot = cor_plot_metric,
       filename = "/Users/duohongrui/Desktop/sim-article/figures/Supp_Fig22.pdf",
       units = "in",
       width = 10,
       height = 7.5)


###--------------------------------------------------------------------------###
###   Supplementary Fig.23-24 (Cor plot for platforms, metrics of each class)
###--------------------------------------------------------------------------###
source("./Chunk8-Data Analysis/utils.plot.R")

cor_plot_platform <- cor_plot_dendro_class(data = function_data)
cor_plot_metric <- cor_plot_dendro_class(data = function_data,
                                         type = "metric",
                                         column_names = colnames(function_data)[c(6:11, 13:19, 21:27, 29:32)])

ggsave(plot = wrap_plots(cor_plot_platform, ncol = 2),
       filename = "/Users/duohongrui/Desktop/sim-article/figures/Supp_Fig23.pdf",
       units = "in",
       width = 16,
       height = 16)

ggsave(plot = wrap_plots(cor_plot_metric, ncol = 2),
       filename = "/Users/duohongrui/Desktop/sim-article/figures/Supp_Fig24.pdf",
       units = "in",
       width = 16,
       height = 16)
