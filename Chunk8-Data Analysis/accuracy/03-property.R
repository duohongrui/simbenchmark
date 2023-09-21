library(tidyverse)
library(ggpmisc)

overall_data <- readRDS("./Chunk8-Data Analysis/overall_data.rds") %>% 
  select(1:62)
colnames(overall_data)[38:62] <- str_split(colnames(overall_data)[38:62], "_", simplify = TRUE)[, 2]
property <- colnames(overall_data)[23:37]

source("./Chunk8-Data Analysis/utils.plot.R")
###--------------------------------------------------------------------------###
###             Supplementary Fig.7 (property linear regression)
###--------------------------------------------------------------------------###

property_cor_plot <- cor_plot(data = overall_data,
                              criteria = "accuracy",
                              iterate_items = property,
                              anno_y = 0.7,
                              xlab_prefix = "Accuracy scores of ",
                              ylab = "Overall accuracy scores",
                              ncol = 4)

ggsave(plot = property_cor_plot,
       filename = "/Users/duohongrui/Desktop/sim-article/figures/Supp_Fig7.pdf",
       units = "cm",
       width = 18,
       height = 15)

###--------------------------------------------------------------------------###
###             Supplementary Fig.8 (metrics linear regression)
###--------------------------------------------------------------------------###
metric <- colnames(overall_data)[15:22]
metric_cor_plot <- cor_plot(data = overall_data,
                            criteria = "accuracy",
                            iterate_items = metric,
                            anno_y = 0.7,
                            xlab_prefix = "Accuracy scores of ",
                            ylab = "Overall accuracy scores",
                            ncol = 4)

ggsave(plot = metric_cor_plot,
       filename = "/Users/duohongrui/Desktop/sim-article/figures/Supp_Fig8.pdf",
       units = "cm",
       width = 18,
       height = 9)
