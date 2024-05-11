library(tidyverse)
library(patchwork)
library(ggpmisc)
library(ggpubr)

overall_data <- readRDS("./Chunk8-Data Analysis/overall_data.rds") %>% 
  select(1:62)
colnames(overall_data)[39:62] <- str_split(colnames(overall_data)[39:62], "_", simplify = TRUE)[, 2]
platform <- colnames(overall_data)[39:62]

source("./Chunk8-Data Analysis/2-accuracy/05-utils.plot.R")
###--------------------------------------------------------------------------###
###             Supplementary Fig.5 (platform linear regression)
###--------------------------------------------------------------------------###

platform_cor_plot <- cor_plot(data = overall_data,
                              criteria = "accuracy",
                              iterate_items = platform,
                              xlab_prefix = "Accuracy scores on ",
                              ylab = "Overall accuracy scores",
                              ncol = 4)

ggsave(plot = platform_cor_plot,
       filename = "/Users/duohongrui/Desktop/sim-article/figures/Supp_Fig_5.pdf",
       units = "cm",
       width = 18,
       height = 19)