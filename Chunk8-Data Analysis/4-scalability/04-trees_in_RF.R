source("./Chunk8-Data Analysis/4-scalability/05-utils_functions.R")

data_list <- list.files("../scalability/")

### read data into a list
all_result <- list()
for (i in data_list) {
  result <- readRDS(file.path("../scalability", i))
  all_result[[i]] <- result
}

### turn to a tibble
scalability_data <- purrr::map_dfr(1:length(all_result), .f = function(index){
  all_result[[index]]
})
data <- readRDS("./Chunk8-Data Analysis/4-scalability/execution_detection_data.rds")
colnames(data) <- c("method",
                    "data_id",
                    "cell_num",
                    "gene_num",
                    "estimation_time",
                    "estimation_memory",
                    "simulation_time",
                    "simulation_memory")
data[c(1,2), c(6,8)] <- NA
data <- data %>% 
  select(-2) %>% 
  bind_rows(., scalability_data %>% select(-4)) %>% 
  mutate(
    row = row_number()
  )
data <- data %>% 
  filter(estimation_memory != 0 | is.na(estimation_memory))
data[data == 0] <- 0.01
### scGAN NA for memory
data[data$method == "scGAN", "estimation_memory"] <- NA
data[data$method == "scGAN", "simulation_memory"] <- NA

ntrees <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 800, 1000, 3000, 4000, 5000)
methods <- unique(data$method)
RF_gradient_result <- RF_gradient_function(ntree = ntrees, methods = methods)
saveRDS(RF_gradient_result, "./Chunk8-Data Analysis/4-scalability/RF_gradient_result.rds")

RF_gradient_result$tree <- factor(RF_gradient_result$tree,
                                  levels = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 800, 1000, 3000, 4000, 5000))



#### Plot
line_data <- RF_gradient_result %>% 
  pivot_longer(cols = colnames(.)[5:55], names_to = "method", values_to = "correlation") %>% 
  group_by(step, feature, tree, method) %>% 
  summarise(
    mean_correlation = mean(correlation, na.rm = TRUE)
  )
line_data$tree <- as.numeric(line_data$tree)

################################################################################
###################    Supplementary Figures 21-24   ###########################
################################################################################
### estimation-time
est_time_plots <- RF_tree_plot_function(Step = "estimation",
                                        Feature = "time",
                                        Methods = methods,
                                        line_data = line_data,
                                        RF_gradient_result = RF_gradient_result)
ggsave(plot = wrap_plots(est_time_plots, ncol = 5),
       width = 18,
       height = 22,
       units = "cm",
       filename = "../sim-article/figures/Supp_Fig_21.pdf")

### estimation-memory
est_memory_plots <- RF_tree_plot_function(Step = "estimation",
                                          Feature = "memory",
                                          Methods = methods,
                                          line_data = line_data,
                                          RF_gradient_result = RF_gradient_result)
ggsave(plot = wrap_plots(est_memory_plots, ncol = 5),
       width = 18,
       height = 22,
       units = "cm",
       filename = "../sim-article/figures/Supp_Fig_22.pdf")

### simulation-time
sim_time_plots <- RF_tree_plot_function(Step = "simulation",
                                        Feature = "time",
                                        Methods = methods,
                                        line_data = line_data,
                                        RF_gradient_result = RF_gradient_result)
ggsave(plot = wrap_plots(sim_time_plots, ncol = 5),
       width = 18,
       height = 22,
       units = "cm",
       filename = "../sim-article/figures/Supp_Fig_23.pdf")


### simulation-memory
sim_memory_plots <- RF_tree_plot_function(Step = "simulation",
                                          Feature = "memory",
                                          Methods = methods,
                                          line_data = line_data,
                                          RF_gradient_result = RF_gradient_result)
ggsave(plot = wrap_plots(sim_memory_plots, ncol = 5),
       width = 18,
       height = 22,
       units = "cm",
       filename = "../sim-article/figures/Supp_Fig_24.pdf")
