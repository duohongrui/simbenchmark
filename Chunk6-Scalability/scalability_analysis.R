library(stringr)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

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

### Shape Constrained Addictive Model for every method
classification_shape <- function(
  data,
  method_name,
  value = "gene",
  step = "estimation",
  content = "time"
){
  
  dat <- data %>% 
    filter(method == method_name) %>% 
    mutate(
      y = case_when(
        step == "estimation" & content == "time" ~ estimation_time,
        step == "estimation" & content == "memory" ~ estimation_memory,
        step == "simulation" & content == "time" ~ simulation_time,
        step == "simulation" & content == "memory" ~ simulation_memory
      ),
      x = case_when(
        value == "gene" ~ gene_num,
        value == "cell" ~ cell_num
      ),
      z = case_when(
        value == "gene" ~ cell_num,
        value == "cell" ~ gene_num
      )
    ) %>% 
    filter(z == 1000)
  
  # using glmnet instead of lm because lm would easily overfit
  datm <- stats::model.matrix(y ~ log(x) + sqrt(x) + x + I(x^2) + I(x^3), data = dat)[, -1]
  rr.cv <- glmnet::cv.glmnet(datm, dat$y, alpha = 1)
  rr.fit <- glmnet::glmnet(datm, dat$y, alpha = 1, lambda = rr.cv$lambda.min)
  coeffs <- rr.fit$beta[, 1]
  
  max_term <-
    if (any(coeffs > .25)) {
      # if at least one coefficient is very high,
      # take the most complex one
      which(coeffs > .25) %>% tail(1) %>% names()
    } else {
      # otherwise, just take the largest coefficient
      which.max(coeffs) %>% names()
    }
  
  class_result <- c("log(x)" = "<linear",
                    "sqrt(x)" = "<linear",
                    "x" = "linear",
                    "I(x^2)" = "quadratic",
                    "I(x^3)" = ">quadratic"
  )[max_term] %>% unname()
  return(dplyr::lst(model = rr.fit,
                    classification = class_result,
                    x = dat$x,
                    y = dat$y))
}

method <- unique(scalability_data$method)

scalability_model <- purrr::map(method, .f = function(method){
  map(c("estimation", "simulation"), .f = function(step){
    if(method == "scDesign" | method == "SPsimSeq"){
      map(c("time", "memory"), .f = function(content){
        result <- classification_shape(scalability_data,
                                       method_name = method,
                                       value = "gene",
                                       step = "simulation",
                                       content = content)
      }) %>% setNames(c("time", "memory"))
    }else{
      map(c("time", "memory"), .f = function(content){
        result <- classification_shape(scalability_data,
                                       method_name = method,
                                       value = "gene",
                                       step = step,
                                       content = content)
      }) %>% setNames(c("time", "memory"))
    }
  }) %>% setNames(c("estimation", "simulation"))
}) %>% setNames(method)

scalability_model$SPsimSeq$estimation <- NULL
scalability_model$scDesign$estimation <- NULL

### long tibble
scalability_long_data <- scalability_data %>% 
  pivot_longer(5:8, names_to = "feature", values_to = "value") %>% 
  separate(feature, into = c("step", "feature"), sep = "_")

class <- unlist(map(1:nrow(scalability_long_data), .f = function(index){
  method <- scalability_long_data$method[index]
  step <- scalability_long_data$step[index]
  feature <- scalability_long_data$feature[index]
  result <- scalability_model[[method]][[step]][[feature]][["classification"]]
  if(is.null(result)){
    NA
  }else{
    result
  }
}))

scalability_long_data <- scalability_long_data %>% 
  mutate(
    classification = class
  )

class_palette <- c("<linear" = "#3d87a6",
                   "linear" = "#9bcde1",
                   ">linear" = "#d73027",
                   "quadratic" = "#d78a27",
                   ">quadratic" = "#d73027")

### estimation time -- cell
est_time_cell <- scalability_long_data %>% 
  filter(step == "estimation", feature == "time", gene_num == 1000) %>% 
  ggplot(., aes(x = cell_num, y = value/max(value, na.rm = TRUE)))+
  facet_wrap(~ method, ncol = 1, scales = "free")+
  stat_smooth(se = FALSE, mapping = aes(color = classification))+
  theme(panel.background = element_rect(fill = "black", color = "red", linewidth = 2),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")+
  scale_color_manual("", values = class_palette)
  
### estimation time -- gene
est_time_gene <- scalability_long_data %>% 
  filter(step == "estimation", feature == "time", cell_num == 1000) %>% 
  ggplot(., aes(x = gene_num, y = value/max(value, na.rm = TRUE)))+
  facet_wrap(~ method, ncol = 1, scales = "free")+
  stat_smooth(se = FALSE, mapping = aes(color = classification))+
  theme(panel.background = element_rect(fill = "black", color = "red", linewidth = 2),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")+
  scale_color_manual("", values = class_palette)

### estimation memory -- cell
est_memory_cell <- scalability_long_data %>% 
  filter(step == "estimation", feature == "memory", gene_num == 1000) %>% 
  ggplot(., aes(x = cell_num, y = value/max(value, na.rm = TRUE)))+
  facet_wrap(~ method, ncol = 1, scales = "free")+
  stat_smooth(se = FALSE, mapping = aes(color = classification))+
  theme(panel.background = element_rect(fill = "black", color = "red", linewidth = 2),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")+
  scale_color_manual("", values = class_palette)


### estimation memory -- gene
est_memory_gene <- scalability_long_data %>% 
  filter(step == "estimation", feature == "memory", cell_num == 1000) %>% 
  ggplot(., aes(x = gene_num, y = value/max(value, na.rm = TRUE)))+
  facet_wrap(~ method, ncol = 1, scales = "free")+
  stat_smooth(se = FALSE, mapping = aes(color = classification))+
  theme(panel.background = element_rect(fill = "black", color = "red", linewidth = 2),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")+
  scale_color_manual("", values = class_palette)


### simulation time -- cell
sim_time_cell <- scalability_long_data %>% 
  filter(step == "simulation", feature == "time", gene_num == 1000) %>% 
  ggplot(., aes(x = cell_num, y = value/max(value, na.rm = TRUE)))+
  facet_wrap(~ method, ncol = 1, scales = "free")+
  stat_smooth(se = FALSE, mapping = aes(color = classification))+
  theme(panel.background = element_rect(fill = "black", color = "red", linewidth = 2),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")+
  scale_color_manual("", values = class_palette)

### simulation time -- gene
sim_time_gene <- scalability_long_data %>% 
  filter(step == "simulation", feature == "time", cell_num == 1000) %>% 
  ggplot(., aes(x = gene_num, y = value/max(value, na.rm = TRUE)))+
  facet_wrap(~ method, ncol = 1, scales = "free")+
  stat_smooth(se = FALSE, mapping = aes(color = classification))+
  theme(panel.background = element_rect(fill = "black", color = "red", linewidth = 2),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom")+
  scale_color_manual("", values = class_palette)

### simulation memory -- cell
sim_memory_cell <- scalability_long_data %>% 
  filter(step == "simulation", feature == "memory", gene_num == 1000) %>% 
  ggplot(., aes(x = cell_num, y = value/max(value, na.rm = TRUE)))+
  facet_wrap(~ method, ncol = 1, scales = "free")+
  stat_smooth(se = FALSE, mapping = aes(color = classification))+
  theme(panel.background = element_rect(fill = "black", color = "red", linewidth = 2),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")+
  scale_color_manual("", values = class_palette)

### simulation memory -- gene
sim_memory_gene <- scalability_long_data %>% 
  filter(step == "simulation", feature == "memory", cell_num == 1000) %>% 
  ggplot(., aes(x = gene_num, y = value/max(value, na.rm = TRUE)))+
  facet_wrap(~ method, ncol = 1, scales = "free")+
  stat_smooth(se = FALSE, mapping = aes(color = classification))+
  theme(panel.background = element_rect(fill = "black", color = "red", linewidth = 2),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")+
  scale_color_manual("", values = class_palette)


g <- wrap_plots(
  est_time_cell,
  est_time_gene,
  est_memory_cell,
  est_memory_gene,
  sim_time_cell,
  sim_time_gene,
  sim_memory_cell,
  sim_memory_gene,
  nrow = 1
)

g


