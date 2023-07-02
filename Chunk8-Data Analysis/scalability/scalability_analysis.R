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

### scalability score

#### aggregate by repeat time
score_data <- scalability_data %>% 
  group_by(method, cell_num, gene_num) %>% 
  summarise(
    estimation_time = mean(estimation_time, na.rm = TRUE),
    estimation_memory = mean(estimation_memory, na.rm = TRUE),
    simulation_time = mean(simulation_time, na.rm = TRUE),
    simulation_memory = mean(simulation_memory, na.rm = TRUE)
  ) %>% 
  ungroup()
  
#### aggregate by dataset
score_data <- score_data %>% 
  group_by(cell_num, gene_num) %>% 
  mutate(
    across(all_of(c("estimation_time",
                    "estimation_memory",
                    "simulation_time",
                    "simulation_memory")), ~ pnorm((.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
  ) %>% 
  ungroup()

#### aggregate by method
score_data <- score_data %>% 
  group_by(method) %>% 
  summarise(
    estimation_time = mean(estimation_time, na.rm = TRUE),
    estimation_memory = mean(estimation_memory, na.rm = TRUE),
    simulation_time = mean(simulation_time, na.rm = TRUE),
    simulation_memory = mean(simulation_memory, na.rm = TRUE)
  ) %>% 
  mutate(
    across(all_of(c("estimation_time",
                    "estimation_memory",
                    "simulation_time",
                    "simulation_memory")), ~ 1 - .x)
  ) %>% 
  ungroup()
score_data[42, c(3,5)] <- NaN
### time and memory score
score_data <- score_data %>% 
  mutate(
    time_score = apply(.[, -1], MARGIN = 1, function(x){mean(c(x[1], x[3]), na.rm = TRUE)}),
    memory_score = apply(.[, -1], MARGIN = 1, function(x){mean(c(x[2], x[4]), na.rm = TRUE)}),
    scalability_score = apply(.[, -1], MARGIN = 1, function(x){mean(c(x[1], x[2], x[3], x[4]), na.rm = TRUE)})
  )
saveRDS(score_data, file = "Chunk8-Data Analysis/scalability/score_data.rds")


################################################################################
#######################                Model           #########################
################################################################################
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

scalability_gene_model <- purrr::map(method, .f = function(method){
  map(c("estimation", "simulation"), .f = function(step){
    if(method == "scDesign" | method == "SPsimSeq" | method == "SimBPDD"){
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

scalability_cell_model <- purrr::map(method, .f = function(method){
  map(c("estimation", "simulation"), .f = function(step){
    if(method == "scDesign" | method == "SPsimSeq" | method == "SimBPDD"){
      map(c("time", "memory"), .f = function(content){
        result <- classification_shape(scalability_data,
                                       method_name = method,
                                       value = "cell",
                                       step = "simulation",
                                       content = content)
      }) %>% setNames(c("time", "memory"))
    }else{
      map(c("time", "memory"), .f = function(content){
        result <- classification_shape(scalability_data,
                                       method_name = method,
                                       value = "cell",
                                       step = step,
                                       content = content)
      }) %>% setNames(c("time", "memory"))
    }
  }) %>% setNames(c("estimation", "simulation"))
}) %>% setNames(method)
scalability_gene_model$SPsimSeq$estimation <- NULL
scalability_gene_model$scDesign$estimation <- NULL
scalability_gene_model$SimBPDD$estimation <- NULL
scalability_cell_model$SPsimSeq$estimation <- NULL
scalability_cell_model$scDesign$estimation <- NULL
scalability_cell_model$SimBPDD$estimation <- NULL
saveRDS(scalability_gene_model, file = "Chunk8-Data Analysis/scalability/scalability_gene_model.rds")
saveRDS(scalability_cell_model, file = "Chunk8-Data Analysis/scalability/scalability_cell_model.rds")



### long tibble
# scalability_long_data <- tibble(
#   method = rep(unique(scalability_data$method), each = 4),
#   step = rep(rep(c("estimation", "simulation"), 2), 45),
#   feature = rep(rep(c("time", "memory"), each = 2), 45)
# )
scalability_long_data <- scalability_data %>% 
  pivot_longer(5:8, names_to = "feature", values_to = "value") %>% 
  separate(feature, into = c("step", "feature"), sep = "_")

gene_trend_class <- unlist(map(1:nrow(scalability_long_data), .f = function(index){
  method <- scalability_long_data$method[index]
  step <- scalability_long_data$step[index]
  feature <- scalability_long_data$feature[index]
  result <- scalability_gene_model[[method]][[step]][[feature]][["classification"]]
  if(is.null(result)){
    NA
  }else{
    result
  }
}))

gene_coef <- map(1:nrow(scalability_long_data), .f = function(index){
  method <- scalability_long_data$method[index]
  step <- scalability_long_data$step[index]
  feature <- scalability_long_data$feature[index]
  result <- coefficients(scalability_gene_model[[method]][[step]][[feature]][["model"]])[, 1]
  if(is.null(result)){
    NA
  }else{
    result
  }
})

cell_trend_class <- unlist(map(1:nrow(scalability_long_data), .f = function(index){
  method <- scalability_long_data$method[index]
  step <- scalability_long_data$step[index]
  feature <- scalability_long_data$feature[index]
  result <- scalability_cell_model[[method]][[step]][[feature]][["classification"]]
  if(is.null(result)){
    NA
  }else{
    result
  }
}))

cell_coef <- map(1:nrow(scalability_long_data), .f = function(index){
  method <- scalability_long_data$method[index]
  step <- scalability_long_data$step[index]
  feature <- scalability_long_data$feature[index]
  result <- coefficients(scalability_cell_model[[method]][[step]][[feature]][["model"]])[, 1]
  if(is.null(result)){
    NA
  }else{
    result
  }
})

scalability_long_data <- scalability_long_data %>% 
  mutate(
    gene_trend_class = gene_trend_class,
    cell_trend_class = cell_trend_class,
    gene_coef = gene_coef,
    cell_coef = cell_coef
  )
saveRDS(scalability_long_data, file = "Chunk8-Data Analysis/scalability/scalability_long_data.rds")

##########################################################
###################   SUMMARY FIGURE   ###################
##########################################################

class_palette <- c("<linear" = "#3d87a6",
                   "linear" = "#4DAF4A",
                   ">linear" = "#984EA3",
                   "quadratic" = "#FF7F00",
                   ">quadratic" = "#E41A1C")

### estimation time -- cell
est_time_cell <- map(unique(scalability_long_data$method), function(y){
  tmp <- scalability_long_data %>% 
    filter(step == "estimation", feature == "time", gene_num == 1000, method == y)
  if(y == "SCRIP-GP-commonBCV" |
     y == "SCRIP-GP-trendedBCV" |
     y == "SCRIP-BGP-commonBCV" |
     y == "SCRIP-BGP-trendedBCV"){
    y <- paste0(paste0(str_split(y, pattern = "-", simplify = TRUE)[1, ][1:2], collapse = "-"),
                "-\n",
                paste0(str_split(y, pattern = "-", simplify = TRUE)[1, ][3], collapse = "-"))
  }
  p <- tmp %>% 
    ggplot(aes(x = cell_num, y = value))+
    geom_point(color = "white", size = 0.7)+
    stat_function(fun = function(x){
      tmp2 <- tmp %>% 
        pull(cell_coef)
      tmp2 <- tmp2[[1]]
      tmp2[1] + tmp2[2]*log(x) + tmp2[3]*sqrt(x) + tmp2[4]*x + tmp2[5]*x^2 + tmp2[6]*x^3
    }, aes(color = cell_trend_class)) +
    theme(panel.background = element_rect(fill = "black", linewidth = 2),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.position = "none",
          title = element_text(size = 8))+
    scale_color_manual("", values = class_palette)+
    ggtitle(label = y, subtitle = tmp$cell_trend_class[1])
  p
}) %>% setNames(unique(scalability_long_data$method))
ggsave(plot = wrap_plots(est_time_cell, ncol = 9),
       filename = "./Chunk8-Data Analysis/scalability/est_time_cell.pdf",
       width = 10,
       height = 6,
       units = "in")


### estimation time -- gene
est_time_gene <- map(unique(scalability_long_data$method), function(y){
  tmp <- scalability_long_data %>% 
    filter(step == "estimation", feature == "time", cell_num == 1000, method == y)
  if(y == "SCRIP-GP-commonBCV" |
     y == "SCRIP-GP-trendedBCV" |
     y == "SCRIP-BGP-commonBCV" |
     y == "SCRIP-BGP-trendedBCV"){
    y <- paste0(paste0(str_split(y, pattern = "-", simplify = TRUE)[1, ][1:2], collapse = "-"),
                "-\n",
                paste0(str_split(y, pattern = "-", simplify = TRUE)[1, ][3], collapse = "-"))
  }
  p <- tmp %>% 
    ggplot(aes(x = gene_num, y = value))+
    geom_point(color = "white", size = 0.7)+
    stat_function(fun = function(x){
      tmp2 <- tmp %>% 
        pull(cell_coef)
      tmp2 <- tmp2[[1]]
      tmp2[1] + tmp2[2]*log(x) + tmp2[3]*sqrt(x) + tmp2[4]*x + tmp2[5]*x^2 + tmp2[6]*x^3
    }, aes(color = cell_trend_class)) +
    theme(panel.background = element_rect(fill = "black", linewidth = 2),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.position = "none",
          title = element_text(size = 8))+
    scale_color_manual("", values = class_palette)+
    ggtitle(label = y, subtitle = tmp$cell_trend_class[1])
  p
}) %>% setNames(unique(scalability_long_data$method))
ggsave(plot = wrap_plots(est_time_gene, ncol = 9),
       filename = "./Chunk8-Data Analysis/scalability/est_time_gene.pdf",
       width = 10,
       height = 6,
       units = "in")


### estimation memory -- cell
est_memory_cell <- map(unique(scalability_long_data$method), function(y){
  tmp <- scalability_long_data %>% 
    filter(step == "estimation", feature == "memory", gene_num == 1000, method == y)
  if(y == "SCRIP-GP-commonBCV" |
     y == "SCRIP-GP-trendedBCV" |
     y == "SCRIP-BGP-commonBCV" |
     y == "SCRIP-BGP-trendedBCV"){
    y <- paste0(paste0(str_split(y, pattern = "-", simplify = TRUE)[1, ][1:2], collapse = "-"),
                "-\n",
                paste0(str_split(y, pattern = "-", simplify = TRUE)[1, ][3], collapse = "-"))
  }
  p <- tmp %>% 
    ggplot(aes(x = cell_num, y = value))+
    geom_point(color = "white", size = 0.7)+
    stat_function(fun = function(x){
      tmp2 <- tmp %>% 
        pull(cell_coef)
      tmp2 <- tmp2[[1]]
      tmp2[1] + tmp2[2]*log(x) + tmp2[3]*sqrt(x) + tmp2[4]*x + tmp2[5]*x^2 + tmp2[6]*x^3
    }, aes(color = cell_trend_class)) +
    theme(panel.background = element_rect(fill = "black", linewidth = 2),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.position = "none",
          title = element_text(size = 8))+
    scale_color_manual("", values = class_palette)+
    ggtitle(label = y, subtitle = tmp$cell_trend_class[1])
  p
}) %>% setNames(unique(scalability_long_data$method))
ggsave(plot = wrap_plots(est_memory_cell, ncol = 9),
       filename = "./Chunk8-Data Analysis/scalability/est_memory_cell.pdf",
       width = 10,
       height = 6,
       units = "in")


### estimation memory -- gene
est_memory_gene <- map(unique(scalability_long_data$method), function(y){
  tmp <- scalability_long_data %>% 
    filter(step == "estimation", feature == "memory", cell_num == 1000, method == y)
  if(y == "SCRIP-GP-commonBCV" |
     y == "SCRIP-GP-trendedBCV" |
     y == "SCRIP-BGP-commonBCV" |
     y == "SCRIP-BGP-trendedBCV"){
    y <- paste0(paste0(str_split(y, pattern = "-", simplify = TRUE)[1, ][1:2], collapse = "-"),
                "-\n",
                paste0(str_split(y, pattern = "-", simplify = TRUE)[1, ][3], collapse = "-"))
  }
  p <- tmp %>% 
    ggplot(aes(x = gene_num, y = value))+
    geom_point(color = "white", size = 0.7)+
    stat_function(fun = function(x){
      tmp2 <- tmp %>% 
        pull(cell_coef)
      tmp2 <- tmp2[[1]]
      tmp2[1] + tmp2[2]*log(x) + tmp2[3]*sqrt(x) + tmp2[4]*x + tmp2[5]*x^2 + tmp2[6]*x^3
    }, aes(color = cell_trend_class)) +
    theme(panel.background = element_rect(fill = "black", linewidth = 2),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.position = "none",
          title = element_text(size = 8))+
    scale_color_manual("", values = class_palette)+
    ggtitle(label = y, subtitle = tmp$cell_trend_class[1])
  p
}) %>% setNames(unique(scalability_long_data$method))
ggsave(plot = wrap_plots(est_memory_gene, ncol = 9),
       filename = "./Chunk8-Data Analysis/scalability/est_memory_gene.pdf",
       width = 10,
       height = 6,
       units = "in")



### simulation time -- cell
sim_time_cell <- map(unique(scalability_long_data$method), function(y){
  tmp <- scalability_long_data %>% 
    filter(step == "simulation", feature == "time", gene_num == 1000, method == y)
  if(y == "SCRIP-GP-commonBCV" |
     y == "SCRIP-GP-trendedBCV" |
     y == "SCRIP-BGP-commonBCV" |
     y == "SCRIP-BGP-trendedBCV"){
    y <- paste0(paste0(str_split(y, pattern = "-", simplify = TRUE)[1, ][1:2], collapse = "-"),
                "-\n",
                paste0(str_split(y, pattern = "-", simplify = TRUE)[1, ][3], collapse = "-"))
  }
  p <- tmp %>% 
    ggplot(aes(x = cell_num, y = value))+
    geom_point(color = "white", size = 0.7)+
    stat_function(fun = function(x){
      tmp2 <- tmp %>% 
        pull(cell_coef)
      tmp2 <- tmp2[[1]]
      tmp2[1] + tmp2[2]*log(x) + tmp2[3]*sqrt(x) + tmp2[4]*x + tmp2[5]*x^2 + tmp2[6]*x^3
    }, aes(color = cell_trend_class)) +
    theme(panel.background = element_rect(fill = "black", linewidth = 2),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.position = "none",
          title = element_text(size = 8))+
    scale_color_manual("", values = class_palette)+
    ggtitle(label = y, subtitle = tmp$cell_trend_class[1])
  p
}) %>% setNames(unique(scalability_long_data$method))
ggsave(plot = wrap_plots(sim_time_cell, ncol = 9),
       filename = "./Chunk8-Data Analysis/scalability/sim_time_cell.pdf",
       width = 10,
       height = 6,
       units = "in")



### simulation time -- gene
sim_time_gene <- map(unique(scalability_long_data$method), function(y){
  tmp <- scalability_long_data %>% 
    filter(step == "simulation", feature == "time", cell_num == 1000, method == y)
  if(y == "SCRIP-GP-commonBCV" |
     y == "SCRIP-GP-trendedBCV" |
     y == "SCRIP-BGP-commonBCV" |
     y == "SCRIP-BGP-trendedBCV"){
    y <- paste0(paste0(str_split(y, pattern = "-", simplify = TRUE)[1, ][1:2], collapse = "-"),
                "-\n",
                paste0(str_split(y, pattern = "-", simplify = TRUE)[1, ][3], collapse = "-"))
  }
  p <- tmp %>% 
    ggplot(aes(x = gene_num, y = value))+
    geom_point(color = "white", size = 0.7)+
    stat_function(fun = function(x){
      tmp2 <- tmp %>% 
        pull(cell_coef)
      tmp2 <- tmp2[[1]]
      tmp2[1] + tmp2[2]*log(x) + tmp2[3]*sqrt(x) + tmp2[4]*x + tmp2[5]*x^2 + tmp2[6]*x^3
    }, aes(color = cell_trend_class)) +
    theme(panel.background = element_rect(fill = "black", linewidth = 2),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.position = "none",
          title = element_text(size = 8))+
    scale_color_manual("", values = class_palette)+
    ggtitle(label = y, subtitle = tmp$cell_trend_class[1])
  p
}) %>% setNames(unique(scalability_long_data$method))
ggsave(plot = wrap_plots(sim_time_gene, ncol = 9),
       filename = "./Chunk8-Data Analysis/scalability/sim_time_gene.pdf",
       width = 10,
       height = 6,
       units = "in")


### simulation memory -- cell
sim_memory_cell <- map(unique(scalability_long_data$method), function(y){
  tmp <- scalability_long_data %>% 
    filter(step == "simulation", feature == "memory", gene_num == 1000, method == y)
  if(y == "SCRIP-GP-commonBCV" |
     y == "SCRIP-GP-trendedBCV" |
     y == "SCRIP-BGP-commonBCV" |
     y == "SCRIP-BGP-trendedBCV"){
    y <- paste0(paste0(str_split(y, pattern = "-", simplify = TRUE)[1, ][1:2], collapse = "-"),
                "-\n",
                paste0(str_split(y, pattern = "-", simplify = TRUE)[1, ][3], collapse = "-"))
  }
  p <- tmp %>% 
    ggplot(aes(x = cell_num, y = value))+
    geom_point(color = "white", size = 0.7)+
    stat_function(fun = function(x){
      tmp2 <- tmp %>% 
        pull(cell_coef)
      tmp2 <- tmp2[[1]]
      tmp2[1] + tmp2[2]*log(x) + tmp2[3]*sqrt(x) + tmp2[4]*x + tmp2[5]*x^2 + tmp2[6]*x^3
    }, aes(color = cell_trend_class)) +
    theme(panel.background = element_rect(fill = "black", linewidth = 2),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.position = "none",
          title = element_text(size = 8))+
    scale_color_manual("", values = class_palette)+
    ggtitle(label = y, subtitle = tmp$cell_trend_class[1])
  p
}) %>% setNames(unique(scalability_long_data$method))
ggsave(plot = wrap_plots(sim_memory_cell, ncol = 9),
       filename = "./Chunk8-Data Analysis/scalability/sim_memory_cell.pdf",
       width = 10,
       height = 6,
       units = "in")


### simulation memory -- gene
sim_memory_gene <- map(unique(scalability_long_data$method), function(y){
  tmp <- scalability_long_data %>% 
    filter(step == "simulation", feature == "memory", cell_num == 1000, method == y)
  if(y == "SCRIP-GP-commonBCV" |
     y == "SCRIP-GP-trendedBCV" |
     y == "SCRIP-BGP-commonBCV" |
     y == "SCRIP-BGP-trendedBCV"){
    y <- paste0(paste0(str_split(y, pattern = "-", simplify = TRUE)[1, ][1:2], collapse = "-"),
                "-\n",
                paste0(str_split(y, pattern = "-", simplify = TRUE)[1, ][3], collapse = "-"))
  }
  p <- tmp %>% 
    ggplot(aes(x = gene_num, y = value))+
    geom_point(color = "white", size = 0.7)+
    stat_function(fun = function(x){
      tmp2 <- tmp %>% 
        pull(cell_coef)
      tmp2 <- tmp2[[1]]
      tmp2[1] + tmp2[2]*log(x) + tmp2[3]*sqrt(x) + tmp2[4]*x + tmp2[5]*x^2 + tmp2[6]*x^3
    }, aes(color = cell_trend_class)) +
    theme(panel.background = element_rect(fill = "black", linewidth = 2),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.position = "none",
          title = element_text(size = 8))+
    scale_color_manual("", values = class_palette)+
    ggtitle(label = y, subtitle = tmp$cell_trend_class[1])
  p
}) %>% setNames(unique(scalability_long_data$method))
ggsave(plot = wrap_plots(sim_memory_gene, ncol = 9),
       filename = "./Chunk8-Data Analysis/scalability/sim_memory_gene.pdf",
       width = 10,
       height = 6,
       units = "in")

# g <- wrap_plots(
#   est_time_cell,
#   est_time_gene,
#   est_memory_cell,
#   est_memory_gene,
#   sim_time_cell,
#   sim_time_gene,
#   sim_memory_cell,
#   sim_memory_gene,
#   nrow = 1
# )
# 
# ggsave(plot = g, filename = "../g.plot.pdf", width = 20, height = 30, units = "in")


################################################################################
###################    Shape constrained additive models   #####################
################################################################################
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

# #### aggregate by repeat time
# score_data <- scalability_data %>% 
#   group_by(method, cell_num, gene_num) %>% 
#   summarise(
#     estimation_time = mean(estimation_time, na.rm = TRUE),
#     estimation_memory = mean(estimation_memory, na.rm = TRUE),
#     simulation_time = mean(simulation_time, na.rm = TRUE),
#     simulation_memory = mean(simulation_memory, na.rm = TRUE)
#   ) %>% 
#   ungroup()

### execution data
scam_data_list <- list.files("/Volumes/Elements/sim_bench/simulation_data/")
scam_data <- map_dfr(1:length(scam_data_list), .f = function(x){
  print(x)
  data_tmp <- readRDS(paste0("/Volumes/Elements/sim_bench/simulation_data/", scam_data_list[x]))
  method_name_split <- str_split(scam_data_list[x], "_", simplify = TRUE)
  method_name <- method_name_split[1]
  if(method_name_split[2] == "ERCC"){
    data_id <- method_name_split[3]
  }else{
    data_id <- method_name_split[2]
  }
  data <- c(method_name,
            data_id,
            data_tmp$sim_data_info$cell_num,
            data_tmp$sim_data_info$gene_num,
            data_tmp$sim_data_info$estimate_time,
            data_tmp$sim_data_info$estimate_memory,
            data_tmp$sim_data_info$simulation_time,
            data_tmp$sim_data_info$simulation_memory)
  names(data) <- c("method",
                   "data_id",
                   "cell_num",
                   "gene_num",
                   "estimate_time",
                   "estimate_memory",
                   "simulation_time",
                   "simulation_memory")
  data
})
scam_data$cell_num <- as.numeric(scam_data$cell_num)
scam_data$gene_num <- as.numeric(scam_data$gene_num)
scam_data$estimate_time <- as.numeric(scam_data$estimate_time)
scam_data$estimate_memory <- as.numeric(scam_data$estimate_memory)
scam_data$simulation_time <- as.numeric(scam_data$simulation_time)
scam_data$simulation_memory <- as.numeric(scam_data$simulation_memory)
saveRDS(scam_data, file = "./Chunk8-Data Analysis/scalability/execution_detection_data.rds")

colnames(scam_data) <- c("method",
                         "data_id",
                         "cell_num",
                         "gene_num",
                         "estimation_time",
                         "estimation_memory",
                         "simulation_time",
                         "simulation_memory")
scam_data[c(1,2), c(6,8)] <- NA
scam_data <- scam_data %>% 
  select(-2) %>% 
  bind_rows(., scalability_data %>% select(-4)) %>% 
  mutate(
    row = row_number()
  )
scam_data <- scam_data %>% 
  filter(estimation_memory != 0 | is.na(estimation_memory))
scam_data[scam_data == 0] <- 0.01

# scam_data <- scam_data %>% 
#   select(-2) %>% 
#   bind_rows(., score_data) %>% 
#   mutate(
#     row = row_number()
#   )
### scGAN NA for memory
scam_data[scam_data$method == "scGAN", "estimation_memory"] <- NA
scam_data[scam_data$method == "scGAN", "simulation_memory"] <- NA


### model function for every method
scam_function <- function(step, feature, model_selected = "scam"){
  methods <- unique(train_data$method)
  id <- paste0(step, "_", feature)
  predict_result <- map(methods, .f = function(x){
    if(x == "SPsimSeq" & step == "estimation" |
       x == "scDesign" & step == "estimation" | 
       x == "SimBPDD" & step == "estimation" | 
       x == "scGAN" & feature == "memory"){
      message(x)
      cat(paste0("Passing modeling the ", feature, " in ", step, " for ", x, "\n"))
      lst(model = NULL,
          predict_result = NULL)
    }else{
      message(x)
      ### training
      tmp <- train_data %>% 
        filter(method == x)
      ### BASiCS and TedSim NA
      if(x == "BASiCS" & feature == "memory" | x == "TedSim"){
        tmp <- tmp %>% 
          drop_na()
      }
      if(model_selected == "scam"){
        model <- scam::scam(get(id) ~ s(cell_num, gene_num, bs = "tedmi"), data = tmp) %>%
          strip::strip(keep = "predict")
      }
      if(model_selected == "RF"){
        model <- randomForest::randomForest(get(id) ~ cell_num + gene_num, data = tmp, ntree = 600)
      }
      
      ### prediction
      tmp2 <- test_data %>% 
        filter(method == x)
      if(x == "BASiCS" & feature == "memory" | x == "TedSim"){
        tmp2 <- test_data %>% 
          filter(method == x) %>% 
          drop_na()
      }
      predicted_result <- predict(model, tmp2 %>% select(2,3))
      
      predict_result <- tibble("method" = x,
                               "real" = tmp2 %>% pull(id),
                               "predicted" = predicted_result)
      colnames(predict_result) <- c("method",
                                    paste0("real_", id),
                                    paste0("predicted_", id))
      lst(model,
          predict_result)
    }
  }) %>% setNames(unique(train_data$method))
  return(predict_result)
}

bind_prediction_result <- function(model){
  tmp_data <- map_dfr(1:length(model), .f = function(x){
    model[[x]][["predict_result"]]
  })
  col_name <- colnames(tmp_data)
  colnames(tmp_data) <- c("method", "real_data", "predict_data")
  tmp_data <- tmp_data %>% 
    filter(predict_data > 0)
  colnames(tmp_data) <- col_name
  tmp_data
}


correlation_plot <- function(table){
  #################### all results
  ### step and feature
  step <- str_split(colnames(table)[2], "_", simplify = TRUE)[2]
  feature <- str_split(colnames(table)[2], "_", simplify = TRUE)[3]
  
  colnames(table) <- c("method", "x", "y")
  table$x <- log10(table$x)
  table$y <- log10(table$y)
  
  ### correlation
  cor_value <- round(cor(table$x, table$y), digits = 2)
  
  if(feature == "time"){
    feature_id <- "execution time"
  }else{
    feature_id <- "memory usage"
  }
  
  p <- ggplot(table, aes(x = x, y = y)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    theme_bw() +
    xlab(paste0("Actual log10 ", feature_id, " in ", step)) +
    ylab(paste0("Predicted log10 ", feature_id, " in ", step)) +
    ggtitle(paste0("Correlation: ", cor_value))
  
  #################### for every method
  
  method_prediction_result <- map(unique(table$method), .f = function(z){
    ### method data
    tmp <- table %>% 
      filter(method == z)
    
    ### correlation
    cor_value_method <- round(cor(tmp$x, tmp$y), digits = 2)
    
    ### Plot
    p_method <- ggplot(tmp, aes(x = x, y = y)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1, color = "red") +
      theme_bw() +
      xlab(paste0("Actual log10 ", feature_id, " in ", step)) +
      ylab(paste0("Predicted log10 ", feature_id, " in ", step)) +
      ggtitle(paste0("Correlation: ", cor_value_method), subtitle = paste0("Method: ", z))
    
    list(cor_value = cor_value_method,
         plot = p_method)
    
  }) %>% setNames(unique(table$method))
  
  ### return all results
  list(correlation = cor_value,
       correlation_plot = p,
       method_prediction_result = method_prediction_result)
}


set.seed(65)
train_data <- scam_data %>% 
  group_by(method) %>% 
  sample_frac(0.8) %>% 
  ungroup()
test_data <- scam_data[scam_data$row %in% setdiff(1:nrow(scam_data), train_data$row), ]
saveRDS(train_data, file = "./Chunk8-Data Analysis/scalability/scam_train_data.rds")
saveRDS(test_data, file = "./Chunk8-Data Analysis/scalability/scam_test_data.rds")

########################## scam
#### time of estimation
estimation_time_model <- scam_function(step = "estimation", feature = "time")
est_time_table <- bind_prediction_result(estimation_time_model)
est_time_correlation <- correlation_plot(est_time_table)
est_time_correlation$correlation_plot
#### memory of estimation
estimation_memory_model <- scam_function(step = "estimation", feature = "memory")
est_memory_table <- bind_prediction_result(estimation_memory_model)
est_memory_correlation <- correlation_plot(est_memory_table)
est_memory_correlation$correlation_plot
#### time of simulation
simulation_time_model <- scam_function(step = "simulation", feature = "time")
sim_time_table <- bind_prediction_result(simulation_time_model)
sim_time_correlation <- correlation_plot(sim_time_table)
sim_time_correlation$correlation_plot
#### memory of simulation
simulation_memory_model <- scam_function(step = "simulation", feature = "memory")
sim_memory_table <- bind_prediction_result(simulation_memory_model)
sim_memory_correlation <- correlation_plot(sim_memory_table)
sim_memory_correlation$correlation_plot
scam_model <- list("estimation_time_model" = list("estimation_time_model" = estimation_time_model,
                                                  "est_time_table" = est_time_table,
                                                  "est_time_correlation" = est_time_correlation),
                   "estimation_memory_model" = list("estimation_memory_model" = estimation_memory_model,
                                                  "est_memory_table" = est_memory_table,
                                                  "est_memory_correlation" = est_memory_correlation),
                   "simulation_time_model" = list("simulation_time_model" = simulation_time_model,
                                                  "sim_time_table" = sim_time_table,
                                                  "sim_time_correlation" = sim_time_correlation),
                   "simulation_memory_model" = list("simulation_memory_model" = simulation_memory_model,
                                                  "sim_memory_table" = sim_memory_table,
                                                  "sim_memory_correlation" = sim_memory_correlation))
saveRDS(scam_model, file = "./Chunk8-Data Analysis/scalability/scam_model.rds")

ggsave(plot = wrap_plots(map(scam_model, function(x){x[[3]][[2]]}), ncol = 2),
       filename = "./Chunk8-Data Analysis/scalability/scam_model.pdf",
       height = 8,
       width = 8,
       units = "in")

########################## RF
#### time of estimation
estimation_time_model <- scam_function(step = "estimation", feature = "time", model_selected = "RF")
est_time_table <- bind_prediction_result(estimation_time_model)
est_time_correlation <- correlation_plot(est_time_table)
est_time_correlation$correlation_plot
#### memory of estimation
estimation_memory_model <- scam_function(step = "estimation", feature = "memory", model_selected = "RF")
est_memory_table <- bind_prediction_result(estimation_memory_model)
est_memory_correlation <- correlation_plot(est_memory_table)
est_memory_correlation$correlation_plot
#### time of simulation
simulation_time_model <- scam_function(step = "simulation", feature = "time", model_selected = "RF")
sim_time_table <- bind_prediction_result(simulation_time_model)
sim_time_correlation <- correlation_plot(sim_time_table)
sim_time_correlation$correlation_plot
#### memory of simulation
simulation_memory_model <- scam_function(step = "simulation", feature = "memory", model = "RF")
sim_memory_table <- bind_prediction_result(simulation_memory_model)
sim_memory_correlation <- correlation_plot(sim_memory_table)
sim_memory_correlation$correlation_plot
RF_model <- list("estimation_time_model" = list("estimation_time_model" = estimation_time_model,
                                                "est_time_table" = est_time_table,
                                                "est_time_correlation" = est_time_correlation),
                 "estimation_memory_model" = list("estimation_memory_model" = estimation_memory_model,
                                                  "est_memory_table" = est_memory_table,
                                                  "est_memory_correlation" = est_memory_correlation),
                 "simulation_time_model" = list("simulation_time_model" = simulation_time_model,
                                                "sim_time_table" = sim_time_table,
                                                "sim_time_correlation" = sim_time_correlation),
                 "simulation_memory_model" = list("simulation_memory_model" = simulation_memory_model,
                                                  "sim_memory_table" = sim_memory_table,
                                                  "sim_memory_correlation" = sim_memory_correlation))
saveRDS(RF_model, file = "./Chunk8-Data Analysis/scalability/RF_model.rds")
ggsave(plot = wrap_plots(map(RF_model, function(x){x[[3]][[2]]}), ncol = 2),
       filename = "./Chunk8-Data Analysis/scalability/RF_model.pdf",
       height = 8,
       width = 8,
       units = "in")
################################################################################
#######################            Plot Rect Data      #########################
################################################################################
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

### scalability score

#### aggregate by repeat time
score_data <- scalability_data %>% 
  group_by(method, cell_num, gene_num) %>% 
  summarise(
    estimation_time = mean(estimation_time, na.rm = TRUE),
    estimation_memory = mean(estimation_memory, na.rm = TRUE),
    simulation_time = mean(simulation_time, na.rm = TRUE),
    simulation_memory = mean(simulation_memory, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  filter(cell_num == 100 & gene_num == 1000 |
         cell_num == 1000 & gene_num == 1000 |
         cell_num == 1000 & gene_num == 10000 |
         cell_num == 10000 & gene_num == 1000) %>% 
  mutate(
    estimation_time_text = case_when(
      estimation_time < 1 ~ "<1s",
      estimation_time >= 1 & estimation_time < 60 ~ paste0(as.character(round(estimation_time)), "s"),
      estimation_time >= 60 & estimation_time < 3600 ~ paste0(estimation_time %/% 60, "m"),
      estimation_time > 3600 ~ paste0(">", estimation_time %/% 3600, 'h')
    ),
    simulation_time_text = case_when(
      simulation_time < 1 ~ "<1s",
      simulation_time >= 1 & simulation_time < 60 ~ paste0(as.character(round(simulation_time)), "s"),
      simulation_time >= 60 & simulation_time < 3600 ~ paste0(simulation_time %/% 60, "m"),
      simulation_time > 3600 ~ paste0(">", simulation_time %/% 3600, 'h')
    ),
    estimation_memory_text = case_when(
      estimation_memory < 10 ~ paste0(round(estimation_memory, digits = 1), 'M'),
      estimation_memory >= 10 & estimation_memory < 1024 ~ paste0(round(estimation_memory), 'M'),
      estimation_memory >= 1024 ~ paste0(">", estimation_memory %/%1024, 'G')
    ),
    simulation_memory_text = case_when(
      simulation_memory < 10 ~ paste0(round(simulation_memory, digits = 1), 'M'),
      simulation_memory < 1024 ~ paste0(round(simulation_memory), 'M'),
      simulation_memory >= 1024 ~ paste0(">", simulation_memory %/%1024, 'G')
    )
  ) %>% 
  group_by(cell_num, gene_num) %>% 
  mutate(
    across(all_of(c("estimation_time",
                    "estimation_memory",
                    "simulation_time",
                    "simulation_memory")), ~ 1 - (pnorm((.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE))))
  ) %>% 
  ungroup() %>% 
  mutate(
    estimation_score = (estimation_time + estimation_memory) / 2,
    simulation_score = (simulation_time + simulation_memory) / 2,
    estimation_text = paste0(estimation_time_text, "\n", estimation_memory_text),
    estimation_text = case_when(
      estimation_text == "NA\nNA" ~ "",
      TRUE ~ estimation_text
    ),
    simulation_text = paste0(simulation_time_text, "\n", simulation_memory_text),
    simulation_text = case_when(
      simulation_text == "NA\nNA" ~ "",
      TRUE ~ simulation_text
    ),
    data = paste0(as.character(cell_num), "x", as.character(gene_num))
  )

data_100_1000_est <- score_data %>% 
  filter(cell_num == 100, gene_num == 1000) %>%
  select(method, estimation_score, estimation_text)
colnames(data_100_1000_est) <- c("method", paste0("scalability_100_1000_", colnames(data_100_1000_est)[2:3]))

data_100_1000_sim <- score_data %>% 
  filter(cell_num == 100, gene_num == 1000) %>%
  select(method, simulation_score, simulation_text)
colnames(data_100_1000_sim) <- c("method", paste0("scalability_100_1000_", colnames(data_100_1000_sim)[2:3]))

data_1000_1000_est <- score_data %>% 
  filter(cell_num == 1000, gene_num == 1000) %>%
  select(method, estimation_score, estimation_text)
colnames(data_1000_1000_est) <- c("method", paste0("scalability_1000_1000_", colnames(data_1000_1000_est)[2:3]))

data_1000_1000_sim <- score_data %>% 
  filter(cell_num == 1000, gene_num == 1000) %>%
  select(method, simulation_score, simulation_text)
colnames(data_1000_1000_sim) <- c("method", paste0("scalability_1000_1000_", colnames(data_1000_1000_sim)[2:3]))

data_1000_10000_est <- score_data %>% 
  filter(cell_num == 1000, gene_num == 10000) %>%
  select(method, estimation_score, estimation_text)
colnames(data_1000_10000_est) <- c("method", paste0("scalability_1000_10000_", colnames(data_1000_10000_est)[2:3]))

data_1000_10000_sim <- score_data %>% 
  filter(cell_num == 1000, gene_num == 10000) %>%
  select(method, simulation_score, simulation_text)
colnames(data_1000_10000_sim) <- c("method", paste0("scalability_1000_10000_", colnames(data_1000_10000_sim)[2:3]))

data_10000_1000_est <- score_data %>% 
  filter(cell_num == 10000, gene_num == 1000) %>%
  select(method, estimation_score, estimation_text)
colnames(data_10000_1000_est) <- c("method", paste0("scalability_10000_1000_", colnames(data_10000_1000_est)[2:3]))

data_10000_1000_sim <- score_data %>% 
  filter(cell_num == 10000, gene_num == 1000) %>%
  select(method, simulation_score, simulation_text)
colnames(data_10000_1000_sim) <- c("method", paste0("scalability_10000_1000_", colnames(data_10000_1000_sim)[2:3]))

score_data <- data_100_1000_est %>% 
  full_join(data_100_1000_sim, by = "method")%>% 
  full_join(data_1000_1000_est, by = "method")%>% 
  full_join(data_1000_1000_sim, by = "method")%>% 
  full_join(data_1000_10000_est, by = "method")%>% 
  full_join(data_1000_10000_sim, by = "method")%>% 
  full_join(data_10000_1000_est, by = "method")%>% 
  full_join(data_10000_1000_sim, by = "method")

saveRDS(score_data, file = "./Chunk8-Data Analysis/scalability/scalability_plot_data.rds")
