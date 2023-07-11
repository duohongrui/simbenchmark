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