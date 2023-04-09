library(stringr)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)

### usability data
usability_data <- openxlsx::read.xlsx("Chunk7-Usability/method_scoring.xlsx", rowNames = TRUE)
colnames(usability_data) <- gsub(x = colnames(usability_data), pattern = "[.]", replacement = " ")
cite_score <- 1 - seq(0, 1, 1/42)
cite_score <- cite_score[order(usability_data$Cite, decreasing = TRUE)]
usability_data$Cite <- cite_score

### weight data
weight_data <- openxlsx::read.xlsx("Chunk7-Usability/Usability-scoringsheet.xlsx")

### summarize usability score
usability_score <- as.matrix(usability_data) %*% matrix(weight_data$Weight, ncol = 1)/100

usability <- tibble(method = rownames(usability_data),
                    usability_score = usability_score[,1])
saveRDS(usability, file = "Chunk8-Data Analysis/usability/usability.rds")

##########################################################
###################   SUMMARY FIGURE   ###################
##########################################################

########################   detailed plot   ##################################
usability_long_data <- usability_data %>% 
  tibble::rownames_to_column(var = "method") %>% 
  pivot_longer(2:ncol(.), names_to = "content", values_to = "score") %>% 
  mutate(
    category = rep(c(rep("Avalability", 8),
                     rep("Code", 8),
                     rep("Evaluation", 4),
                     rep("Maintenance", 5),
                     rep("Documentation", 5),
                     rep("Paper", 4)), 42)
  )

rect_width <- 0.5
big_space <- 1
tiny_space <- 0.2

usability_long_data <- usability_long_data %>% 
  mutate(
    x = rep(cumsum(rep(rect_width, 42)), each = 34),
    xmin = x - rect_width/2,
    xmax = x + rect_width/2,
    spacing = c(0, diff(as.numeric(factor(category)))) != 0,
    width = case_when(
      spacing ~ big_space,
      TRUE ~ tiny_space
    ),
    ymax = rep(c(1:34) * rect_width + cumsum(width[1:34]), 42),
    ymin = ymax - 0.5,
    y = ymin + 0.25,
  )

### text positions

text_data <- tibble(
  text = c(unique(usability_long_data$method),
           unique(usability_long_data$content),
           unique(usability_long_data$category)),
  x = c(unique(usability_long_data$x),
        rep(-0.1, 34),
        rep(-6, 6)),
  y = c(rep(-0.5, 42),
        unique(usability_long_data$y),
        c(2.55, 9.0, 14, 18.3, 22.55, 26.3)),
  angle = c(rep(45, 42),
            rep(0, 34),
            rep(90, 6)),
  # size = c(rep(5, 42), rep(10, 6)),
  hjust = c(rep("right", 42),
            rep("right", 34),
            rep("middle", 6))
)

### method segment
segment_data <- tibble('x' = c(unique(usability_long_data$x), rep(-5, 6)),
                       'xend' = c(unique(usability_long_data$x), rep(-5, 6)),
                       'y' = c(rep(-0.01, 42), unique(usability_long_data$ymax)[c(8, 16, 20, 25, 30, 34)]),
                       'yend' = c(rep(-0.40, 42), unique(usability_long_data$ymin)[c(1, 9, 17, 21, 26, 31)]))


g <- ggplot()+
  geom_rect(aes(xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = ymax,
                fill = score,),
            usability_long_data,
            size = 0.5,
            color = "white")+
  geom_text(aes(x = x,
                y = y,
                label = text,
                hjust = hjust,
                angle = angle),
            text_data)+
  geom_segment(aes(x = x,
                   xend = xend,
                   y = y,
                   yend = yend),
               size = 0.5,
               color = "black",
               linetype = "solid",
               segment_data)+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "bottom")+
  scale_fill_gradient(low = "#4B1248", high = "#F0C27B")+
  expand_limits(x = c(min(usability_long_data$xmin) - 0.2,
                      max(usability_long_data$xmax) + 0.2),
                y = c(min(usability_long_data$ymin) - 2,
                      max(usability_long_data$ymax) + 1))


########################   usability barplot   ##################################

usability_score_data <- data.frame("method" = rownames(usability_score),
                              "score" = usability_score[, 1])

barplot <- ggplot(usability_score_data, aes(x = reorder(method, -score), y = score, fill = score))+
  geom_bar(stat = "identity", width = 0.8)+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, 0.25), limits = c(0, 1))+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank()
  )+
  ylab(label = "Usability Score")+
  scale_fill_gradient(low = "#4B1248", high = "#F0C27B")


########################   plot patch   ##################################

layout <- "
#AAAAAAAAAAA#
#BBBBBBBBBBBB
"

usability_plot <- wrap_plots(barplot, g, 
                             ncol = 1,
                             heights = c(1,3),
                             design = layout)+
  plot_annotation(tag_levels = 'A')
usability_plot
ggsave(usability_plot, filename = "Chunk8-Data Analysis/usability/usability.pdf", width = 15, height = 15, units = "in")


