library(ggplot2)
library(tidyverse)

data_info <- openxlsx::read.xlsx("./Chunk1-Data preparation/evaluation_datasets.xlsx", rowNames = TRUE, sep.names = " ") %>% 
  mutate(
    Data = str_split(`Dataset ID`, "_", simplify = TRUE)[, 1],
    Element = `Cell/Spot Number` * `Gene Number`
  ) %>% 
  select(Data, Platform, Element)
size_quantile <- quantile(data_info$Element) %>% unname()

accuracy_long_data <- readRDS("./Chunk8-Data Analysis/accuracy/accuracy_long_data.rds") %>% 
  left_join(data_info, by = "Data") %>% 
  mutate(
    Size = case_when(
      Element >= size_quantile[1] & Element < size_quantile[2] ~ "small",
      Element >= size_quantile[2] & Element < size_quantile[3] ~ "mid-to-small",
      Element >= size_quantile[3] & Element < size_quantile[4] ~ "mid-to-large",
      Element >= size_quantile[4] ~ "large"
    )
  )

data_size_score <- accuracy_long_data %>% 
  group_by(Method, metric, Data, Size) %>% 
  summarise(
    value = mean(value, na.rm = TRUE),
  ) %>% 
  ungroup() %>% 
  group_by(Method, Data, Size) %>% 
  summarise(
    value = mean(value, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  group_by(Method, Size) %>% 
  summarise(
    value = mean(value, na.rm = TRUE)
  ) %>% 
  ungroup()


data_size_score$Size <- factor(data_size_score$Size, levels = c("small", "mid-to-small", "mid-to-large", "large"))
### barplot
data_size_score_plot <- ggplot(data_size_score) +
  geom_col(aes(x = Size, y = value, fill = Size), width = 0.7) +
  facet_wrap(Method ~ .,ncol = 5) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.ticks = element_line(linewidth = unit(0.15, "cm")),
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.position = "top"
  ) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(4, "Set2")) +
  scale_y_continuous(expand = c(0, 0), labels = seq(0, 0.8, 0.2), limits = c(0, 0.85)) +
  guides(fill = guide_legend(title = "Data Size")) +
  ylab("Aggregated scores")
ggsave(plot = data_size_score_plot,
       filename = "../sim-article/figures/Supp_Fig.pdf",
       width = 9,
       height = 12)







