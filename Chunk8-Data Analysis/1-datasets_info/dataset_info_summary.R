################################################################################
##########################   Supplementary Figure 1  ###########################
################################################################################
data <- openxlsx::read.xlsx("/Users/duohongrui/Desktop/simbenchmark/Chunk1-Data preparation/evaluation_datasets.xlsx", sep.names = " ")

library(dplyr)
library(stringr)
library(dplyr)
library(ggplot2)
library(patchwork)

colors <- c("#ec671e", "#2e8d3a", "#1d5e9b", "#774d9c", "#29a8b9", "#9cb5db", "#f4a965",
            "#ee7e80", "#b38880", "#f5c290", "#838384", "#cdce79", "#be1a22", "#f4c33f",
            "#37b380", "#4c98c8", "#36546d", "#b98519", "#dc4530", "#B3DE69", "#5d5098",
            "#edec73", "#f18a1c", "#0f86a9")

### single-cell platforms
singlecell_data <- data %>% 
  slice(1:101) %>% 
  mutate(
    Platform = case_when(
      Platform == "Smart-seq2\r\n10X Genomics" ~ "Mix sources1",
      Platform == "CEL-seq\r\nCEL-seq2" ~ "Mix sources2",
      TRUE ~ Platform
    )
  )

platforms <- table(singlecell_data$Platform) %>% as.data.frame()
platforms <- platforms %>% 
  mutate(proportion = Freq/sum(Freq))

p1 <- ggplot(platforms, mapping = aes(x = "", y = Freq, fill = Var1))+
  geom_bar(stat = "identity", position = "stack", width = 1, color = "white")+
  coord_polar(theta = "y")+
  scale_fill_manual(values = colors)+
  labs(x = "", y = "", title = "")+
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none",
        legend.title = element_blank())+
  geom_text(aes(x = sum(Freq)/65, y = Freq/2 + c(0, cumsum(Freq)[-length(Freq)]),
                label = paste0(Var1, "\n",
                               "(n=", Freq, ", ", sprintf("%0.2f", proportion*100), "%)")))

### spatial platforms
spatial_data <- data %>% 
  slice(102:152)

platforms <- table(spatial_data$Platform) %>% as.data.frame()
platforms <- platforms %>% 
  mutate(proportion = Freq/sum(Freq))

p2 <- ggplot(platforms, mapping = aes(x = "", y = Freq, fill = Var1))+
  geom_bar(stat = "identity", position = "stack", width = 1, color = "white")+
  coord_polar(theta = "y")+
  scale_fill_manual(values = colors[12:23])+
  labs(x = "", y = "", title = "")+
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none",
        legend.title = element_blank())+
  geom_text(aes(x = sum(Freq)/30, y = Freq/2 + c(0, cumsum(Freq)[-length(Freq)]),
                label = paste0(Var1, "\n",
                               "(n=", Freq, ", ", sprintf("%0.2f", proportion*100), "%)")))

### quantification of counts for single-cell technique
quanti <- table(data$`Quantification Strategy`) %>% as.data.frame()
p3 <- ggplot(quanti, mapping = aes(x = Var1, y = Freq, fill = Var1))+
  geom_bar(stat = "identity", color = "black", alpha = 0.9, width = 0.6)+
  scale_fill_manual(values = colors[5:10])+
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  theme_classic()+
  geom_text(aes(label = paste0("n=", Freq, "\n(", sprintf("%0.2f", (Freq/sum(Freq))*100), "%)")), nudge_y = 5)+
  theme(axis.text = element_text(color = "black"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        legend.title = element_blank()) +
  xlab("Quantification strategy of \ncounts for scRNA-seq data") +
  ylab("Number of methods") +
  ylim(0, 62)


### cell or gene number
p4 <- ggplot(data, aes(x = `Cell/Spot Number`, y = `Gene Number`)) +
  geom_point(aes(color = Platform)) +
  scale_color_manual(values = colors) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text = element_text(color = "black"),
        panel.grid = element_blank()) +
  guides(color = guide_legend(ncol = 6))



### method category
method <- openxlsx::read.xlsx("Chunk1-Data preparation/methods.xlsx", sheet = 1, colNames = TRUE, sep.names = " ") %>% 
  select(-2) %>% 
  slice(1:49)
##### 1. model
model <- table(method$`Model Category`) %>% as.data.frame()
p5 <- ggplot(model, mapping = aes(x = "", y = Freq, fill = Var1))+
  geom_bar(stat = "identity", position = "stack", width = 1, color = "white")+
  coord_polar(theta = "y")+
  scale_fill_manual(values = colors)+
  labs(x = "", y = "", title = "")+
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none",
        legend.title = element_blank())+
  geom_text(aes(x = sum(Freq)/30, y = Freq/2 + c(0, cumsum(Freq)[-length(Freq)]),
                label = paste0(Var1, "\n",
                               "(n=", Freq, ", ", sprintf("%0.2f", (Freq/sum(Freq))*100), "%)")))


##### 2. functionality
category <- data.frame("y" = c(36, 22, 19, 15),
                       "label" = c("Cell group", "DEGs", "Cell batch", "Cellular differentiation \ntrajectory")) %>% 
  mutate(prop = y/49)

p6 <- ggplot(category, mapping = aes(x = label, y = y, fill = label))+
  geom_bar(stat = "identity", color = "black", alpha = 0.9, width = 0.6)+
  scale_fill_manual(values = colors[5:10])+
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  theme_classic()+
  geom_text(aes(label = paste0("n=", y, "\n(", sprintf("%0.2f", prop*100), "%)")), nudge_y = 5)+
  theme(axis.text = element_text(color = "black"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        legend.title = element_blank()) +
  xlab("Method Functionality") +
  ylab("Number of methods") +
  ylim(0, 42)



ggsave(
  filename = "../sim-article/figures/Supp_Fig_1.pdf",
  plot = (p1+p2+p3)/(p5+p6)+plot_layout(heights = c(1,1,1)),
  height = 15,
  width = 10,
  units = "in"
)

ggsave(
  filename = "../sim-article/figures/Supp_Fig_2.pdf",
  plot = p4,
  height = 10,
  width = 9,
  units = "in"
)
