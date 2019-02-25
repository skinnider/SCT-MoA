library(ggplot2)
library(ggridges)
library(ggsci)
library(viridis)
library(cowplot)
library(colorspace)
library(Polychrome)
library(RColorBrewer)
library(ggalt)
library(scales)
library(drlib)
library(scico)
library(patchwork)

# Define theme
theme_sc = theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        panel.grid = element_blank(),
        panel.border = element_blank(), 
        axis.line.y = element_line(colour = "grey50"),
        axis.line.x = element_line(colour = "grey50"), 
        axis.ticks = element_line(colour="grey50"),
        legend.position = "top",
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(size = 10, hjust = 0.5))

clean_theme = theme_bw() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        panel.grid = element_blank(),
        panel.border = element_blank(), 
        strip.text = element_text(size = 8),
        strip.background = element_rect(color = "grey50", size = 0),
        axis.line.y = element_line(colour = "grey50"),
        axis.line.x = element_line(colour = "grey50"), 
        axis.ticks = element_line(colour = "grey50"),
        legend.position = "top",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(size = 12, hjust = 0.5))

darken = function(color, factor=1.4){
  col = col2rgb(color)
  col = col/factor
  col = rgb(t(col), maxColorValue=255)
  col
}

lighten = function(color, factor=1.4){
  col = col2rgb(color)
  col = col*factor
  col = rgb(t(col), maxColorValue=255)
  col
}
