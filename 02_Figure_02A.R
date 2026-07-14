# =====================================================================
# Code for
# A coherent signature of disturbance in successional forests
# Code written and maintained by Edwin Pos and Erica Newman
# E.T.Pos@uu.nl, ericanewman.sci@gmail.com
# Updated July 13, 2026
# =====================================================================

# =================================================
# FIGURE 2A (OVERVIEW BIOMASS PREDICTION ALL PLOTS)
# =================================================

#Clear environment, set working directory here and load libraries
#rm(list=ls())
#rm(list = ls(all.names = TRUE))

#install.packages("ggraph")

# Load libraries and dependencies
library(ggplot2)
library(ggraph)
library(dplyr)
library(stringr)

# Load the data and inspect
trees <- read.csv("trees.csv", header=T)

# cleaned file should come back with these numbers
dim(trees) # 136 12

unique_plots <- unique(trees$Plot)
unique_plots # should not have BEJ and should instead have FEB, LSUR_new should not be included either

# Custom color palette

my_colors <-     c("CR" =       "#51310B", 
                   "FEB" =      "#8E4F00", 
                   "JE" =       "#DEC27C",
                   "LEP" =      "#999933",
                   "LEPviejo" = "#117733", 
                   "LSUR" =     "#7CCCBE",  
                   "SV" =       "#33BBEE", 
                   "TIR" =      "#332288")


my_shapes <- c("CR" =       21, 
               "FEB" =      8, 
               "JE" =        24, 
               "LEP" =      25, 
               "LEPviejo" = 12, 
               "LSUR" =     23, 
               "SV" =        22, 
               "TIR" =       3)


trees$Plot <- factor(trees$Plot)

# Reorder the 'Plot' factor levels based on the 'B'  
# (using -value to get descending order) highest maximum to lowest maximum

trees$Plot <- reorder(trees$Plot, -trees$B, FUN = max, na.rm = TRUE)

figure_02A <- ggplot(
  trees,
  aes(log(pBnum_exact), log(B), color = Plot, shape = Plot)
) +
  geom_point(size = 1, stroke = 1, alpha = 0.7) +
  
  # To merge into ONE legend, 'name' must be exactly the same in both scales
  scale_color_manual(name = "Plot", values = my_colors) +
  scale_shape_manual(name = "Plot", values = my_shapes) +
  
  coord_cartesian(xlim = c(10.3, 12), ylim = c(10.3, 12)) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray20", alpha = 0.5) +
  
  theme_bw(base_size = 11) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "plain") 
  ) +
  labs(x = "Ln(predicted biomass)", y = "Ln(observed biomass)")

#figure_02A

ggsave("figure_02A.png", 
       plot = figure_02A, 
       width = 5,    # Width of the plot
       height = 3.5,   # Height of the plot
       units = "in", # Units: "in", "cm", "mm", or "px"
       dpi = 300)    # Resolution
