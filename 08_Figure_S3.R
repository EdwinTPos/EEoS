# =====================================================================
# Code for
# A coherent signature of disturbance in successional forests
# Code written and maintained by Edwin Pos and Erica Newman
# E.T.Pos@uu.nl, ericanewman.sci@gmail.com
# Updated July 13, 2026
# =====================================================================

# =========================================
# FIGURE S3 state variables S, N, E by plot
# =========================================

#Clear environment
#rm(list=ls())
#rm(list = ls(all.names = TRUE))

#Load libraries and dependencies
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

# Load the data
trees <- read.csv("trees.csv", header=T)  

#Global Setup
desired_order <- c("SV", "LEPviejo", "CR", "FEB", "LEP", "TIR", "LSUR", "JE")

#Prepare Data
trees$Plot <- factor(trees$Plot, levels = desired_order)
trees$S <- as.numeric(trees$S)
trees$N <- as.numeric(trees$N)
trees$E <- as.numeric(trees$E)
trees$Evenness <- as.numeric(trees$Evenness)

#Custom color palette

my_colors <-     c("CR" =       "#51310B", 
                   "FEB" =      "#8E4F00", 
                   "JE" =       "#DEC27C",
                   "LEP" =      "#999933",
                   "LEPviejo" = "#117733", 
                   "LSUR" =     "#7CCCBE",  
                   "SV" =       "#33BBEE", 
                   "TIR" =      "#332288")
custom_colors <- my_colors

my_shapes <- c("CR" =       21, 
               "FEB" =      8, 
               "JE" =        24, 
               "LEP" =      25, 
               "LEPviejo" = 12, 
               "LSUR" =     23, 
               "SV" =        22, 
               "TIR" =       3)

plot_list <- trees %>%
  pivot_longer(cols = c(S, N, E, Evenness), names_to = "Variable", values_to = "Value") %>%
  mutate(
    Value = case_when(
      Variable == "E" ~ Value,
      Variable == "S" ~ Value,
      TRUE            ~ Value
    ),
    Plot = factor(Plot, levels = desired_order)
  ) %>%
  group_split(Variable) %>%
  lapply(function(df_var) {
    var_name <- unique(df_var$Variable)
    ggplot(df_var, aes(x = Year, y = Value, color = Plot, fill = Plot, shape = Plot)) +
      geom_line(aes(color = Plot), linewidth = 1, alpha = 0.5) +
      geom_point(size = 2, stroke = 0.5, alpha = .8) +
      scale_color_manual(values = my_colors) +
      scale_fill_manual(values = my_colors) +
      scale_shape_manual(values = my_shapes) +
      scale_x_continuous(labels = function(x) round(x, 0)) +
      labs(
        title = paste("State Variable:", var_name),
        x = "Year",
        y = var_name
      ) +
      theme_bw(base_size = 11) +
      theme(
        aspect.ratio = 1,
        legend.position = "right",
        legend.title = element_text(size = 12, face = "plain") 
      )
  })

# 2x2 layout: plots 1 & 2 on top, plot 3 + collected legend on bottom
fig_s3 <- (plot_list[[1]] + plot_list[[2]] +
                        plot_list[[3]] + plot_list[[4]]) + 
                       plot_layout(ncol = 2, guides = "collect")
# View and Save
#fig_s3
ggsave("figure_S3.png", plot = fig_s3, width = 7.5, height = 6.5, dpi = 300)
