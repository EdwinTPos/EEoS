# =====================================================================
# Code for
# A coherent signature of disturbance in successional forests
# Code written and maintained by Edwin Pos and Erica Newman
# E.T.Pos@uu.nl, ericanewman.sci@gmail.com
# Updated October 6, 2026
# =====================================================================

#Clear environment
#rm(list=ls())
#rm(list = ls(all.names = TRUE))

# Load libraries and dependencies
library(ggplot2)
library(ggraph)
library(dplyr)
library(stringr)
library(ggrepel)
library(patchwork)

########################################################################
### FIGURE 2B: Distance between prediction and observation over time ###
########################################################################

#Read cleaned data from 01_EEoS_Analysis2026 file
trees <- read.csv("trees.csv", header=T)

# Biomasss over time, labeled by foresttype (with TIR standing out)
trees$group_col <- ifelse(
  trees$Plot == "TIR", "TIR",
  ifelse(trees$Foresttype == "primary", "primary", "secondary"))

trees$Category = paste(trees$Foresttype,","," ",trees$group_col, sep="")
trees$Category[trees$Category == "primary, primary"] <- "primary forest"
trees$Category[trees$Category == "secondary, secondary"] <- "secondary forest"

# Plot change in state variables over time, first set variable (B,E,S or N) and then run plot code

# Custom color palette

my_colors <-     c("CR" =       "#51310B", 
                   "FEB" =      "#8E4F00", 
                   "JE" =       "#DEC27C",
                   "LEP" =      "#999933",
                   "LEPviejo" = "#117733", 
                   "LSUR" =     "#7CCCBE",  
                   "SV" =       "#33BBEE", 
                   "TIR" =      "#332288")
custom_colors <- my_colors

my_shapes <- c("CR" =       16,   # filled circle, secondary
               "FEB" =      16,   # filled circle, secondary
               "JE" =       16,   # filled circle, secondary
               "LEP" =      16,   # filled circle, secondary
               "LEPviejo" = 15,   # filled square -> primary
               "LSUR" =     16,   # filled circle, secondary
               "SV" =       15,   # filled square -> primary
               "TIR" =      16)   # filled circle, secondary

# Reorder the 'Plot' factor levels based on the (StatVar) column
# (using -value to get descending order) highest to lowest

trees$Plot <- reorder(trees$Plot, -trees[["dist"]], FUN = max, na.rm = TRUE)
ages <- data.frame(Plot=c("JE","FEB","LSUR","TIR","LEP","CR","LEPviejo","SV"),Age=c(10,10,12,15,20,20,200,200))

# Add ages for each forest plot location 
trees <- trees %>%
  left_join(ages, by = "Plot") %>%
  group_by(Plot) %>%
  mutate(
    first_year = min(Year),
    Age_new = Age + (Year - first_year)
  ) %>%
  ungroup() %>%
  select(-first_year)

# Tag first and last year within each plot
age_labels <- trees %>%
  group_by(Plot) %>%
  filter(Year == min(Year) | Year == max(Year)) %>%
  ungroup()

age_labels[age_labels$Foresttype == "Primary","Age_new"] <- NA

figure2B <- ggplot(trees, aes(Year, abs(.data[["dist"]]),
                               color = Plot, shape = Plot,
                               fill = Foresttype)) +
  geom_line(aes(color = Plot), linewidth = 2, alpha = 0.5) +
  geom_point(size = 2, stroke = 1, alpha = 0.7) +
  geom_text_repel(
  data = age_labels,
  aes(label = Age_new),
  size = 3, fontface = "bold",
  direction = "y",
  nudge_x = -0.5,              # shift labels left; adjust to taste
  hjust = -1,                   # right-align text so it ends near the nudged spot
  min.segment.length = Inf,    # never draw connector lines
  show.legend = FALSE
)+
   scale_color_manual(name = "Plot", values = my_colors,
                     guide = guide_legend(
                       order = 1,
                       override.aes = list(shape = my_shapes, size = 3)
                     )) +
  scale_shape_manual(name = "Plot", values = my_shapes, guide = "none") +
  scale_fill_manual(
    name = "Forest type",
    values = c(Primary = "grey20", Secondary = "grey20"),
    guide = guide_legend(
      order = 2,
      override.aes = list(shape = c(15,16), color = "black",
                          fill = c("grey20", "grey20")))
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "plain") 
  ) +
  labs(x = "Year", y = "[predicted-observed]ln(B)")

#Save figure
ggsave("Figure_02B.png", 
       plot = figure2B, 
       width = 5,    # Width of the plot
       height = 4,   # Height of the plot
       units = "in", # Units: "in", "cm", "mm", or "px"
       dpi = 300)    # Resolution

# Combine figures 2A and 2B if needed (note that this only works if figure_02A is loaded in workspace)
# Remove the legend from figure 2b for easy plotting.
figure2B <- figure2B + guides(colour = "none", shape = "none", linetype = "none")

combined_figure2 <- (
  figure_02A / figure2B
) +
  plot_layout(guides = "collect") +                             
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 12, face = "bold"),
    legend.position = "right"                                   
  )

ggsave("Figure_02_combined.png", 
       plot = combined_figure2, 
       width = 5,    # Width of the plot
       height = 5,   # Height of the plot
       units = "in", # Units: "in", "cm", "mm", or "px"
       dpi = 300)    # Resolution
