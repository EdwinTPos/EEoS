# =====================================================================
# Code for
# A coherent signature of disturbance in successional forests
# Code written and maintained by Edwin Pos and Erica Newman
# E.T.Pos@uu.nl, ericanewman.sci@gmail.com
# Updated July 13, 2026
# =====================================================================

# ====================================
# FIGURE 3: Biomass increase over time
# ====================================

#Clear environment, set working directory here and load libraries
#rm(list=ls())
#rm(list = ls(all.names = TRUE))

# Load libraries and dependencies
library(ggplot2)
library(dplyr)
library(ggrepel)

#Read cleaned data from 01_EEoS_Analysis2026 file
trees <- read.csv("trees.csv", header=T)

# Biomasss over time, labeled by foresttype (with TIR standing out)
trees$group_col <- ifelse(
  trees$Plot == "TIR", "TIR",
  ifelse(trees$Foresttype == "Primary", "Primary", "Secondary"))

trees$Category = paste(trees$Foresttype,","," ",trees$group_col, sep="")
trees$Category[trees$Category == "Primary, Primary"] <- "primary forest"
trees$Category[trees$Category == "Secondary, Secondary"] <- "secondary forest"

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

#Choose which state variable to plot; manuscript shows "B" for Figure 03
StatVar = "B"
#StatVar = "S"
#StatVar = "N"
#StatVar = "E"

# Reorder the 'Plot' factor levels based on the 'B' (StatVar) column
# (using -value to get descending order) highest to lowest

trees$Plot <- reorder(trees$Plot, -trees[[StatVar]], FUN = max, na.rm = TRUE)
ages <- data.frame(Plot=c("JE","FEB","LSUR","TIR","LEP","CR","LEPviejo","SV"),Age=c(10,10,12,15,20,20,200,200))

trees2 <- trees %>%
  left_join(ages, by = "Plot") %>%
  group_by(Plot) %>%
  mutate(
    first_year = min(Year),
    Age_new = Age + (Year - first_year)
  ) %>%
  ungroup() %>%
  select(-first_year)

# Tag first and last year within each plot
age_labels <- trees2 %>%
  group_by(Plot) %>%
  filter(Year == min(Year) | Year == max(Year)) %>%
  ungroup()

age_labels[age_labels$Foresttype == "Primary","Age_new"] <- NA

figure3 <- ggplot(trees2, aes(Year, log(.data[[StatVar]]), color = Plot, shape = Plot,
                             fill = Foresttype)) +
  geom_line(aes(color = Plot), linewidth = 2, alpha = 0.5) +
  geom_point(size = 2, stroke = 1, alpha = 0.7) +
  geom_text_repel(
    data = age_labels,
    aes(label = Age_new),
    size = 3, fontface = "bold",
    direction = "y",
    min.segment.length = 0,
    segment.alpha = 0.4, segment.size = 0.3,
    show.legend = FALSE
  ) +
  scale_color_manual(name = "Plot", values = my_colors,
                     guide = guide_legend(
                       order = 1,
                       override.aes = list(shape = my_shapes, size = 2)
                     )) +
  scale_shape_manual(name = "Plot", values = my_shapes, guide = "none") +
  scale_fill_manual(
    name = "Forest type",
    values = c(Primary = "grey20", Secondary = "grey20"),
    guide = guide_legend(
      order = 2,
      keywidth = unit(0.9, "line"), keyheight = unit(0.9, "line"),
      override.aes = list(shape = c(15,16), color = "black",
                          fill = c("grey20", "grey20")))
  ) +
  xlab("Year") +
  ylab("Ln(observed biomass)") +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = .5, vjust = 0.5, size = 12),
    axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.key.size = unit(0.8, "line"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "plain"),
    legend.spacing.y = unit(0.5, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.position = "right"
  )
#figure3

#Save figure
ggsave("Figure_03.png", 
       plot = figure3, 
       width = 5,    # Width of the plot
       height = 4,   # Height of the plot
       units = "in", # Units: "in", "cm", "mm", or "px"
       dpi = 300)    # Resolution
