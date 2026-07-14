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

my_shapes <- c("CR" =       21, 
               "FEB" =      8, 
               "JE" =        24, 
               "LEP" =      25, 
               "LEPviejo" = 12, 
               "LSUR" =     23, 
               "SV" =        22, 
               "TIR" =       3)

#Choose which state variable to plot; manuscript shows "B" for Figure 03
StatVar = "B"
#StatVar = "S"
#StatVar = "N"
#StatVar = "E"

# Reorder the 'Plot' factor levels based on the 'B' (StatVar) column
# (using -value to get descending order) highest to lowest

trees$Plot <- reorder(trees$Plot, -trees[[StatVar]], FUN = max, na.rm = TRUE)

figure3 <- ggplot(trees, aes(Year, log(.data[[StatVar]]), label = Category)) +
  geom_line(aes(color = Plot), linewidth = 2, alpha = 0.5) +
  geom_point(aes(fill = Plot, shape = Foresttype), alpha = 0.75, size = 2) +  
  scale_shape_manual(
    name = "Forest Type", # Giving this a name helps distinguish the legends
    values = c("Primary" = 22, "Secondary" = 21)
  ) +
  scale_fill_manual(
    name = "Plot Name", # Your Plot legend title
    values = custom_colors,
    # Note: If using custom labels, they must align with the new factor order
    # If custom_colors is a named vector, it will map correctly automatically.
    # labels = ... 
  ) +
  scale_color_manual(
    values = custom_colors,
    guide = "none" 
  ) +
  xlab("Year") +
  ylab(paste("Ln(", StatVar, ")", sep = "")) +  #general
  #ylab("Ln(observed biomass)") +               #used for manuscript
  guides(
    shape = guide_legend(override.aes = list(size = 3)),
    fill = guide_legend(override.aes = list(shape = 21)) # Forces circles in the fill legend
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = .5, vjust = 0.5, size = 12),
    axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.key.size = unit(1.5, "line"),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "plain") 
  )

  theme(
    
    
  )

#figure3

#Save figure
ggsave("Figure_03.png", 
       plot = figure3, 
       width = 5,    # Width of the plot
       height = 4,   # Height of the plot
       units = "in", # Units: "in", "cm", "mm", or "px"
       dpi = 300)    # Resolution
