# =====================================================================
# Code for
# A coherent signature of disturbance in successional forests
# Code written and maintained by Edwin Pos and Erica Newman
# E.T.Pos@uu.nl, ericanewman.sci@gmail.com
# Updated July 13, 2026
# =====================================================================

#Clear environment
#rm(list=ls())
#rm(list = ls(all.names = TRUE))

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

my_shapes <- c("CR" =       21, 
               "FEB" =      8, 
               "JE" =        24, 
               "LEP" =      25, 
               "LEPviejo" = 12, 
               "LSUR" =     23, 
               "SV" =        22, 
               "TIR" =       3)

# Reorder the 'Plot' factor levels based on the (StatVar) column
# (using -value to get descending order) highest to lowest

trees$Plot <- reorder(trees$Plot, -trees[["dist"]], FUN = max, na.rm = TRUE)

figure2B <- ggplot(trees, aes(Year, abs(.data[["dist"]]), color = Plot, shape = Plot)) +
  geom_line(aes(color = Plot), linewidth = 2, alpha = 0.5) +
  geom_point(size = 1, stroke = 1, alpha = 0.7) +  
  # To merge into ONE legend, 'name' must be exactly the same in both scales
  scale_color_manual(name = "Plot", values = my_colors) +
  scale_shape_manual(name = "Plot", values = my_shapes) +
  xlab("Year") +
  #ylab(paste("Ln(", StatVar, ")", sep = "")) +
  ylab("Predicted - observed ln(B)") +
  guides(
    shape = guide_legend(override.aes = list(size = 3)),
    fill = guide_legend(override.aes = list(shape = 21)) # Forces circles in the fill legend
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = .5, vjust = 0.5, size = 12),
    axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.key.size = unit(1.5, "line"),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "plain")
  )

#figure2B

#Save figure
ggsave("Figure_02B.png", 
       plot = figure2B, 
       width = 5,    # Width of the plot
       height = 4,   # Height of the plot
       units = "in", # Units: "in", "cm", "mm", or "px"
       dpi = 300)    # Resolution
