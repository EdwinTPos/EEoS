# =====================================================================
# Code for
# A coherent signature of disturbance in successional forests
# Code written and maintained by Edwin Pos and Erica Newman
# E.T.Pos@uu.nl, ericanewman.sci@gmail.com
# Updated July 13, 2026
# =====================================================================

# ================================================
# FIGURE 4: KL divergence RAD AND IPD FOR TREEDATA
# ================================================

#Clear environment, set working directory here and load libraries
#rm(list=ls())
#rm(list = ls(all.names = TRUE))

# Load libraries and dependencies
library(meteR)
library(scales)
library(stringr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(dplyr)

source("12_functions_EoS.R")

# Load data and run function to compare IPD and RAD results between predicted by METE 
# and observed from empirical data
trees <- read.csv("trees.csv", header=T)
treedata <- read.csv("treedata.csv", header=T)

#Save all unique samples to iterate over and indicate which metric to use:
#Can be "AIC", "LL" (loglikelihood), "LLZ" (z-score transformed loglikelihood) or "KL" (Kullback-Leibler)
Samples <- unique(treedata$PlotCode)
metric <- ("KL") 

#Compare the IPD and SAD
output_table <- METE_comparison(Samples, metric, treedata)

write.csv(output_table, "output_table.csv", row.names = FALSE)

#When continuing session and file is provided
#output_table <- read.csv("output_table.csv", header=T) 

# Custom color palette and custom shapes

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

#Set category (combination of plot and foresttype) for nice plotting
output_table$Category = paste(output_table$Foresttype,","," ",
                              output_table$Plot)


# ===========================================================
# OPTION 1: two graphs (IPD and SAD) with all plots combines
# ===========================================================

ipd.graph <-
ggplot(output_table, aes(
  x = Year,
  y = metabol_KL,
  group = Category,
  linetype = Foresttype,
  color = Plot
)) +
  # Smooth model (hidden)
  geom_smooth(
    method = "loess",
    span = 0.75,
    se = FALSE,
    linewidth = 0,   
    fullrange = TRUE
  ) +
  # Overlay semi-transparent fitted lines
  stat_smooth(
    method = "loess",
    span = 0.75,
    se = FALSE,
    geom = "line",
    linewidth = 1.2,
    alpha = 0.8     
  ) +
  # applying colors 
  scale_color_manual(values = custom_colors) +
  labs(
    x = "Year",
    y = "Kullback-Liebler divergence"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.key.size = unit(2, "line")
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1)))

sad.graph <-
  ggplot(output_table, aes(
  Year,
  abundan_KL,
  group = Category,
  linetype = Foresttype,
  color = Plot
)) +
  # Smooth model, hidden line (only compute data)
  geom_smooth(
    aes(color = Plot),
    method = "loess",
    span = 0.75,
    se = FALSE,
    linewidth = 0,   # hide the internal line
    fullrange = TRUE
  ) +
  # Overlay semi-transparent fitted lines manually
  stat_smooth(
    aes(color = Plot),
    method = "loess",
    span = 0.75,
    se = FALSE,
    geom = "line",
    linewidth = 1.5,
    alpha = 0.8      # <— transparency for line
  ) +
  scale_color_manual(values = custom_colors) +
  xlab("Year") +
  ylab("Kullback-Liebler divergence") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.key.size = unit(2, "line")
  )

#ipd.graph
#sad.graph

# ===============================================================================
# OPTION 2: separate figures for each location -> this was used in the manuscript
# ===============================================================================

desired_order <- c("SV", "LEPviejo", "CR", "FEB", "LEP", "TIR", "LSUR", "JE")

make_plot <- function(plot_name) {

  plot_data <- output_table %>% filter(Plot == plot_name)
  if (nrow(plot_data) == 0) return(NULL)
  
  year_range <- range(plot_data$Year, na.rm = TRUE)
  
df_long <- plot_data %>%
  pivot_longer(
    cols = c(metabol_KL, abundan_KL),
    names_to = "Metric_Label",
    values_to = "KL_Value"
  ) %>%
  mutate(
    Metric_Label = factor(
      recode(Metric_Label,
             "metabol_KL" = "IPD",
             "abundan_KL" = "SAD"),
      levels = c("IPD", "SAD") # Explicitly set the order here
    )
  )

  ggplot(df_long, aes(
    x = as.numeric(Year),
    y = as.numeric(KL_Value),
    color = Plot,  # Color now controls the lines and the open shape outlines
    group = interaction(Metric_Label, Category)
  )) +
    # Wide lines for trends
    geom_smooth(
      aes(linetype = Metric_Label),
      method = "loess", 
      span = 0.75,
      se = FALSE,
      linewidth = 1.8, 
      alpha = 0.8
    ) +
    # Open points (0 = square, 1 = circle)
    geom_point(
      aes(shape = Metric_Label), 
      size = 2, 
      stroke = 1, # Controls the thickness of the open shape borders
      alpha = 0.6
    ) + 
    scale_color_manual(values = custom_colors, drop = FALSE) +
    # IPD (metabol_KL) should be solid, SAD (abundan_KL) should be dashed
    scale_linetype_manual(values = c("IPD" = "solid", "SAD" = "dashed")) +
    scale_shape_manual(values = c("IPD" = 0, "SAD" = 1)) + # 0: Open Square, 1: Open Circle
    scale_x_continuous(breaks = seq(as.numeric(year_range)[1],as.numeric(year_range)[2], by = 5)) +
    ylim(0,.5) +
    labs(title = plot_name, x = NULL, y = NULL) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 18),
      legend.title = element_text(size = 12, face = "plain") 
    )
}

# Generate and combine as before
plot_list <- lapply(desired_order, make_plot)
plot_list <- plot_list[!sapply(plot_list, is.null)]

# Build the composite plot
composite_plot <- wrap_plots(plot_list, ncol = 4, nrow = 2) +
  plot_annotation(
    title = "KL Divergence: IPD (Solid/Open Square) vs SAD (Dotted/Open Circle)"
  )

# Save it
ggsave("Figure_04.png", 
       plot = composite_plot, 
       width = 16,    # Adjust as needed
       height = 10,   # Adjust as needed
       units = "in", 
       dpi = 300)
