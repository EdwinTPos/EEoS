# =====================================================================
# Code for
# A coherent signature of disturbance in successional forests
# Code written and maintained by Edwin Pos and Erica Newman
# E.T.Pos@uu.nl, ericanewman.sci@gmail.com
# Updated October 6, 2026
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

source("13_functions_EoS.R")

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

custom_colors <-     c("CR" =       "#51310B", 
                   "FEB" =      "#8E4F00", 
                   "JE" =       "#DEC27C",
                   "LEP" =      "#999933",
                   "LEPviejo" = "#117733", 
                   "LSUR" =     "#7CCCBE",  
                   "SV" =       "#33BBEE", 
                   "TIR" =      "#332288")

my_shapes <- c("CR" =       16,   # filled circle, secondary
               "FEB" =      16,   # filled circle, secondary
               "JE" =       16,   # filled circle, secondary
               "LEP" =      16,   # filled circle, secondary
               "LEPviejo" = 15,   # filled square -> primary
               "LSUR" =     16,   # filled circle, secondary
               "SV" =       15,   # filled square -> primary
               "TIR" =      16)   # filled circle, secondary
               
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

desired_order <- c("SV", "LEPviejo", "CR", "FEB", "LEP", "LSUR", "TIR" , "JE")

make_plot <- function(plot_name) {

  plot_data <- output_table %>% filter(Plot == plot_name)
  if (nrow(plot_data) == 0) return(NULL)

  forest_type <- unique(plot_data$Foresttype)
  if (length(forest_type) != 1) forest_type <- paste(forest_type, collapse = "/")
  title_label <- paste(plot_name, " (", forest_type, ")", sep = "")

  year_range <- range(plot_data$Year, na.rm = TRUE)

  df_long <- plot_data %>%
    pivot_longer(
      cols      = c(metabol_KL, abundan_KL),
      names_to  = "Metric_Label",
      values_to = "KL_Value"
    ) %>%
    mutate(
      Metric_Label = factor(
        recode(Metric_Label,
               "metabol_KL" = "IPD",
               "abundan_KL" = "SAD"),
        levels = c("IPD", "SAD")
      )
    )

  ggplot(df_long, aes(
    x     = as.numeric(Year),
    y     = as.numeric(KL_Value),
    color = Plot,
    group = interaction(Metric_Label, Category)
  )) +
    geom_smooth(
      aes(linetype = Metric_Label),
      method = "loess", span = 0.75, se = FALSE,
      linewidth = 2.2,          # thicker trend lines
      alpha = 0.8
    ) +
    geom_point(
      aes(shape = Metric_Label),
      size = 3.8,               # larger symbols
      stroke = 1.3,             # thicker symbol strokes
      alpha = 0.7
    ) +
    # hide the per-plot colour legend so only the metric legend is collected
    scale_color_manual(values = custom_colors, drop = FALSE, guide = "none") +
    # same name (NULL) on both -> they merge into ONE legend key
    scale_linetype_manual(name = NULL,
      values = c("IPD" = "solid", "SAD" = "dashed")) +
    scale_shape_manual(name = NULL,
      values = c("IPD" = 8, "SAD" = 3)) +   # 8 = asterisk, 3 = cross (+); use 4 for x
    # make the line/symbol in the legend key large & thick (must match on both to merge)
    guides(
      shape    = guide_legend(override.aes = list(linewidth = 2.2, size = 5, stroke = 1.3)),
      linetype = guide_legend(override.aes = list(linewidth = 2.2, size = 5, stroke = 1.3))
    ) +
    scale_x_continuous(
      breaks = seq(as.numeric(year_range)[1], as.numeric(year_range)[2], by = 5)
    ) +
    ylim(0, .5) +
    # identical titles on every panel -> collected into one shared pair
    labs(title = title_label, x = "Year", y = "KL Divergence") +
    theme_bw() +
    theme(
      plot.title  = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title  = element_text(size = 20),
      axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 18),
      legend.text = element_text(size = 18)
      # note: NO legend.position = "none" here — the legend must exist to be collected
    )
}

plot_list <- lapply(desired_order, make_plot)
plot_list <- plot_list[!sapply(plot_list, is.null)]

composite_plot <- (
  wrap_plots(plot_list, ncol = 4, nrow = 2) +
    plot_layout(guides = "collect", axis_titles = "collect")
) &
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 18),
    legend.key.width = unit(2.2, "cm")   # widen keys so solid vs dashed reads clearly
  )

composite_plot

# Save it
ggsave("Figure_04.png", 
       plot = composite_plot, 
       width = 16,    # Adjust as needed
       height = 10,   # Adjust as needed
       units = "in", 
       dpi = 300)
