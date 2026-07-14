# =====================================================================
# Code for
# A coherent signature of disturbance in successional forests
# Code written and maintained by Edwin Pos and Erica Newman
# E.T.Pos@uu.nl, ericanewman.sci@gmail.com
# Updated July 13, 2026
# =====================================================================

# ================================================================================
# Figure S5: quadrat accumulation curve and precision of the plot-level B estimate
# ================================================================================

#Clear environment, set working directory here and load libraries
#rm(list=ls())
#rm(list = ls(all.names = TRUE))

#Loading libraries and dependencies
library(dplyr)
library(ggplot2)

source("12_functions_EoS.R")

# Set column in the data holding the forest-location codes
location_col <- "PlotName"         

# Custom colour palette and shapes
my_colors <- c("CR"       = "#51310B",
               "FEB"      = "#8E4F00",
               "JE"       = "#DEC27C",
               "LEP"      = "#999933",
               "LEPviejo" = "#117733",
               "LSUR"     = "#7CCCBE",
               "SV"       = "#33BBEE",
               "TIR"      = "#332288")
my_shapes <- c("CR"       = 21,
               "FEB"      = 8,
               "JE"       = 24,
               "LEP"      = 25,
               "LEPviejo" = 12,
               "LSUR"     = 23,
               "SV"       = 22,
               "TIR"      = 3)

# Set warning if a location code is missing from the palette
palette_check <- function(x) {
  miss <- setdiff(unique(as.character(x)), names(my_colors))
  if (length(miss) > 0) {
    warning("Locations not found in the palette (check `location_col`): ",
            paste(miss, collapse = ", "))
  }
}

# ========================================================
# FIGURE S5A - quadrat accumulation curve for B
# Bootstrap is repeated over a gradient of quadrat counts
# ========================================================

# S5A accumulation analysis settings
acc_quad_grid <- c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
acc_n_boot    <- 500       # iterations per quadrat level
acc_seed      <- 20240101  # fixed seed -> reproducible S5A

treedata <- read.csv("treedata.csv", header=T)

# Run accumulation curve, first check if data is available, set seed then run function itself
if (!exists("treedata")) {
  stop("`treedata` not found - needed for the S5A accumulation analysis.")
}

set.seed(acc_seed)
acc <- accumulation_curve_B(
  data      = treedata,
  quad_grid = acc_quad_grid,
  n_boot    = acc_n_boot,
  normalize = FALSE        
)

#Check if all locations are available and set levels
palette_check(acc[[location_col]])
acc$location <- factor(as.character(acc[[location_col]]),
                       levels = names(my_colors))

#Plot figure 5A
p_S5A <- ggplot(acc, aes(x = n_quadrats, y = log(B_estimate),
                         colour = location, fill = location,
                         group = interaction(PlotCode, Year))) +
  geom_ribbon(aes(ymin = log(B_lower), ymax = log(B_upper)),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(aes(shape = location), size = 1, stroke = 0.7) +
  scale_colour_manual(values = my_colors, name = "Location", drop = FALSE) +
  scale_fill_manual(values = my_colors, name = "Location", drop = FALSE) +
  scale_shape_manual(values = my_shapes, name = "Location", drop = FALSE) +
  scale_x_continuous(breaks = acc_quad_grid) +
  labs(
    x = "Number of quadrats sampled",
    y = "Estimated plot total biomass, ln(B) (full-plot scale)",
    title = "Figure S5A. Sampling sufficiency of plot-level total biomass",
    subtitle = "Bootstrap estimate (line) and 95% CI (ribbon) vs. sampling intensity"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, colour = "grey30"),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
    axis.text.x = element_text(size = 8),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "plain")
  )

# All eight overlaid ribbons was very crowded, use a facet instead for separate plots
p_S5A <- p_S5A + facet_wrap(~ location, scales = "free_y")

#NOTE: it can be that there is a warning that in log(B_lower) there are NaN's produced. This happens 
#when B_estimate < 1.96×B_est_se, caused either by the estimate itself being small in magnitude 
# (close to zero) which is unlikely here or the standard error being large relative to the estimate
# indicating substantial uncertainty or imprecision in the estimate, which is likely the case here.
# July 13 this was the case for 9 out of 1496 observations, ± .6%

ggsave("figure_S5A.png", plot = p_S5A,
       width = 8.5, height = 5.5, dpi = 300, bg = "white")

# ===============================================================
# FIGURE S5B - precision of the plot-level B estimate
# CV = 100 * B_se / B_mean ; lower CV = more representative
# ===============================================================

# Retrieve main bootstrap output 
boot_summary <- read.csv("bootstrap_summary.csv", header=T)

#Run palette_check to see if all locations are in the data
palette_check(boot_summary[[location_col]])

#Use the boot_summary for CV calculations
prec <- boot_summary %>%
  mutate(
    CV       = 100 * B_se / B_mean,
    location = factor(as.character(.data[[location_col]]),
                      levels = names(my_colors)),
    plot_id  = PlotCode
  ) %>%
  arrange(CV) %>%
  mutate(plot_id = factor(plot_id, levels = plot_id))

median_cv <- median(prec$CV, na.rm = TRUE)

#Plot figure 5B
p_S5B <- ggplot(prec, aes(x = plot_id, y = CV,
                          colour = location, fill = location)) +
  geom_hline(yintercept = median_cv, linetype = "dashed",
             colour = "grey55", linewidth = 0.4) +
  geom_segment(aes(xend = plot_id, y = 0, yend = CV),
               linewidth = 0.6) +
  geom_point(aes(shape = location), size = 2, stroke = 0.7) +
  coord_flip() +
  scale_colour_manual(values = my_colors, name = "Location", drop = FALSE) +
  scale_fill_manual(values = my_colors, name = "Location", drop = FALSE) +
  scale_shape_manual(values = my_shapes, name = "Location", drop = FALSE) +
  annotate("text", x = 0.7, y = median_cv,
           label = sprintf("median = %.1f%%", median_cv),
           hjust = -0.05, vjust = 0, size = 3, colour = "grey45") +
  labs(
    x = NULL,
    y = "Bootstrap coefficient of variation of B (%)",
    title = "Figure S5B. Precision of single-plot biomass estimates",
    subtitle = "Within-plot bootstrap CV per plot; lower = more representative"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, colour = "grey30"),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
    axis.text.y = element_text(size = 5),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "plain")
  )

#View and save figure
#p_S5B

ggsave("figure_S5B.png", plot = p_S5B,
       width = 7.5, height = 9, dpi = 300, bg = "white")
