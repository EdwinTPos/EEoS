# =====================================================================
# Code for
# A coherent signature of disturbance in successional forests
# Code written and maintained by Edwin Pos and Erica Newman
# E.T.Pos@uu.nl, ericanewman.sci@gmail.com
# Updated October 6, 2026
# =====================================================================

# ======================================
# FIGURE S2 biomass of each plot by year
# ======================================

#Clear environment
#rm(list=ls())
#rm(list = ls(all.names = TRUE))

#install.packages("ggraph")

#Load libraries and dependencies
library(ggplot2)
library(ggraph)
library(ggrepel)
library(patchwork)
library(dplyr)

# Load the data
trees <- read.csv("trees.csv", header=T)

# Custom color palette
my_colors <-     c("CR" =       "#51310B", 
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

# 1. Define specific order of calculation
desired_order <- c("SV", "LEPviejo", "CR", "FEB", "LEP", "TIR", "LSUR", "JE")

# 2. Convert Plot to a factor to enforce the order
trees$Plot <- factor(trees$Plot, levels = desired_order)

# 3. Create the subsets for the two subfigures A and B based on the new plot factor levels
trees_subset1 <- trees[trees$Plot %in% c("SV", "LEPviejo", "CR", "FEB"),]
trees_subset2 <- trees[trees$Plot %in% c("LEP", "TIR", "LSUR", "JE"),]

# 4. Define general plotting function
make_square_plot <- function(sub_df) {
  lx <- log(sub_df$pBnum_exact)
  ly <- log(sub_df$B)
  raw_range <- range(c(lx, ly), na.rm = TRUE)
  padded_range <- c(raw_range[1] - 0.1, raw_range[2] + 0.1)

  ggplot(sub_df, aes(x = log(pBnum_exact), y = log(B),
                     color = Plot, fill = Plot, shape = Plot,
                     alpha = Year)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                colour = "#595959", alpha = 0.8, lwd = 0.8) +
    geom_point(size = 2.5, stroke = 1.2) +
    coord_fixed(ratio = 1, xlim = raw_range, ylim = raw_range) +
    scale_color_manual(values = my_colors, guide = "none") +
    scale_fill_manual(values = my_colors, guide = "none") +
    scale_shape_manual(values = my_shapes, guide = "none") +
    scale_alpha_continuous(
      range = c(0.25, 1),
      breaks = scales::breaks_pretty(n = 4),
      guide = guide_legend(
        title = "Year",
        override.aes = list(
          color = unique(sub_df$Plot) %>% as.character() %>% (\(x) my_colors[x])(),
          fill  = unique(sub_df$Plot) %>% as.character() %>% (\(x) my_colors[x])(),
          shape = unique(sub_df$Plot) %>% as.character() %>% (\(x) my_shapes[x])(),
          size  = 2.5,
          stroke = 1.2
        )
      )
    ) +
    labs(title = unique(sub_df$Plot), x = NULL, y = NULL) +
    theme_bw() +
    theme(
      legend.position = c(0.22, 0.61),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
      legend.key = element_blank(),
      legend.text = element_text(size = 5),
      legend.key.size = unit(0.4, "cm"),
      legend.margin = margin(2, 4, 2, 4),
      aspect.ratio = 1,
      legend.title = element_text(size = 12, face = "plain") 
    )
}

# 5. Create the first 2x2 graphic (Subset 1: SV, LEPviejo, CR, FEB)
# group_split respects factor levels, so this will be in order 1-4
list_1 <- trees_subset1 %>%
  group_split(Plot) %>%
  lapply(make_square_plot)

fig_s2_a <- wrap_plots(list_1, ncol = 2, nrow = 2) +
  plot_annotation(title = "Biomass Predictions: Group A")

# 6. Create the second 2x2 graphic (Subset 2: LEP, TIR, LSUR, JE)
# This will be in order 5-8
list_2 <- trees_subset2 %>%
  group_split(Plot) %>%
  lapply(make_square_plot)

fig_s2_b <- wrap_plots(list_2, ncol = 2, nrow = 2) +
  plot_annotation(title = "Biomass Predictions: Group B")

# 7. View and Save
#fig_s2_a
#fig_s2_b

ggsave("figure_S2A.png", plot = fig_s2_a, width = 7.5, height = 6.5, dpi = 300)
ggsave("figure_S2B.png", plot = fig_s2_b, width = 7.5, height = 6.5, dpi = 300)
