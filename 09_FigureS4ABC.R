# =====================================================================
# Code for
# A coherent signature of disturbance in successional forests
# Code written and maintained by Edwin Pos and Erica Newman
# E.T.Pos@uu.nl, ericanewman.sci@gmail.com
# Updated October 6, 2026
# =====================================================================

# ==================================================
# FIGURE S4 Error propagation of biomass calculation
# ==================================================

#Clear environment
#rm(list=ls())
#rm(list = ls(all.names = TRUE))

#Load libraries and dependencies
library(dplyr)
library(purrr)
library(BIOMASS)
library(ggplot2)
library(cowplot)

# ================================================
# Error propagation analysis of biomass estimation
# ================================================

plantdata <- read.csv("plantdata_for_error_propagation.csv", header=T)

# Residual standard error for the height model, on the log scale.
# Value (0.243) comes from Chave et al. 2014 (Global Change Biology 20:3177-3190),
# equation 6a: ln(H) = 0.893 - E + 0.760*ln(D) - 0.0340*[ln(D)]^2, fitted on
# 4004 trees across 58 pantropical sites (RSE = 0.243, df = 3998).
# This is the same height model used internally by retrieveH() when called with
# 'coord', so 0.243 is the chosen uncertainty to propagate through AGBmonteCarlo().
errH_value <- 0.243

# Monte Carlo error propagation per plot × year, adapted from vignette of package biomass
run_mc <- function(df, n = 1000) {
  ok <- complete.cases(as.numeric(df$DBH), as.numeric(df$WD_used), as.numeric(df$WD_sd), df$H) &
        df$DBH > 0 & df$WD_used > 0
  if (sum(ok) < 2) return(NULL)

  mc <- AGBmonteCarlo(
    D       = as.numeric(df$DBH)[ok],
    WD      = as.numeric(df$WD_used[ok]),
    errWD   = as.numeric(df$WD_sd[ok]),
    H       = as.numeric(df$H[ok]),
    errH    = errH_value,
    Dpropag = "chave2004",
    n       = n
  )
  data.frame(
    PlotName  = unique(df$PlotName),
    Year      = unique(df$Year),
    n_trees   = sum(ok),
    AGB_mean  = mc$meanAGB,
    AGB_med   = mc$medAGB,
    AGB_sd    = mc$sdAGB,
    AGB_lwr   = mc$credibilityAGB[1],
    AGB_upr   = mc$credibilityAGB[2],
    CV_pct    = 100 * mc$sdAGB / mc$meanAGB
  )
}

groups <- plantdata |>
  group_by(PlotName, Year) |>
  group_split()

mc_results <- groups |>
  map(run_mc, .progress = TRUE) |>
  bind_rows()

#Match AGB ouput in KG
mc_results <- mc_results |>
  mutate(across(c(AGB_mean, AGB_med, AGB_sd, AGB_lwr, AGB_upr),
                ~ .x * 1000))

save(mc_results, file = "mc_results.Rdata")

#Load data if continuing from previous run
load("mc_results.Rdata")

biomass_precision <- mc_results %>%
  mutate(Year = suppressWarnings(as.integer(Year))) %>%
  group_by(PlotName) %>%
  summarise(
    n_years        = n(),
    years          = paste0(min(Year), "-", max(Year)),
    n_trees_mean   = mean(n_trees,  na.rm = TRUE),  # typical plot size
    AGB_mean       = mean(AGB_mean, na.rm = TRUE),  # typical estimate
    AGB_SD         = mean(AGB_sd,   na.rm = TRUE),  # typical bootstrap SD
    AGB_lwr        = mean(AGB_lwr,  na.rm = TRUE),  # typical 2.5% bound
    AGB_upr        = mean(AGB_upr,  na.rm = TRUE),  # typical 97.5% bound
    CV_pct         = mean(CV_pct,   na.rm = TRUE),  # typical precision
    CV_min_pct     = min(CV_pct,    na.rm = TRUE),  # best  year
    CV_max_pct     = max(CV_pct,    na.rm = TRUE),  # worst year
    .groups = "drop"
  ) %>%
  arrange(PlotName)
 
# rounded version for potential use in manuscript
biomass_precision_disp <- biomass_precision %>%
  mutate(
    n_trees_mean = round(n_trees_mean, 0),
    across(c(AGB_mean, AGB_SD, AGB_lwr, AGB_upr), ~round(.x, 0)),
    across(ends_with("_pct"), ~round(.x, 1))
  )

print(biomass_precision_disp)

#Table of precision results -> not used for the manuscript
write.csv(biomass_precision_disp, "table_S6_biomass_precision.csv", row.names = FALSE) 
#biomass_precision_disp <- read.csv("table_S6_biomass_precision.csv")
 
# ===================================
# Visualizations of error propagation
# ===================================

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

# =======================================================
# FIGURE S4A: Plot-level AGB with 95% CI, faceted by plot
# =======================================================
Plot_levelAGB_certainty <- ggplot(mc_results, aes(x = Year, y = log(AGB_mean), 
                                                 group = PlotName, color = PlotName, 
                                                 fill = PlotName, shape = as.character(PlotName))) +
  geom_ribbon(aes(ymin = log(AGB_lwr), ymax = log(AGB_upr)), alpha = 0.15, color = NA) +
  geom_line() +
  geom_pointrange(aes(ymin = log(AGB_lwr), ymax = log(AGB_upr)), size = 0.25) +
  facet_wrap(~ PlotName, scales = "free_y") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(x = "Year", y = "Plot log(B) ± 95% CI", title = "Plot-level B with propagated uncertainty") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none",
        legend.title = element_text(size = 9, face = "plain")) 

# 2. Add text in space for "9th" panel so the figure is a nice square
Final_Plot <- ggdraw(Plot_levelAGB_certainty) +
  draw_label(
    label = " ", 
    x = 0.84,       # Horizontal position (0 = left, 1 = right)
    y = 0.20,       # Vertical position (0 = bottom, 1 = top)
    hjust = 0.5,    # Centre text horizontally
    vjust = 0.5,    # Centre text vertically
    size = 10,      # Font size
    fontface = "italic", # Style ("plain", "bold", "italic")
    color = "black" # Color
  )

# 3. View and Save plot
Final_Plot

ggsave("figure_S4A.png", 
       plot = Final_Plot, 
       width = 5,    # Width of the plot
       height = 3.5,   # Height of the plot
       units = "in", # Units: "in", "cm", "mm", or "px"
       dpi = 300)    # Resolution

# ====================================================
# Figure S4B: coefficient of variation per plot × year
# ====================================================

Error_CV_biomass <- ggplot(mc_results, aes(x = Year, y = CV_pct,
                       color = PlotName, 
                       fill = PlotName, 
                       shape = as.character(PlotName))) +
    # Fix: Added group = PlotName so R knows how to connect the lines
    geom_line(aes(group = PlotName), linewidth = 1, alpha = 0.5) + 
    geom_point(size = 2, stroke = 0.5, alpha = .8) +
  scale_color_manual(name = "Plot", values = my_colors) +
  scale_fill_manual(name = "Plot", values = my_colors) + 
  scale_shape_manual(name = "Plot", values = my_shapes) +
  labs(x = "Year", 
       y = "CV of AGB estimate (%)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "right",
    legend.title = element_text(size = 12, face = "plain"))+
  theme_bw(base_size = 11)

# Save the plot
ggsave("figure_S4B.png", 
       plot = Error_CV_biomass, 
       width = 5,    
       height = 3.5,   
       units = "in", 
       dpi = 300)

# ======================================================
# Figure S4C: Forest plot of most recent census per plot
# ======================================================

latest <- mc_results |>
  group_by(PlotName) |>
  filter(Year == max(Year)) |>
  ungroup()

range_recent_ABG <- ggplot(latest, aes(x = reorder(PlotName, log(AGB_mean)), y = log(AGB_mean))) +
  geom_pointrange(aes(ymin = log(AGB_lwr), ymax = log(AGB_upr)),
                  color = "steelblue") +
  coord_flip() +
  labs(x = NULL, y = "Plot log(B) (Kg), 95% CI",
       title = "Most recent census per plot")+
  theme_bw(base_size = 11) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "plain") 
  )

  ggsave("figure_S4C.png", 
       plot = range_recent_ABG, 
       width = 5,    # Width of the plot
       height = 3.5,   # Height of the plot
       units = "in", # Units: "in", "cm", "mm", or "px"
       dpi = 300)    # Resolution
