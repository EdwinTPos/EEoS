# =====================================================================
# Code for
# A coherent signature of disturbance in successional forests
# Code written and maintained by Edwin Pos and Erica Newman
# E.T.Pos@uu.nl, ericanewman.sci@gmail.com
# Updated October 6, 2026
# =====================================================================

# =================================================================================
# FIGURE S6 Rate of change of state variables S, N, E against KL divergence by plot
# =================================================================================

#Clear environment
#rm(list=ls())
#rm(list = ls(all.names = TRUE))

#Load libraries and dependencies
library(ggplot2)
library(patchwork)
library(dplyr)
library(tibble) 

#Load data
trees <- read.csv("trees.csv", header=T) 
output_table <- read.csv("output_table.csv", header=T) # made in file "05_Figure_04.R"

for(i in 1:length(trees$PlotCode)){
sample = trees[i,"PlotCode"]
trees[i,"IPD KL"] <-  output_table[output_table$Samples == sample,"metabol_KL"]
trees[i,"SAD KL"] <-  output_table[output_table$Samples == sample,"abundan_KL"]
}

#Set custom colors
custom_colors <- c("CR"="#51310B","FEB"="#8E4F00","JE"="#DEC27C","LEP"="#999933",
                   "LEPviejo"="#117733","LSUR"="#7CCCBE","SV"="#33BBEE","TIR"="#332288")
my_shapes <- c("CR"=21,"FEB"=8,"JE"=24,"LEP"=25,"LEPviejo"=12,
               "LSUR"=23,"SV"=22,"TIR"=3)

# Compute growth rates for each state variable, uses function "lag" to identify previous value
df_deriv <- trees %>%
  group_by(Plot) %>%
  arrange(Year) %>%
  mutate(
    dB_dt = (B - lag(B)) / (Year - lag(Year)),
    dN_dt = (N - lag(N)) / (Year - lag(Year)),
    dS_dt = (S - lag(S)) / (Year - lag(Year))
  ) %>%
  filter(!is.na(dB_dt) & !is.na(dN_dt) & !is.na(dS_dt)) %>%
  ungroup()

# Set shared theme for plotting
shared_theme <- theme_bw() +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title   = element_text(size = 12),
    plot.margin = margin(10, 10, 10, 10),
    plot.subtitle = element_text(size = 9, hjust = 0, colour = "grey25",
                                 lineheight = 0.9,
                                 margin = margin(b = 6)),
    legend.position = "right"
  )

# Define panel specifications
model_list <- list(
  list(outcome = "dB_dt", KL = "IPD KL",
       title = "Biomass (B)",
       xlab = expression(paste("dB/dt (yr"^-1, ")"))),
  list(outcome = "dB_dt", KL = "SAD KL",
       title = "Biomass (B)",
       xlab = expression(paste("dB/dt (yr"^-1, ")"))),
  list(outcome = "dN_dt", KL = "IPD KL",
       title = "Individuals (N)",
       xlab = expression(paste("dN/dt (yr"^-1, ")"))),
  list(outcome = "dN_dt", KL = "SAD KL",
       title = "Individuals (N)",
       xlab = expression(paste("dN/dt (yr"^-1, ")"))),
  list(outcome = "dS_dt", KL = "IPD KL",
       title = "Species Richness (S)",
       xlab = expression(paste("dS/dt (yr"^-1, ")"))),
  list(outcome = "dS_dt", KL = "SAD KL",
       title = "Species Richness (S)",
       xlab = expression(paste("dS/dt (yr"^-1, ")")))
)

plot_order <- names(custom_colors)

#Build 6 panels
plots_list <- vector("list", length(model_list))

for (i in seq_along(model_list)) {
  spec <- model_list[[i]]
  outcome_var <- spec$outcome
  KL_var <- spec$KL

  # Per-plot r-values for subtitle
  corr_plot <- df_deriv %>%
    group_by(Plot) %>%
    summarise(
      r = cor(.data[[outcome_var]], .data[[KL_var]], use = "complete.obs",
              method = "pearson"),
      .groups = "drop"
    ) %>%
    mutate(Plot = factor(Plot, levels = plot_order)) %>%
    arrange(Plot)

  line1 <- paste(sprintf("%s: r=%.2f", corr_plot$Plot[1:3], corr_plot$r[1:3]),
                 collapse = "  |  ")
  line2 <- paste(sprintf("%s: r=%.2f", corr_plot$Plot[4:6], corr_plot$r[4:6]),
                 collapse = "  |  ")
  line3 <- paste(sprintf("%s: r=%.2f", corr_plot$Plot[7:8], corr_plot$r[7:8]),
                 collapse = "  |  ")
  subtitle_str <- paste0(line1, "\n", line2, "\n", line3)

  # Subset: primary forest only (for the overall primary LOESS)
  df_primary <- filter(df_deriv, Foresttype == "Primary")

  p <- ggplot(df_deriv, aes(x = .data[[outcome_var]], y = .data[[KL_var]] )) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +

    # Points first (background layer)
    geom_point(
      aes(colour = Plot, fill = Plot, shape = Plot),
      size = 2, alpha = 0.5
    ) +

    # Per-plot LOESS curves (coloured by plot)
    geom_smooth(
      aes(colour = Plot, group = Plot),
      method = "loess", se = FALSE,
      formula = y ~ x, span = 0.75,
      linewidth = 1.0, alpha = 0.8
    ) +

    scale_color_manual(values = custom_colors, name = "Plot") +
    scale_fill_manual(values = custom_colors, name = "Plot") +
    scale_shape_manual(values = my_shapes, name = "Plot") +

    labs(
      title = spec$title,
      subtitle = subtitle_str,
      y = gsub("_", " ", KL_var),
      x = spec$xlab
    ) +
    shared_theme

  plots_list[[i]] <- p
}

# Arrange with single shared legend
fig_s6 <- wrap_plots(plots_list, ncol = 2, nrow = 3, guides = "collect") +
  plot_annotation(tag_levels = 'A') &
  theme(legend.position = "right") &
  guides(
    colour = guide_legend(override.aes = list(alpha = 1, size = 3)),
    fill   = guide_none(),
    shape  = guide_none()
  )

# View and Save
#fig_s6
ggsave("figure_S6.png", plot = fig_s6, width = 10, height = 10, dpi = 300)
