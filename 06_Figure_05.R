# =====================================================================
# Code for
# A coherent signature of disturbance in successional forests
# Code written and maintained by Edwin Pos and Erica Newman
# E.T.Pos@uu.nl, ericanewman.sci@gmail.com
# Updated July 13, 2026
# =====================================================================

# ================================================
# FIGURE 5: ORDINATION RESULTS FOR COMPOSITION ANALYSES
# ================================================

#Clear environment, set working directory here and load libraries
#rm(list=ls())
#rm(list = ls(all.names = TRUE))

# Load libraries and dependencies
library(vegan)
library(labdsv)
library(dplyr)
library(tibble)
library(ggplot2)
library(indicspecies)
library(patchwork)

source("12_functions_EoS.R")

# Load full composition data
plantdata <- read.csv("LivingTrees1997-2017.csv", header=T) #Original data from dryad: https://doi.org/10.5061/dryad.ncjsxksvr

#Remove all non-trees
plantdata <- plantdata[plantdata$LifeForm %in% c("tree","midtree","palm","understory palm","Tree","treelet","arbol","Treelet","Palm"),]

#We remove LSUR_new and change BEJ to FEB based on personal communication with R. Chazdon
plantdata <- plantdata[plantdata$PlotName != "LSURnew",]
plantdata[plantdata$PlotName == "BEJ",]$PlotName <- "FEB"
unique(plantdata$PlotName) #check unique names; should not include BEJ or LSUR_new but should include FEB

#Create Count column & unique plot identifier by combining year and plot
plantdata$Count = rep(1,length(rownames(plantdata))) #Set count 
plantdata$PlotCode = paste(plantdata$PlotName,"_",plantdata$Year, sep="")

newdata <- data.frame(plantdata$PlotCode, plantdata$SpeciesCode, plantdata$Count)
colnames(newdata) <- c("PlotCode","SpeciesCode","Count")

newdata %>%
  group_by(PlotCode, SpeciesCode) %>%
  summarise(across(c(Count), sum))

#Now we make the contigency table
vegtab <- matrify(newdata)

#Check if there are any species or sites that have all-zero status
rowSums(vegtab);colSums(vegtab)

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

my_shapes <- c("CR" =       16,   # filled circle, secondary
               "FEB" =      16,   # filled circle, secondary
               "JE" =       16,   # filled circle, secondary
               "LEP" =      16,   # filled circle, secondary
               "LEPviejo" = 15,   # filled square -> primary
               "LSUR" =     16,   # filled circle, secondary
               "SV" =       15,   # filled square -> primary
               "TIR" =      16)   # filled circle, secondary

# ========================================================
# PREPARING THE DATA FOR ANALYSIS
# ========================================================

# Build a per-plot metadata table from the raw records
# Each PlotCode should have consistent PlotName, Foresttype, Year
plot_meta <- plantdata %>%
  group_by(PlotCode) %>%
  summarise(
    PlotName = first(na.omit(PlotName)),
    Foresttype = first(na.omit(Foresttype)),
    Year = first(na.omit(as.numeric(as.character(Year)))),
    .groups = "drop"
  )

# Check for inconsistencies within PlotCodes
duplicates <- plot_meta %>%
  group_by(PlotCode) %>%
  filter(n() > 1)
if (nrow(duplicates) > 0) {
  cat("WARNING: Some PlotCodes have inconsistent metadata. Using first value.\n")
  print(duplicates)
}

# Ensure vegtab row names match PlotCodes in plot_meta
cat("vegtab rows:", nrow(vegtab), "\n")
cat("Unique PlotCodes in plantdata:", nrow(plot_meta), "\n")

# Find the intersection — PlotCodes present in both
common_codes <- intersect(rownames(vegtab), plot_meta$PlotCode)
cat("Common PlotCodes:", length(common_codes), "\n")

# Subset and align both to common codes just in case
vegtab <- vegtab[common_codes, , drop = FALSE]
plot_meta <- plot_meta %>%
  filter(PlotCode %in% common_codes) %>%
  arrange(match(PlotCode, rownames(vegtab)))

# Final alignment check
stopifnot(rownames(vegtab) == plot_meta$PlotCode)

#Split using plot_meta, NOT plantdata
primary_idx   <- which(plot_meta$Foresttype == "Primary")
secondary_idx <- which(plot_meta$Foresttype == "Secondary")

veg_primary   <- vegtab[primary_idx, , drop = FALSE]
veg_secondary <- vegtab[secondary_idx, , drop = FALSE]

meta_primary   <- plot_meta[primary_idx, ]
meta_secondary <- plot_meta[secondary_idx, ]

# ===========================================================
# RUNNING ORDINATION: NMDS
# ===========================================================

set.seed(123) # start of NMDS is stochastic — fix the seed for reproducibility

nmds <- metaMDS(vegtab, distance = "bray", binary=T)

# 1. Extract site scores from the NMDS object
# Assuming 'nmds' is your object and you want Dimension 1 and 2
site_scores <- vegan::scores(nmds, display = "sites") %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "PlotCode")

# 2. Merge with metadata
# Ensure 'newdata' has columns: PlotCode, Plot (Location), Year
plot_data <- left_join(site_scores, plantdata, by = "PlotCode")

# Check for NAs if the merge fails
if(any(is.na(plot_data$Plot))) {
  warning("Some sites in DCA were not found in metadata. Check PlotCode spelling.")
}

# 3. Define factors to ensure correct ordering and matching with your vectors
# It is crucial that the levels of the factor match the names/indices of your color/shape vectors
plot_data$PlotName <- factor(plot_data$PlotName, levels = unique(plot_data$PlotName))

# 4. Create the NMDS plot
nmds_plot <- ggplot(plot_data, aes(x = NMDS1, y = NMDS2, shape = PlotName, color = PlotName)) + 
   geom_point(size = 4, stroke = 1, alpha = 0.7) +
  scale_color_manual(name = "Plotname", values = my_colors) +
  scale_shape_manual(name = "Plotname", values = my_shapes) +
theme_bw(base_size = 15) +
  theme(
    axis.text.x = element_text(size = 20, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.key.size = unit(1.5, "line"),
    legend.position = "right",
    legend.title = element_text(size = 15, face = "plain"),
    plot.title = element_text(size = 20, face = "bold"), 
    plot.subtitle = element_text(size = 15)
  )+
    # Add year labels and add nudge_x/nudge_y to prevent text overlapping the point
  geom_text(aes(label = Year), 
             nudge_y = 0.05, 
            size = 7, 
            check_overlap = TRUE) +
  labs(
    x = "NMDS Axis 1",
    y = "NMDS Axis 2",
    )

ggsave("Figure_05A_NMDS.png", 
       plot = nmds_plot, 
       width = 16, height = 10, units = "in", dpi = 300)

# ===========================================
# INDICATOR SPECIES ANALYSIS OVER TIME
# ===========================================

# Step 1: Identify primary forest indicator species
set.seed(123)      # permutation tests are stochastic — fix the seed for reproducibility

isa_result <- multipatt(
  x       = vegtab,
  cluster = plot_meta$Foresttype,
  func    = "IndVal",
  duleg   = TRUE,                 # 2 groups → only single-group tests exist anyway
  control = how(nperm = 9999)     
)

# Recover classic IndVal and Holm p-values
isa_result$sign$IndVal <- isa_result$sign$stat^2 #multipatt reports the square root of the IndVal
isa_result$sign$p.holm <- p.adjust(isa_result$sign$p.value, method = "holm")

# Attach the A (specificity) and B (fidelity) components for the Primary group
isa_result$sign$A <- isa_result$A[rownames(isa_result$sign), "Primary"]
isa_result$sign$B <- isa_result$B[rownames(isa_result$sign), "Primary"]

# Primary indicators: assigned to Primary, IndVal >= 0.5, Holm p < 0.05
sig_idx <- which(
  isa_result$sign$s.Primary == 1    &
  isa_result$sign$IndVal    >= 0.5  &
  isa_result$sign$p.holm    <  0.05
)

indicator_species <- rownames(isa_result$sign)[sig_idx]

primary_indicators <- isa_result$sign[sig_idx, c("IndVal", "A", "B", "p.holm")]
primary_indicators <- primary_indicators[order(-primary_indicators$IndVal), ]

cat("\nPrimary-forest indicators (IndVal >= 0.5, Holm p < 0.05):",
    length(indicator_species), "\n\n")
print(primary_indicators)

# DIAGNOSTIC: verify they exist in veg_secondary
cat("\nMatching in veg_secondary columns:", 
    sum(indicator_species %in% colnames(veg_secondary)), 
    "of", length(indicator_species), "\n")

# Step 2: Count indicators per secondary plot per year
total_primary_indicators <- length(indicator_species)

indicator_accumulation <- data.frame()

for (i in seq_len(nrow(veg_secondary))) {
  yr <- meta_secondary$Year[i]
  plot_name <- meta_secondary$PlotName[i]
  
  abund_ind <- as.numeric(veg_secondary[i, indicator_species])
  n_indicators <- sum(abund_ind > 0)
  
  indicator_accumulation <- rbind(indicator_accumulation, data.frame(
    PlotName = plot_name,
    Year = yr,
    N_Indicators = n_indicators,
    Prop_Indicators = n_indicators / total_primary_indicators
  ))
}

# Plot the indicator accumulation
isa_accumulation <- ggplot(indicator_accumulation, 
                           aes(x = factor(Year), y = Prop_Indicators, 
                               group = PlotName, color = PlotName, shape=PlotName)) +
  geom_line(aes(color = PlotName), linewidth = 4, alpha = 0.5) +
  geom_point(size = 4, stroke = 1, alpha = 0.7) +
  scale_color_manual(name = "Plot", values = my_colors) +
  scale_shape_manual(name = "Plot", values = my_shapes) +
theme_bw(base_size = 15) +
  theme(
    axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 20),
    legend.key.size = unit(1.5, "line"),
    legend.position = "right",
    legend.title = element_text(size = 15, face = "plain"),
 )+
  labs(y = "Proportion of Primary Forest indicators",
       x = "Year")

ggsave("Figure_05B_ISA.png", 
       plot = isa_accumulation, 
       width = 16, height = 10, units = "in", dpi = 300)

# ========================================================
# BETA-COMP ANALYSIS; TURNOVER, NESTEDNESS AND SIMILARITY
# =========================================================

# Convert to binary matrices explicitly
veg_primary_bin <- ifelse(as.matrix(veg_primary) > 0, 1L, 0L)
veg_secondary_bin <- ifelse(as.matrix(veg_secondary) > 0, 1L, 0L)

# ---- Helper: Sørensen (total) and Simpson (turnover) dissimilarity between two binary vectors
sor_sim <- function(v1, v2) {
  a <- sum(v1 == 1L & v2 == 1L)   # shared
  b <- sum(v1 == 1L & v2 == 0L)   # v1-only
  c <- sum(v2 == 1L & v1 == 0L)   # v2-only
  den_sor <- 2 * a + b + c
  if (den_sor == 0) return(list(sor = NA_real_, sim = NA_real_))  # both plots empty
  den_sim <- a + min(b, c)
  list(
    sor = (b + c) / den_sor,
    sim = if (den_sim > 0) min(b, c) / den_sim else 0    # one plot empty -> pure nestedness
  )
}

# ---- Input checks
stopifnot(
  identical(colnames(veg_primary_bin), colnames(veg_secondary_bin)),  # same species, same order
  !anyNA(veg_primary_bin), !anyNA(veg_secondary_bin),
  nrow(veg_primary_bin)   == nrow(meta_primary),
  nrow(veg_secondary_bin) == nrow(meta_secondary)
)

# ---- Main loop
results <- list()

for (yr in sort(unique(meta_secondary$Year))) {
  prim_idx <- which(meta_primary$Year == yr)
  if (length(prim_idx) < 2) next   # need >= 2 primary plots for a baseline

  # Baseline: mean dissimilarity among primary plots in this year
  pp <- apply(combn(prim_idx, 2), 2, function(pr)
    unlist(sor_sim(veg_primary_bin[pr[1], ], veg_primary_bin[pr[2], ])))
  beta_pp_sor <- mean(pp["sor", ], na.rm = TRUE)
  beta_pp_sim <- mean(pp["sim", ], na.rm = TRUE)
  beta_pp_nes <- beta_pp_sor - beta_pp_sim

  for (i in which(meta_secondary$Year == yr)) {
    # Secondary plot vs each primary plot, averaged
    sp <- sapply(prim_idx, function(p)
      unlist(sor_sim(veg_secondary_bin[i, ], veg_primary_bin[p, ])))
    beta_sp_sor <- mean(sp["sor", ], na.rm = TRUE)
    beta_sp_sim <- mean(sp["sim", ], na.rm = TRUE)
    beta_sp_nes <- beta_sp_sor - beta_sp_sim

    # Excess over the primary-primary baseline, and relative to the remaining room
    excess_sor <- beta_sp_sor - beta_pp_sor #sorensen of secondary vs primary in excess of primary vs. primary
    excess_sim <- beta_sp_sim - beta_pp_sim
    excess_nes <- beta_sp_nes - beta_pp_nes
    rel_sor <- if (1 - beta_pp_sor > 0) excess_sor / (1 - beta_pp_sor) else NA_real_
    rel_sim <- if (1 - beta_pp_sim > 0) excess_sim / (1 - beta_pp_sim) else NA_real_

    results[[length(results) + 1]] <- data.frame(
      PlotName    = meta_secondary$PlotName[i],
      Year        = yr,
      Beta_SP_SOR = beta_sp_sor,
      Beta_PP_SOR = beta_pp_sor,
      Excess_SOR  = excess_sor,
      Rel_SOR     = rel_sor,
      Beta_SP_SIM = beta_sp_sim,
      Beta_PP_SIM = beta_pp_sim,
      Excess_SIM  = excess_sim,
      Rel_SIM     = rel_sim,
      Beta_SP_NES = beta_sp_nes,
      Beta_PP_NES = beta_pp_nes,
      Excess_NES  = excess_nes
    )
  }
}

beta_results <- bind_rows(results)

# ---- Total Sørensen dissimilarity with primary-primary baseline ----
beta_SOR_plot <- ggplot(beta_results, aes(x = factor(Year), y = Beta_SP_SOR,
                          group = PlotName, color = PlotName, shape = PlotName)) +
  geom_line(aes(color = PlotName), linewidth = 4, alpha = 0.5) +
  geom_point(size = 4) +
  geom_step(data = beta_results, inherit.aes = FALSE,
            aes(x = factor(Year), y = Beta_PP_SOR, group = 1),
            linewidth = 1.5, linetype = "longdash", color = "black") +
  scale_color_manual(name = "Plot", values = my_colors) +
  scale_shape_manual(name = "Plot", values = my_shapes) +
theme_bw(base_size = 15) +
  theme(
    axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 20),
    legend.key.size = unit(1.5, "line"),
    legend.position = "right",
    legend.title = element_text(size = 15, face = "plain"),
    plot.title = element_text(size = 20, face = "bold"), 
    plot.subtitle = element_text(size = 15)
  )+
  labs(y = "Beta-SOR (secondary vs primary)",
       x = "Year")

ggsave("Figure_05D_beta_sor_baseline.png", 
       plot = beta_SOR_plot, 
       width = 16, height = 10, units = "in", dpi = 300)


# ---- Turnover component with primary-primary baseline ----
beta_sim_plot <- ggplot(beta_results, aes(x = factor(Year), y = Beta_SP_SIM,
                          group = PlotName, color = PlotName, shape = PlotName)) +
  geom_line(aes(color = PlotName), linewidth = 4, alpha = 0.5) +
  geom_point(size = 4) +
  geom_step(data = beta_results, inherit.aes = FALSE,
            aes(x = factor(Year), y = Beta_PP_SIM, group = 1),
            linewidth = 1.5, linetype = "longdash", color = "black") +
  scale_color_manual(name = "Plot", values = my_colors) +
  scale_shape_manual(name = "Plot", values = my_shapes) +
theme_bw(base_size = 15) +
  theme(
    axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 20),
    legend.key.size = unit(1.5, "line"),
    legend.position = "right",
    legend.title = element_text(size = 15, face = "plain"),
    plot.title = element_text(size = 20, face = "bold"), 
    plot.subtitle = element_text(size = 15)
  )+
  labs(y = "Beta-SIM (secondary vs primary)",
       x = "Year")

ggsave("Figure_05E_beta_sim_baseline.png", 
       plot = beta_sim_plot, 
       width = 16, height = 10, units = "in", dpi = 300)


# ---- Nestedness component with primary-primary baseline ----
nestedness <- ggplot(beta_results, aes(x = factor(Year), y = Beta_SP_NES,
                          group = PlotName, color = PlotName, shape = PlotName)) +
  geom_line(aes(color = PlotName), linewidth = 4, alpha = 0.5) +
  geom_point(size = 4) +
  geom_step(data = beta_results, inherit.aes = FALSE,
            aes(x = factor(Year), y = Beta_PP_NES, group = 1),
            linewidth = 1.5, linetype = "longdash", color = "black") +
  scale_color_manual(name = "Plot", values = my_colors) +
  scale_shape_manual(name = "Plot", values = my_shapes) +
theme_bw(base_size = 15) +
  theme(
    axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 20),
    legend.key.size = unit(1.5, "line"),
    legend.position = "right",
    legend.title = element_text(size = 15, face = "plain"),
    plot.title = element_text(size = 20, face = "bold"), 
    plot.subtitle = element_text(size = 15)
  )+
  labs(y = "Beta-NES (nestedness)",
       x = "Year")

ggsave("Figure_05C_Nestedness_baseline.png", 
       plot = nestedness, 
       width = 16, height = 10, units = "in", dpi = 300)

# ---- Relative Sørensen: 0 = within primary-primary variation ----
beta_rel_plot <- ggplot(beta_results, aes(x = factor(Year), y = Rel_SOR,
                          group = PlotName, color = PlotName, shape = PlotName)) +
  geom_line(aes(color = PlotName), linewidth = 4, alpha = 0.5) +
  geom_point(size = 4) +
  geom_hline(yintercept = 0, linewidth = 1.5, linetype = "longdash") +
    scale_color_manual(name = "Plot", values = my_colors) +
  scale_shape_manual(name = "Plot", values = my_shapes) +
theme_bw(base_size = 15) +
  theme(
    axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 20),
    legend.key.size = unit(1.5, "line"),
    legend.position = "right",
    legend.title = element_text(size = 15, face = "plain"),
    plot.title = element_text(size = 20, face = "bold"), 
    plot.subtitle = element_text(size = 15)
  )+
  labs(y = "Relative Beta-SOR (secondary vs primary)",
       x = "Year")

ggsave("Figure_05F_beta_rel_sor.png", 
       plot = beta_rel_plot, 
       width = 16, height = 10, units = "in", dpi = 300)


# ===================
# Combine all plots
# ===================

# ---- Assemble: 2 x 2 grid, panels labeled A–D in top-left corners ----
# Keep legend only on the bottom-right panel (D), hide everywhere else
combined_figure <- (isa_accumulation + theme(legend.position = "none") + beta_SOR_plot + theme(legend.position = "none")) /
                   (beta_sim_plot + theme(legend.position = "none") + nestedness) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag.position = c(0.02, 0.98),
    plot.tag = element_text(size = 24, face = "bold"),
    plot.title = element_text(hjust = 0.5, size=20),
  plot.subtitle = element_text(hjust = 0.5),
  legend.title = element_text(size = 22, face = "bold"),   # Slightly larger legend title
    legend.text = element_text(size = 18),                 # Slightly larger legend text
    legend.key.size = unit(2, "cm")                      # Slightly larger legend symbols
  )
 
ggsave("Figure_05_Combined_Beta_Analyses.png", 
       plot = combined_figure, 
       width = 24, height = 18, units = "in", dpi = 300)
