########################
# Code for
# A coherent signature of disturbance in successional forests
# Code written and maintained by Edwin Pos and Erica Newman
# E.T.Pos@uu.nl, ericanewman.sci@gmail.com
# Updated July 13, 2026
########################

#Clear environment
rm(list=ls())
rm(list = ls(all.names = TRUE))

#install.packages(c("BIOMASS","meteR","labdsv","philentropy", "factoextra",
#                 "vegan","gridExtra", "untb", "lmtest", "CausalImpact",
#                 "gifski","ggplot2","reshape2","cowplot","maps",
#                 "tidyverse","ggrepel","sf","rnaturalearth", "Rmpfr",
#                 "bsts","dplyr"
#                 ))

#Load libraries and dependencies
library(BIOMASS)
library(labdsv)
library(meteR)
library(philentropy)
library(factoextra)
library(vegan)
library(gridExtra)
library(untb)
library(lmtest)
library(ggplot2)
library(reshape2)
library(maps)
library(tidyverse)
library(ggrepel)
library(sf)
library(Rmpfr)
library(dplyr)
library(purrr)

source("12_functions_EoS.R")

# ============
# LOADING DATA
# ============

#Load prepared dataset and remove LSURnew (see main text of manuscript); original data used for manuscript
plantdata <- read.csv("LivingTrees1997-2017.csv", header=T) #Original data from dryad: https://doi.org/10.5061/dryad.ncjsxksvr

#Remove all non-trees
plantdata <- plantdata[plantdata$LifeForm %in% c("tree","midtree","palm","understory palm","Tree","treelet","arbol","Treelet","Palm"),]

#We remove LSUR_new and change BEJ to FEB based on personal communication with R. Chazdon
plantdata <- plantdata[plantdata$PlotName != "LSURnew",]
plantdata[plantdata$PlotName == "BEJ",]$PlotName <- "FEB"
unique(plantdata$PlotName) #check unique names; should not include BEJ or LSUR_new but should include FEB

#Append columns for longitude and latitude of sampling locations and provide coordinates per plot, 
#necessary for the map and the allometric equation for biomass estimation.
plantdata$long <- NA
plantdata$lat <- NA

plot_data <- tibble(Latitude = c(10.45443 , 10.45445, 10.41317, 10.40268, 
                                 10.43142, 10.38987, 10.42397, 10.44070),
                    Longitude = c(-84.06558, -84.06535, -84.02792, 
                                  -84.11167, -84.03382, -84.12835, 
                                  -84.03928, -84.06662),
                    Name = c("JE", "FEB","LSUR","TIR","LEPviejo","CR",
                             "LEP", "SV"),
                    Foresttype = c("Secondary","Secondary",
                                   "Secondary","Secondary",
                                   "Primary","Secondary",
                                   "Secondary","Primary"))

for(i in 1:length(unique(plantdata$PlotName))){
plantdata[plantdata$PlotName == unique(plantdata$PlotName)[i], "long"] <- plot_data[plot_data$Name == unique(plantdata$PlotName)[i],"Longitude"] 
plantdata[plantdata$PlotName == unique(plantdata$PlotName)[i], "lat"] <- plot_data[plot_data$Name == unique(plantdata$PlotName)[i],"Latitude"] 
}

# create Count column
plantdata$Count = rep(1,length(rownames(plantdata))) #Set count 

# =====================================================
# Transforming tree DBH to biomass and metabolic energy
# =====================================================

# Step 1. Wood density: first harmonize taxonomy then retrieve known densities ----------------------
Taxo <- correctTaxo(genus   = plantdata$Genus,
                    species = plantdata$SpeciesName)

WDdata <- getWoodDensity(genus   = Taxo$genusCorrected,
                         species = Taxo$speciesCorrected,
                         stand   = plantdata$PlotName)

#Step 2: create filtered plantdataset to have only unique rows for the species
plantdata_filtered <- plantdata %>% 
                            distinct(PlotName,SpeciesCode, .keep_all = TRUE) %>%
                            select(PlotName,Family, Genus, SpeciesName, WD)
colnames(plantdata_filtered) <- c("PlotName","family", "genus","species","meanWD")

# Step 3A: check for each plot how complete the WD data is to get a sense of how much is estimated
completeness <- plantdata_filtered %>%
  group_by(PlotName) %>%
  summarise(
    total = n(),
    # Count rows where WD is missing OR outside the 0-1 range
    invalid_WD = sum(is.na(meanWD) | meanWD < 0 | meanWD > 1),
    # Count rows where WD is present AND within 0-1
    valid_WD = sum(!is.na(meanWD) & meanWD >= 0 & meanWD <= 1),
    # Valid / Total
    completeness = round((valid_WD / total) * 100, 1)
  ) %>%
  arrange(completeness)

completeness #show tibble of completeness -> with original data completenes equals zero as no WD was provided

#Step 3B: make a full comparison table to show difference between own and retrieved WD
comparison <- plantdata_filtered %>%
  inner_join(WDdata, by = c("family", "genus", "species"), suffix = c("_obs", "_ref"), relationship="many-to-many") #join filtered plantdata and WDdata for comparison
comparison <- na.omit(comparison) #remove all rows for which there are no WD observed -> NOTE: this breaks the code here because we had no WD data at all
comparison$code <- paste(comparison$genus,"_",comparison$species) #create species code
comparison <- comparison[!duplicated(comparison$code), ] #only look at unique species
comparison <- comparison[comparison$meanWD_obs<1, ] #only look at values for WD between 0 and 1
hist(abs(as.numeric(comparison$meanWD_obs)-as.numeric(comparison$meanWD_ref))) #plot differences between observed and estimated WD as histogram -> NOTE: breaks because meanWD is all NA in our case

# Step 3C: choose strategy, either using provided WD and filling or retrieving all WD from scratch
use_measured_WD <- FALSE   # TRUE  = use plantdata$WD where available, fill missing with getWoodDensity
                           # FALSE = use getWoodDensity estimate for every tree

if (use_measured_WD) {
  # Measured WD where available, estimate elsewhere
  plantdata$WD_used <- ifelse(!is.na(plantdata$WD) & plantdata$WD<1,
                              as.numeric(plantdata$WD),
                              WDdata$meanWD)
  # SD = 0 for measured (assumed perfect), sdWD from database for estimated
  plantdata$WD_sd   <- ifelse(!is.na(plantdata$WD) & plantdata$WD<1,
                              0,
                              WDdata$sdWD)
} else {
  # Use estimate for every tree, with its associated uncertainty
  plantdata$WD_used <- WDdata$meanWD
  plantdata$WD_sd   <- WDdata$sdWD
}
  
# Step 4: estimate height for every tree
coord <- cbind(plantdata$long, plantdata$lat)
Hfeld <- retrieveH(D = as.numeric(plantdata$DBH), coord = coord)
plantdata$H <- Hfeld$H

# Step 5: calculate per-tree AGB, specifying new column B for plantdata object
plantdata$B <- computeAGB(D  = as.numeric(plantdata$DBH),
                          WD = as.numeric(plantdata$WD_used),
                          H  = plantdata$H) * 1000 # multiply by 1000 because BIOMASS estimates Mg 

#Step 6: compute E for each individual tree from B using 3/4 power
plantdata$E <- plantdata$B^(3/4) 

#Step 7: save dataset to be used for the error propagation steps in file "10_FigureS4ABC.R"
write.csv(plantdata, file = "plantdata_for_error_propagation.csv", row.names = FALSE)

# =======================================================
# Create dataset to work with for the subsequent analyses
# =======================================================

plantdata$PlotCode <- paste(plantdata$PlotName,
							plantdata$Year, sep="_")

treedata <- plantdata[,c("PlotCode",
                          "Foresttype","SpeciesCode","Count","B","E",
                          "PlotName", "Year","Quadrat")]

# Normalize E and B so smallest is 1 (following original paper)
treedata$E <- treedata$E/min(treedata$E)
treedata$B <- treedata$B/min(treedata$B)

#Save treedata to use for generating Figure 04 using file "06_Figure_04.R"
write.csv(treedata, file = "treedata.csv", row.names = FALSE)

# =============================================================
# BOOTSTRAPPING BIOMASS ESTIMATIONS
#
# Generating bootstrapped data for estimates of ecological 
# metrics and state variables with bootstrap means and SE
#
# Metrics computed per PlotCode + Year, per iteration:
#   S        species richness
#   N        total abundance
#   H        Shannon diversity, -sum(p_i * ln p_i)
#   Evenness Hill evenness, exp(H) / S
#   B, E     summed (or robustly aggregated) plot metrics
# =============================================================

final_aggdata <- bootstrap_treedata(
  data        = treedata,
  n_boot      = 2000,
  sample_frac = 1,
  replace     = TRUE,
  normalize   = FALSE,
  return_replicates = TRUE,
  output_file = "treedata_bootstrap.Rdata"
)

#Save bootstrap summary data file needed for figure S4B
write.csv(final_aggdata$summary, file = "bootstrap_summary.csv", row.names = FALSE

# Tidying up the final dataframe
aggdata <- data.frame(final_aggdata$summary[,c("PlotCode","Foresttype","PlotName","Year","B_mean","E_mean","S_mean","N_mean","Evenness_mean")])
colnames(aggdata) <- c("PlotCode","Foresttype","Plot","Year","B","E","S","N","Evenness")

#Save aggregated data file if using for the first time or data was not provided
write.csv(aggdata, file = "state_variables_CR.csv", row.names = FALSE)

# =============================
# CALCULATING EEoS PREDICTIONS 
# =============================

# Pull the B values for all individual plots of the final measurement 2017 and order them in decreasing way
# Select all numeric columns
#biomass_last <- filter(aggdata,Year=="2017")
#df_sorted_desc <- arrange(biomass_last, desc(B))

# load data from here if provided
aggdata <- read.csv("state_variables_CR.csv", header=T)

# Make new dataset from working file
df <- aggdata

# Add columns for predicted numerical biomass data
df$pBnum_exact <- rep(0, nrow(df))
df$pBnum_approx <- rep(0, nrow(df))

# Iterate through each row and append the biomass information
power <- 4/3 # used by default, but can be changed to test other scaling relationships

## Exact calculation
pb <- txtProgressBar(min = 0, max = nrow(df), style = 3)
for (i in 1:nrow(df)) {
  s <- df[i, c("B", "E", "S", "N")]   # state variables
  l <- mete_lambdas(s)                # get lambdas
  z <- zfun(l, s)$value               # normalization
  b <- s$S * integrate(function(loge) {
         exp((power + 1) * loge) * nRsum(exp(loge), l, s)
       },
       lower = 0, upper = log(s$E))$value
  df$pBnum_exact[i] <- (b / z)

  setTxtProgressBar(pb, i)
}
close(pb)

#b <- log((s$N)^2/(s$S))*s$E^(4/3)/(s$S^(1/3)*log(1/l[1]+l[2])) #new
# exact is summation over N and E of R times metabolic rate to the 4/3 

# Approximation for the biomass equation with power=4/3. 
# Change order = 1 if you want to include correction.
order <- 0

pb <- txtProgressBar(min = 0, max = nrow(df), style = 3)
for(i in 1:nrow(df)){
  s <- df[i,c("B","E","S","N")] #change to state
  l <- mete_lambdas(s)
 
  b <- 4.17*s$E^(4/3)/(s$S^(1/3)*log(1/l[1]+l[2])) #the approximation

    if (order == 1) {
  b <- b * (1-1.16*(l[1]+l[2])^(1/3))}
  df$pBnum_approx[i] <- b
  
  setTxtProgressBar(pb, i)
  }
close(pb)

# ===================================================
# SETTING UP FINAL DATASET FOR GENERATING ALL FIGURES
# ===================================================

# Create clean dataset to use for plotting
trees <- df

# Set levels for plotting (color and symbol)
trees$Year <- factor(trees$Year)

# Calculate distance between prediction and obs. for all plots
trees$dist <- ((trees$pBnum_exact)-(trees$B))/(trees$B)

# write out to file for rest of analyses
write.csv(trees,"trees.csv", row.names = FALSE)
