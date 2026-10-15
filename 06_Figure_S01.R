# =====================================================================
# Code for
# A coherent signature of disturbance in successional forests
# Code written and maintained by Edwin Pos and Erica Newman
# E.T.Pos@uu.nl, ericanewman.sci@gmail.com
# Updated October 6, 2026
# =====================================================================

# ===============
# FIGURE S1 (MAP)
# ===============

## Code derived from the following sources:
#https://www.r-bloggers.com/2022/10/map-any-region-in-the-world-with-r-part-i-the-basic-map/
# https://gis.stackexchange.com/questions/222799/creating-inset-map-using-r

#remotes::install_github("ropensci/rnaturalearthhires")

# install.packages(
#  "rnaturalearthhires",
#  repos = "https://ropensci.r-universe.dev",
#  type = "source"
#  )

#install.packages("sf")
#install.packages("cowplot")
#install.packages("rnaturalearth")
#install.packages("ggrepel")
#install.packages("jpeg") 

library(sf)
library(cowplot)
library(rnaturalearth)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(jpeg)
library(stringr)

# Get data for map
world_df <- map_data("world")
cr_df <- map_data('world')[map_data('world')$region == "Costa Rica",]

# Create plot data
plot_data <- tibble(Latitude = c(10.45443 , 10.45445, 10.41317, 10.40268, 
                                 10.43142, 10.38987, 10.42397, 10.44070),
                    Longitude = c(-84.06558, -84.06535, -84.02792, 
                                  -84.11167, -84.03382, -84.12835, 
                                  -84.03928, -84.06662),
                    Name = c("JE", "BEJ","LSUR","TIR","LEPVIEJO","CR",
                             "LEP", "SV"),
                    Foresttype = c("Secondary","Secondary",
                                   "Secondary","Secondary",
                                   "Secondary","Secondary",
                                   "Primary","Primary"))

# Transform data points into geographic objects
Sites_geo <- plot_data %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

# Get countries outline
pays <- ne_countries(
  country = c("Costa Rica"),
  scale = "large", returnclass = "sf")

# Base plot
main <- ggplot(data = pays) +
  geom_sf() +
  xlim(c(86,82.5))+
  ylim(c(7.9,11.1))+
  geom_sf(data = Sites_geo, aes(fill = NA), 
          shape = c(3,3,3,3,3,3,3,3),size = 2,
          alpha = 0.7,
          show.legend=FALSE) +
  scale_fill_viridis_d() +
  labs(x = "", y = "") +
  theme_bw()
main

# Add scale and North arrow
main_final = main +
  ggspatial::annotation_scale(
    location = "br",
    bar_cols = c("grey60", "white")
  ) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    pad_x = unit(0.4, "in"), pad_y = unit(0.4, "in"),
    style = ggspatial::north_arrow_minimal(
      #fill = c("grey40", "white"),
      line_col = "grey20"
    )
  )

main_final

##################
inset_1 <- ggplot(data = pays) +
  geom_sf() +
  geom_sf(data = Sites_geo, aes(fill = NA), shape = c(21,21,21,21,21,21,22,22), size=2.5, show.legend=F) +
  geom_text_repel(
    data = plot_data, aes(x = Longitude, y = Latitude, label = Name),
    force = 4, size=3) +
  coord_sf(xlim = c(-84.14, -84.02), ylim = c(10.365, 10.47), expand = FALSE) +
  ggspatial::annotation_scale(
    location = "bl",
    bar_cols = c("grey60", "white")
  ) +
  theme_void() +
  theme(panel.background = element_rect(colour = "black", 
                                        size = 1))

#inset_1

##################
inset_2 <- 
  ggplot(data = world_df) + 
  geom_polygon( aes(x = long, 
                   y = lat, 
                   group = group)
              ) +
  # fill the country in with green
  geom_polygon(data = cr_df, 
               aes(x = long, 
                   y = lat, 
                   group = group),
               fill = "#F67A31"# "#B2C054"
               ) +
  coord_sf(xlim = c(-50, -100), ylim = c(-10, 30), 
           expand = FALSE) +
  theme_void() +
  theme(panel.background = element_rect(colour = "#595959",
                                        size = 1))
#inset_2

#print to png
#png("costarica_overview.png")
#inset_2
#dev.off()

##################

# Place the insets
map_with_inset <-
  ggdraw() +
  draw_plot(main_final) +
  draw_plot(inset_1, x = 0.12, y = .375, width = .25, height = .25)+
  draw_plot(inset_2, x = 0.12, y = 0.10, width = .3, height = .3)

ggsave("Figure_S01.png", 
       plot = map_with_inset, 
       dpi = 300)    # Resolution
