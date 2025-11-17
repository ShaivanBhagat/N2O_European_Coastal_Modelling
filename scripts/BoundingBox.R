## Packages --------------------------------------------------------------
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(RColorBrewer)

## Basemap (Europe) ------------------------------------------------------
world  <- ne_countries(scale = "medium", returnclass = "sf")

# Crop to the area you had in the figure
world_eur <- st_crop(
  world,
  xmin = -40, xmax = 45,
  ymin =  30, ymax = 82
)

## Study-region bounding boxes -------------------------------------------
regions <- tribble(
  ~region,               ~xmin, ~xmax, ~ymin, ~ymax,
  "Northeast Atlantic",  -34,     -5,    48,    75,
  "Scandinavian Seas",   -5,    45,    66,    80,
  "North Sea",           -5,     12,    48,    66,
  "Baltic Sea",          12,     31,    53,    66,
  "Western Mediterranean",-15,   12,    34,    48,
  "Eastern Mediterranean",12,    27,    34,    45,
  "Eastern Mediterranean",15,    47,    30,    40,
  "Black Sea",           27,     42,    40,    47
) |>
  mutate(
    lon = (xmin + xmax) / 2,
    lat = (ymin + ymax) / 2
  )

## Plot ------------------------------------------------------------------
ggplot() +
  # grey basemap
  geom_sf(data = world_eur,
          fill = "grey90", colour = "grey70", linewidth = 0.3) +
  
  # coloured bounding boxes
  geom_rect(
    data  = regions,
    aes(xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax,
        fill = region),
    colour = "black",
    alpha  = 0.35,
    linewidth = 0.4
  ) +
  
  # region labels
  geom_text(
    data = regions,
    aes(x = lon, y = lat, label = region),
    size = 3.5
  ) +
  
  coord_sf(
    xlim = c(-40, 45),
    ylim = c(30, 82),
    expand = FALSE
  ) +
  
  # pastel colours similar to your figure
  scale_fill_brewer(palette = "Pastel1", guide = "none") +
  
  labs(
    x = "Longitude (°E)",
    y = "Latitude (°N)",
    title = "Study Regions (Approximate Bounding Boxes)"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )





