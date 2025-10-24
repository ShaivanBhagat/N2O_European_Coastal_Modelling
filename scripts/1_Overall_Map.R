# ============================================================
# Script 1 : Overall_Map.R
# Purpose: Overall European Plots
# ============================================================

# --- Libraries ---
library(readxl)
library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(maps)
library(cowplot)
library(stringr)
library(MASS)
library(patchwork)

# --- 1. Loading database and ensuring decimal format ---
df <- read_excel("Data_N2O_Guided Research_Shaivan_v1.xlsx")

df <- df %>%
  mutate(`N2O Emission (μmol N₂O m⁻² d⁻¹)` = gsub(",", ".", `N2O Emission (μmol N₂O m⁻² d⁻¹)`))

# --- 2. Cleaning & renaming  ---
df_clean <- df %>%
  rename(lat = `Latitude (°)`,
         lon = `Longitude (°)`,
         flux_raw = `N2O Emission (μmol N₂O m⁻² d⁻¹)`,
         type = `Point or Region`) %>%
  mutate(
    lat = as.numeric(lat),
    lon = as.numeric(lon),
    # normalize flux
    flux = flux_raw %>%
      as.character() %>%                     # make sure it's text
      str_replace_all("\u00A0", "") %>%      # remove non-breaking spaces
      str_replace_all("−", "-") %>%          # replace unicode minus with ASCII -
      str_replace_all(",", ".") %>%          # swap comma to dot
      str_trim() %>%                         # remove stray spaces
      as.numeric()
  )

# ========== 3) Keep rows inside a broader Europe bounding box ==========
# (keeping rows with NA coords for the df itself;  only drop them for mapping later)
df_europe <- df_clean %>%
  filter(
    is.na(lon) | (lon >= -35 & lon <= 50),
    is.na(lat) | (lat >= 20  & lat <= 80)
  )

# ========== 4) Spatialize only rows with coordinates ==========
df_map <- df_europe %>% filter(!is.na(lat), !is.na(lon))
df_sf  <- st_as_sf(df_map, coords = c("lon", "lat"), crs = 4326)

# ========== 5) Load coastline & compute nearest distance (meters) ==========
# Use a SINGLE unioned coastline geometry to avoid a huge distance matrix.
coastline <- ne_download(scale = 10, type = "coastline",
                         category = "physical", returnclass = "sf")
coastline_u <- st_union(coastline)  # single geometry

# st_distance vs a single geometry returns a vector (fast, memory-light)
dist_vec <- as.numeric(st_distance(df_sf, coastline_u))
df_sf$dist_to_coast_m <- dist_vec

# ========== 6) Keep only coastal points (<45 km). ==========
df_coastal_sf <- df_sf %>% filter(dist_to_coast_m < 45000)

# ========== 7) Bringing coordinates back for ggplot ==========
coords <- st_coordinates(df_coastal_sf)

df_coastal <- cbind(df_coastal_sf, coords) %>%
  rename(lon = X, lat = Y) %>%
  st_drop_geometry()

# ========== 8) Flux bins for binned map & pie ==========
df_coastal <- df_coastal %>%
  mutate(
    flux_cat = cut(
      flux,
      breaks = c(-Inf, -10, -5, -2, -1, -0.5, 0, 0.5, 1, 2, 5, 10, 50, Inf),
      labels = c("< -10", "-10 to -5", "-5 to -2", "-2 to -1",
                 "-1 to -0.5", "-0.5 to 0", "0 to 0.5", "0.5 to 1",
                 "1 to 2", "2 to 5", "5 to 10", "10 to 50", "> 50")
    )
  )

# ========== 9) Basemap ==========
world_map <- map_data("world")

# ========== 10) Map with bins + shape by data type ==========
map_plot <- ggplot() +
  geom_polygon(data = world_map,
               aes(x = long, y = lat, group = group),
               fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_point(
    data = df_coastal,
    aes(x = lon, y = lat, fill = flux_cat, shape = type),
    size = 2, stroke = 0.2, alpha = 0.9
  ) +
  scale_shape_manual(
    values = c(
      "Point"      = 21,  # circle (filled)
      "Region"     = 22,  # square (filled)
      "Aggregated" = 24   # triangle (filled)
    ),
    na.translate = FALSE
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, color = "black"))
  ) +
  scale_fill_manual(
    values = c(
      "< -10"    = "#2171B5",
      "-10 to -5"= "#4292C6",
      "-5 to -2" = "#6BAED6",
      "-2 to -1" = "#9ECAE1",
      "-1 to -0.5" = "#C6DBEF",
      "-0.5 to 0"  = "#DEEBF7",
      "0 to 0.5"   = "#FFF7BC",
      "0.5 to 1"   = "#FEE391",
      "1 to 2"     = "#FEC44F",
      "2 to 5"     = "#FE9929",
      "5 to 10"    = "#EC7014",
      "10 to 50"   = "#CC4C02",
      "> 50"       = "#993404"
    ),
    na.value = "black"
  ) +
  scale_shape_manual(
    values = c(
      "Point"      = 21, # circle
      "Region"     = 22,  # square
      "Aggregrated" = 24   # triangle
    ),
    na.translate = FALSE
  ) +
  coord_fixed(xlim = c(-35, 50), ylim = c(20, 80)) +
  theme_bw() +   # instead of theme_void()
  labs(
    fill = expression(N[2]*O~flux~bgroup("(", µmol~m^-2~d^-1, ")")),
shape = "Data type", title = "Observed N2O Fluxes in European Coastal Systems",
x = "Longitude (°E)",
y = "Latitude (°N)"
) +
  theme(
    panel.grid = element_line(),   # clean background
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

# ========== 11) Pie chart of flux bins ==========
pie_data <- df_coastal %>%
  count(flux_cat) %>%
  mutate(percent = 100 * n / sum(n))

pie_plot <- ggplot(pie_data, aes(x = "", y = percent, fill = flux_cat)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c(
    "< -10"        = "#084594",  # deep blue
    "-10 to -5"    = "#2171B5",
    "-5 to -2"     = "#4292C6",
    "-2 to -1"     = "#6BAED6",
    "-1 to -0.5"   = "#9ECAE1",
    "-0.5 to 0"    = "#C6DBEF",
    "0 to 0.5"     = "#FDD49E",
    "0.5 to 1"     = "#FDBB84",
    "1 to 2"       = "#FC8D59",
    "2 to 5"       = "#E34A33",
    "5 to 10"      = "#B30000",
    "10 to 50"     = "#7F0000",
    "> 50"         = "#4D0000"
  )) +
  theme_void() +
  theme(legend.position = "none")

# ========== 12) Combine map + pie inset ==========
final_plot <- cowplot::ggdraw() +
  draw_plot(map_plot) +
  draw_plot(pie_plot, x = 0.78, y = 0.10, width = 0.18, height = 0.18)

print(map_plot)
# print(final_plot)

# ========== 13) Save outputs ==========
# write.csv(df_coastal, "European_Coastal_N2O_Fluxes.csv", row.names = FALSE)
# ggsave("European_Coastal_N2O_Fluxes_Map.png", map_plot, width = 10, height = 7, dpi = 300)


# quick summaries
print(
  df_coastal %>%
    group_by(type) %>%
    summarize(
      n = n(),
      mean_flux = mean(flux, na.rm = TRUE),
      median_flux = median(flux, na.rm = TRUE),
      .groups = "drop"
    )
)


# ========== Aggregrated Plot ==========
map_plot_aggregrated <- ggplot() +
  geom_polygon(data = world_map,
               aes(x = long, y = lat, group = group),
               fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_point(
    data = df_coastal %>% filter(type == "Aggregrated"),
    aes(x = lon, y = lat, fill = flux_cat),
    size = 2, color = "black", shape = 21, stroke = 0.2, alpha = 0.9
  ) +
  scale_fill_manual(
    values = c(
      "< -10"    = "#2171B5",
      "-10 to -5"= "#4292C6",
      "-5 to -2" = "#6BAED6",
      "-2 to -1" = "#9ECAE1",
      "-1 to -0.5" = "#C6DBEF",
      "-0.5 to 0"  = "#DEEBF7",
      "0 to 0.5"   = "#FFF7BC",
      "0.5 to 1"   = "#FEE391",
      "1 to 2"     = "#FEC44F",
      "2 to 5"     = "#FE9929",
      "5 to 10"    = "#EC7014",
      "10 to 50"   = "#CC4C02",
      "> 50"       = "#993404"
    ),
    na.value = "black"
  ) +
  coord_fixed(xlim = c(-10, 40), ylim = c(35, 70)) +
  theme_bw() +   # <-- keep axes & grid
  labs(
    fill = expression(N[2]*O~flux~bgroup("(", µmol~m^-2~d^-1, ")")),
    title = "Observed N2O Fluxes (For Aggregrated Values)",
    x = "Longitude (°E)",
    y = "Latitude (°N)"
  ) +
  theme(
    panel.grid.major = element_line(color = "grey80", linewidth = 0.4), # main gridlines
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.2), # finer gridlines
    axis.text  = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

print(map_plot_aggregrated)
# ggsave("Aggregreated_Plot.png", map_plot_aggregrated, width = 10, height = 8, dpi = 600)

# ========== Region Plot ==========
map_plot_region <- ggplot() +
  geom_polygon(data = world_map,
               aes(x = long, y = lat, group = group),
               fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_point(
    data = df_coastal %>% filter(type == "Region"),
    aes(x = lon, y = lat, fill = flux_cat),
    size = 2, color = "black", shape = 21, stroke = 0.2, alpha = 0.9
  ) +
  scale_fill_manual(
    values = c(
      "< -10"    = "#2171B5",
      "-10 to -5"= "#4292C6",
      "-5 to -2" = "#6BAED6",
      "-2 to -1" = "#9ECAE1",
      "-1 to -0.5" = "#C6DBEF",
      "-0.5 to 0"  = "#DEEBF7",
      "0 to 0.5"   = "#FFF7BC",
      "0.5 to 1"   = "#FEE391",
      "1 to 2"     = "#FEC44F",
      "2 to 5"     = "#FE9929",
      "5 to 10"    = "#EC7014",
      "10 to 50"   = "#CC4C02",
      "> 50"       = "#993404"
    ),
    na.value = "black"
  ) +
  coord_fixed(xlim =  c(-10, 40), ylim = c(35, 70)) +
  theme_bw() +   # <-- keep axes & grid
  labs(
    fill = expression(N[2]*O~flux~bgroup("(", µmol~m^-2~d^-1, ")")),
    title = "Observed N2O Fluxes (For Region values)",
    x = "Longitude (°E)",
    y = "Latitude (°N)"
  ) +
  theme(
    panel.grid.major = element_line(color = "grey80", linewidth = 0.4), # main gridlines
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.2), # finer gridlines
    axis.text  = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

print(map_plot_region)
# ggsave("Region_plot.png", map_plot_region, width = 10, height = 7, dpi = 300)

# ========== Point Plot ==========
map_plot_point <- ggplot() +
  geom_polygon(data = world_map,
               aes(x = long, y = lat, group = group),
               fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_point(
    data = df_coastal %>% filter(type == "Point"),
    aes(x = lon, y = lat, fill = flux_cat),
    size = 2, color = "black", shape = 21, stroke = 0.2, alpha = 0.9
  ) +
  scale_fill_manual(
    values = c(
      "< -10"    = "#2171B5",
      "-10 to -5"= "#4292C6",
      "-5 to -2" = "#6BAED6",
      "-2 to -1" = "#9ECAE1",
      "-1 to -0.5" = "#C6DBEF",
      "-0.5 to 0"  = "#DEEBF7",
      "0 to 0.5"   = "#FFF7BC",
      "0.5 to 1"   = "#FEE391",
      "1 to 2"     = "#FEC44F",
      "2 to 5"     = "#FE9929",
      "5 to 10"    = "#EC7014",
      "10 to 50"   = "#CC4C02",
      "> 50"       = "#993404"
    ),
    na.value = "black"
  ) +
  coord_fixed(xlim = c(-35, 50), ylim = c(35, 80)) +
  theme_bw() +   # <-- keep axes & grid
  labs(
    fill = expression(N[2]*O~flux~bgroup("(", µmol~m^-2~d^-1, ")")),
    title = "Observed N2O Fluxes (For Individual Points)",
    x = "Longitude (°E)",
    y = "Latitude (°N)"
  ) +
  theme(
    panel.grid.major = element_line(color = "grey80", linewidth = 0.4), # main gridlines
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.2), # finer gridlines
    axis.text  = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

print(map_plot_point)
# ggsave("Point_plot.png", map_plot_point, width = 10, height = 8, dpi = 600)


Plot1 <- map_plot_aggregrated + map_plot_region + 
  plot_layout(ncol = 1)
print(Plot1)
# ggsave("Plot1.png", Plot1, width = 10, height = 8, dpi = 600)


# ========== Density Plots ==========
# Keep only points with coordinates and flux
# plots for all regions, points, aggregated
df_heat <- df_coastal %>%
  filter(!is.na(lon), !is.na(lat), !is.na(flux))

# ---- Kernel density heatmap ----
density_plot <- ggplot(df_heat, aes(x = lon, y = lat)) +
  stat_density2d_filled(aes(fill = after_stat(level)), alpha = 0.8, contour_var = "ndensity") +
  geom_point(data = df_heat, aes(x = lon, y = lat), size = 0.5, color = "black", alpha = 0.5) +
  borders("world", colour = "grey40", fill = "grey90") +
  coord_fixed(xlim = c(-20, 40), ylim = c(35, 70)) +
  theme_bw() +
  labs(
    title = "N2O Flux Hotspots (Density Map)",
    x = "Longitude (°E)", y = "Latitude (°N)",
    fill = "Relative Density"
  )

print(density_plot)
# ggsave("Density_plot_1.png", density_plot, width = 10, height = 7, dpi = 600)

# ----- Updated (Colour Scheme) Density Heatmap (of Observations)-----
# This map shows the relative sampling density of N₂O flux observations across Europe.
# The scale is unitless, ranging from low density (yellow, few observations) to 
# high density (dark purple/black, many observations clustered).
# It highlights hotspots of data collection, not hotspots of flux magnitude.
# This map shows the spatial density of N₂O flux observations across European coastal systems. 
# The scale is unitless and reflects where measurements have been most frequently collected,
# highlighting data-rich regions such as the North Sea and Baltic Sea while revealing areas with limited coverage.

density_plot_2 <- ggplot(df_heat, aes(x = lon, y = lat)) +
  stat_density_2d(
    aes(fill = after_stat(ndensity)),
    geom = "raster", contour = FALSE, n = 200
  ) +
  geom_point(data = df_heat, aes(x = lon, y = lat),
             color = "black", size = 0.5, alpha = 0.3) +
  borders("world", colour = "grey40", fill = "grey90") +
  coord_fixed(xlim = c(-20, 40), ylim = c(35, 70)) +
  theme_bw() +
  labs(
    title = "N2O Flux Observation Density",
    x = "Longitude (°E)", y = "Latitude (°N)",
    fill = "Relative Density"
  ) +
  scale_fill_viridis_c(option = "ndensity")

print(density_plot_2)
# ggsave("Density_plot_2.png", density_plot_2, width = 10, height = 7, dpi = 600)



# Regions where both flux values and point density combine. In the kernel density estimation (KDE), 
# I have weighted the distribution by emission values, so grid cells with higher flux contribute more strongly
# Compute KDE weighted by flux

# Flux-Weighted Hotspot Map
# This map shows flux-weighted kernel density estimates of N₂O emissions, 
# where both the number of observations and the magnitude of flux values are combined.
# The resulting hotspots identify regions that contribute disproportionately to emissions, 
#  revealing environmental importance rather than just data availability.
df_heat <- df_coastal %>%
  filter(!is.na(lon), !is.na(lat), !is.na(flux))

dens <- kde2d(df_heat$lon, df_heat$lat,
              n = 200,
              lims = c(-35, 50, 20, 80))  # bounding box

# Expand to dataframe
dens_df <- expand.grid(lon = dens$x, lat = dens$y)
dens_df$z <- as.vector(dens$z)

# Plot: flux-weighted "hotspot" map
# Scale: “Flux-weighted density” — not a direct physical unit, but proportional to where higher emissions cluster.
flux_density_plot <- ggplot(dens_df, aes(x = lon, y = lat, fill = z)) +
  geom_raster(alpha = 0.9, interpolate = TRUE) +
  geom_contour(aes(z = z), color = "black", size = 0.2, alpha = 0.4) +
  geom_point(data = df_heat, 
             mapping = aes(x = lon, y = lat),   # <- explicitly re-map here
             color = "grey90", size = 0.5, alpha = 0.5, inherit.aes = FALSE) +
  borders("world", colour = "grey40", fill = "grey90") +
  coord_fixed(xlim = c(-20, 40), ylim = c(35, 70)) +
  scale_fill_viridis_c(option = "magma", trans = "sqrt",
                       name = expression("Flux-weighted density")) +
  theme_bw() +
  labs(
    title = "Flux-Weighted Hotspots of N2O Emissions",
    x = "Longitude (°E)", y = "Latitude (°N)"
  ) +
  theme(panel.grid = element_line(color = "grey90"))

print(flux_density_plot)

# ggsave("density.png", flux_density_plot, width = 10, height = 8, dpi = 600)


