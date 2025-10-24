# ============================================================
# Script 3: Source_Sink.R
# Purpose: Source vs Sink Script
# ============================================================

# --- Libraries ---
library(readxl)
library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# --- 1) Load your database ---
df <- read_excel("Data_N2O_Guided Research_Shaivan_v1.xlsx")

# Clean and rename
df_clean <- df %>%
  rename(lat = `Latitude (°)`,
         lon = `Longitude (°)`,
         flux = `N2O Emission (μmol N₂O m⁻² d⁻¹)`,
         type = `Point or Region`) %>%
  mutate(
    lat  = as.numeric(lat),
    lon  = as.numeric(lon),
    flux = as.numeric(gsub(",", ".", flux))
  )

# --- 2) Classify as Source / Sink / Neutral ---
df_clean <- df_clean %>%
  mutate(SourceSink = case_when(
    flux > 0  ~ "Source",
    flux < 0  ~ "Sink",
    flux == 0 ~ "Neutral",
    TRUE ~ NA_character_
  ))

# --- 3) Keep rows inside broader Europe bounding box ---
df_europe <- df_clean %>%
  filter(
    is.na(lon) | (lon >= -35 & lon <= 50),
    is.na(lat) | (lat >= 20  & lat <= 80)
  )

# --- 4) Spatialize only rows with coordinates ---
df_map <- df_europe %>% filter(!is.na(lat), !is.na(lon))
df_sf  <- st_as_sf(df_map, coords = c("lon", "lat"), crs = 4326)

# --- 5) Load coastline & compute nearest distance (meters) ---
coastline <- ne_download(scale = 10, type = "coastline",
                         category = "physical", returnclass = "sf")
coastline_u <- st_union(coastline)

dist_vec <- as.numeric(st_distance(df_sf, coastline_u))
df_sf$dist_to_coast_m <- dist_vec

# --- 6) Keep only coastal points (< 45 km) ---
df_coastal_sf <- df_sf %>% filter(dist_to_coast_m < 45000)

# --- 7) Bring coords back for ggplot ---
coords <- st_coordinates(df_coastal_sf)

df_coastal <- cbind(df_coastal_sf, coords) %>%
  rename(lon = X, lat = Y) %>%
  st_drop_geometry()

# --- 8) World basemap ---
world_map <- map_data("world")

# --- 9) Geospatial plot of Sources & Sinks ---
map_sourcesink <- ggplot() +
  geom_polygon(data = world_map,
               aes(x = long, y = lat, group = group),
               fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_point(
    data = df_coastal %>% filter(!is.na(SourceSink)),
    aes(x = lon, y = lat,
        color = SourceSink, shape = SourceSink,
        size = abs(flux)),      # scale point size by absolute flux
    alpha = 0.9
  ) +
  scale_color_manual(values = c("Source" = "#E31A1C",   # red
                                "Sink"   = "#2171B5",   # blue
                                "Neutral"= "grey40")) +
  scale_shape_manual(values = c("Source" = 16, "Sink" = 17, "Neutral" = 15)) +
  scale_size_continuous(range = c(1, 6), guide = "none") + # control point size
  coord_fixed(xlim = c(-20, 40), ylim = c(35, 70)) +
  theme_bw() +
  labs(
    title = "Geospatial Distribution of Coastal N2O Sources and Sinks",
    x = "Longitude (°E)", y = "Latitude (°N)",
    color = "Flux Category", shape = "Flux Category"
  ) +
  theme(
    panel.grid = element_line(color = "grey90"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

print(map_sourcesink)
# ggsave("4.png", map_sourcesink, width = 12, height = 8, dpi = 600)

# --- 2. Summary table by region ---
df_region <- df_coastal %>%
  mutate(Region = case_when(
    lon >= -5  & lon <= 15  & lat >= 48 & lat <= 65 ~ "NorthSea",
    lon >= 15  & lon <= 31  & lat >= 53 & lat <= 65 ~ "Baltic",
    lon >= -20 & lon <= 12  & lat >= 33 & lat <= 48 ~ "MedWest",
    lon >= 12  & lon <= 37  & lat >= 33 & lat <= 48 ~ "MedEast",
    lon >= -8   & lon <= 45  & lat >= 65 & lat <= 80 ~ "Scand",
    lon >= 27  & lon <= 42  & lat >= 40 & lat <= 48 ~ "BlackSea",
    lon >= -35 & lon <= -5  & lat >= 45 & lat <= 75 ~ "Atlantic",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Region))



source_sink_summary <- df_region %>%
  group_by(Region, SourceSink) %>%
  summarise(n = n(),
            mean_flux = mean(flux, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(Region) %>%
  mutate(percent = round(100 * n / sum(n), 1))

print(source_sink_summary)
