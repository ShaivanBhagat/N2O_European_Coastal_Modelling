# ============================================================
# Script 2 : Zonal_Plot.R
# Purpose: Plot N2O fluxes by coastal zones (7 regions)
# ============================================================

# --- Libraries ---
library(readxl)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(cowplot)
library(stringr)
library(purrr)
library(tidyr)
library(MASS)
library(viridis)
library(patchwork)

# --- 1. Loading database and ensuring decimal format ---
df <- read_excel("Data_N2O_Guided Research_Shaivan_v1.xlsx")

df <- df %>%
  mutate(`N2O Emission (μmol N₂O m⁻² d⁻¹)` = gsub(",", ".", `N2O Emission (μmol N₂O m⁻² d⁻¹)`))

# --- 2. Cleaning & renaming  ---
df <- df %>%
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
df_europe <- df %>%
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

# ========== 6) Keep only coastal points  ==========
df_coastal_sf <- df_sf %>% filter(dist_to_coast_m < 45000)

# ========== 7) Bringing coordinates back for ggplot ==========
coords <- st_coordinates(df_coastal_sf)

df_coastal <- cbind(df_coastal_sf, coords) %>%
  rename(lon = X, lat = Y) %>%
  st_drop_geometry()

# --- 8) Flux bins for binned map & pie ---
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

# --- 9. Filter for plotting if required  ---
df_points <- df_coastal %>% filter(type == "Point")

# --- 10. World basemap ---
world_map <- map_data("world")

# --- 11. Define bounding boxes for zones ---
regions <- list(
  NorthSea   = list(xlim = c(-5, 15),   ylim = c(48, 65)),
  Baltic     = list(xlim = c(15, 31),   ylim = c(53, 65)),
  MedWest    = list(xlim = c(-20, 12),  ylim = c(33, 48)),
  MedEast    = list(xlim = c(12, 37),   ylim = c(33, 48)),
  Scand      = list(xlim = c(-8, 45),    ylim = c(65, 80)),
  BlackSea   = list(xlim = c(27, 42),   ylim = c(40, 48)),
  Atlantic   = list(xlim = c(-35, -5),  ylim = c(45, 75))
)

# --- 12. Function to plot one region ---
plot_region <- function(region_name, df_plot) {
  bbox <- regions[[region_name]]
  
  ggplot() +
    geom_polygon(data = world_map,
                 aes(x = long, y = lat, group = group),
                 fill = "grey95", color = "grey70", linewidth = 0.2) +
    geom_point(
      data = df_plot,  # %>% filter(type == "Point"),
      aes(x = lon, y = lat, fill = flux_cat),
      size = 2, color = "black", shape = 21, stroke = 0.2, alpha = 0.9
    ) +
    scale_fill_manual(
      values = c(
        "< -10"     = "#2171B5",
        "-10 to -5" = "#4292C6",
        "-5 to -2"  = "#6BAED6",
        "-2 to -1"  = "#9ECAE1",
        "-1 to -0.5"= "#C6DBEF",
        "-0.5 to 0" = "#DEEBF7",
        "0 to 0.5"  = "#FFF7BC",
        "0.5 to 1"  = "#FEE391",
        "1 to 2"    = "#FEC44F",
        "2 to 5"    = "#FE9929",
        "5 to 10"   = "#EC7014",
        "10 to 50"  = "#CC4C02",
        "> 50"      = "#993404"
      ),
      na.value = "black"
    ) +
    coord_fixed(xlim = bbox$xlim, ylim = bbox$ylim) +
    theme_bw() +
    labs(
      fill = expression(N[2]*O~flux~' ('*µmol~m^{-2}~d^{-1}*')'),
      title = paste("Observed N2O Fluxes –", region_name, " "),
      x = "Longitude (°E)",
      y = "Latitude (°N)"
    ) +
    theme(
      panel.grid.major = element_line(color = "grey80", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
      axis.text  = element_text(size = 10),
      axis.title = element_text(size = 12)
    )
}

# --- 13. Generate all regional plots ---
map_northsea <- plot_region("a) NorthSea", df_points)
map_baltic   <- plot_region("b)Baltic", df_points)
map_medwest  <- plot_region("f) MedWest", df_points)
map_medeast  <- plot_region("g) MedEast", df_points)
map_scand    <- plot_region("e) Scand", df_points)
map_blacksea <- plot_region("d) BlackSea", df_points)
map_atlantic <- plot_region("c) Atlantic", df_points)

# --- 14. Arrange them in a grid ---
zonal_plot <- cowplot::plot_grid(
  map_northsea, map_baltic, map_medwest,
   map_medeast, map_scand, map_blacksea, map_atlantic,
  ncol = 2
)

# --- 14a. Extract legend from one plot -------
legend <- cowplot::get_legend(
  map_northsea + theme(legend.position = "bottom")
)

# --- 14b. Remove legends from all plots ---
map_northsea <- map_northsea + theme(legend.position = "none")
map_baltic   <- map_baltic   + theme(legend.position = "none")
map_medwest  <- map_medwest  + theme(legend.position = "none")
map_medeast  <- map_medeast  + theme(legend.position = "none")
map_scand    <- map_scand    + theme(legend.position = "none")
map_blacksea <- map_blacksea + theme(legend.position = "none")
map_atlantic <- map_atlantic + theme(legend.position = "none")

# --- 14c. Arrange plots in grid (without legends) ---
zonal_grid <- cowplot::plot_grid(
  map_northsea, map_baltic, map_medwest,
  map_medeast, map_scand, map_blacksea, map_atlantic,
  ncol = 2
)

# --- 14d. Combine grid + single legend ---
zonal_plot <- cowplot::plot_grid(
  zonal_grid, legend
  #rel_widths = c(1, 0.2)  
)

print(zonal_plot)

print(map_northsea)
print(map_baltic)
print(map_blacksea)
print(map_medeast)
print(map_medwest)
print(map_scand)
print(map_atlantic)

map_northsea + map_baltic + map_atlantic + map_scand + legend +
  plot_layout(ncol = 3)

zonal_plot <- (
  map_northsea + map_baltic + map_atlantic + map_scand
) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")  

print(zonal_plot)


zonal_plot_2 <- (map_blacksea + map_medeast + map_medwest) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")
print(zonal_plot_2)

ggsave("3.png", zonal_plot_2, width = 12, height = 8, dpi = 600)

# --- 15. Save output ---
ggsave("European_Coastal_N2O_Fluxes_Zones.png", zonal_plot, width = 12, height = 8, dpi = 600)



## ---- Regional Summaries ------
df_coastal <- df %>% filter(!is.na(lon), !is.na(lat))

# --- Function to filter points in a region ---
filter_region <- function(df, region_bbox) {
  df_coastal %>%
    filter(
      lon >= region_bbox$xlim[1], lon <= region_bbox$xlim[2],
      lat >= region_bbox$ylim[1], lat <= region_bbox$ylim[2]
    )
}

# --- Function to summarize region ---
summarize_region <- function(df, region_name, region_bbox) {
  df_reg <- filter_region(df_coastal, region_bbox)
  
  if (nrow(df_reg) == 0) {
    return(data.frame(
      Region = region_name, n = 0,
      Mean = NA, Median = NA, SD = NA,
      Min = NA, Max = NA
    ))
  }
  
  data.frame(
    Region = region_name,
    n = nrow(df_reg),
    Mean = mean(df_reg$flux, na.rm = TRUE),
    Median = median(df_reg$flux, na.rm = TRUE),
    SD = sd(df_reg$flux, na.rm = TRUE),
    Min = min(df_reg$flux, na.rm = TRUE),
    Max = max(df_reg$flux, na.rm = TRUE)
  )
}

assign_region <- function(df, regions) {
  df$Region <- NA
  for (r in names(regions)) {
    bbox <- regions[[r]]
    in_region <- with(df, lon >= bbox$xlim[1] & lon <= bbox$xlim[2] &
                        lat >= bbox$ylim[1] & lat <= bbox$ylim[2])
    df$Region[is.na(df$Region) & in_region] <- r
  }
  df
}

df_coastal_unique <- assign_region(df_coastal, regions)
summary_table <- df_coastal_unique %>%
  group_by(Region) %>%
  summarise(
    n = n(),
    Mean = mean(flux, na.rm = TRUE),
    Median = median(flux, na.rm = TRUE),
    SD = sd(flux, na.rm = TRUE),
    Min = min(flux, na.rm = TRUE),
    Max = max(flux, na.rm = TRUE)
  )
# --- Apply across all regions ---
summary_list <- lapply(names(regions), function(r) {
  summarize_region(df_coastal, r, regions[[r]])
})

summary_table <- bind_rows(summary_list)
print(summary_table)

nrow(df_europe)

nrow(distinct(df_europe, lon, lat, flux_raw))
sum(is.na(df_europe$flux))
sum(summary_table$n)
nrow(df_coastal)


# --- Save summary ---
# write.csv(summary_table, "Regional_Flux_Summary.csv", row.names = FALSE)

# --- Barplot: mean ± SD per region ---
barplot_flux <- ggplot(summary_table, aes(x = Region, y = Mean, fill = Region)) +
  geom_col(color = "black", alpha = 0.7) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.5) +
  theme_bw() +
  labs(
    y = expression("Mean N"[2]*O~flux~' ('*µmol~m^-2~d^-1*')'),
    x = "Region",
    title = "Mean ± SD of N2O Flux by Region"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

print(barplot_flux)
# ggsave("European_Coastal_N2O_Fluxes_Zones.png", barplot_flux, width = 12, height = 8, dpi = 600)

#### ----- #####
# --- 1) Parse the Year column into usable numbers ---
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


df_region_parsed <- df_region %>%
  mutate(
    Year_chr = as.character(Year),
    # normalize ranges: replace en-dash/em-dash/“to” with hyphen
    Year_chr = str_replace_all(Year_chr, "[\u2012\u2013\u2014\u2015–—]", "-"),
    Year_chr = str_replace_all(Year_chr, "\\bto\\b", "-"),
    # extract all 4-digit years (1900–2099)
    years_found = str_extract_all(Year_chr, "(?<!\\d)(?:19|20)\\d{2}(?!\\d)")
  ) %>%
  mutate(
    Year_start = map_dbl(years_found, ~ if (length(.x) >= 1) as.numeric(.x[1]) else NA_real_),
    Year_end   = map_dbl(years_found, ~ if (length(.x) >= 2) as.numeric(dplyr::last(.x)) else NA_real_),
    # if only one year found, make end = start
    Year_end   = if_else(is.na(Year_end) & !is.na(Year_start), Year_start, Year_end),
    # midpoint year (rounded)
    Year_mid   = if_else(!is.na(Year_start) & !is.na(Year_end),
                         round((Year_start + Year_end) / 2), NA_real_),
    # try direct numeric conversion as fallback for plain "2018" etc.
    Year_num_direct = suppressWarnings(as.numeric(Year)),
    # final usable year
    Year_use = coalesce(Year_mid, Year_num_direct))

# Peek at any rows I still couldn’t parse:
unparsed_examples <- df_region_parsed %>%
  filter(is.na(Year_use)) %>%
  distinct(Year) %>%
  head(20)

print(unparsed_examples)

# --- 2) Build time summaries using Year_use ---
df_time <- df_region_parsed %>%
  filter(!is.na(Year_use)) %>%
  group_by(Year = Year_use, Region) %>%
  summarise(
    mean_flux = mean(flux, na.rm = TRUE),
    sd_flux   = sd(flux,   na.rm = TRUE),
    n         = n(),
    .groups = "drop"
  )

overall_time <- df_region_parsed %>%
  filter(!is.na(Year_use)) %>%
  group_by(Year = Year_use) %>%
  summarise(
    mean_flux = mean(flux, na.rm = TRUE),
    sd_flux   = sd(flux,   na.rm = TRUE),
    n         = n(),
    .groups = "drop"
  )

# --- 3) Plots ---
p_overall <- ggplot(overall_time, aes(x = Year, y = mean_flux)) +
  geom_line(color = "darkred", linewidth = 1) +
  geom_ribbon(aes(ymin = mean_flux - sd_flux, ymax = mean_flux + sd_flux),
              fill = "pink", alpha = 0.3) +
  theme_bw() +
  labs(
    title = "Overall Temporal Trend of N2O Fluxes",
    x = "Year",
    y = expression(N[2]*O~flux~' ('*µmol~m^{-2}~d^{-1}*')')
  )

p_region <- ggplot(df_time, aes(x = Year, y = mean_flux, color = Region)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  theme_bw() +
  labs(
    title = "Regional Temporal Trends of N2O Fluxes",
    x = "Year",
    y = expression(N[2]*O~flux~' ('*µmol~m^{-2}~d^{-1}*')')
  )

# --- 4) Show / save ---
print(p_overall)
print(p_region)


#### ------ By Season/Month --------
# ---- 1) Helpers to parse Month/Season ----
month_map <- setNames(1:12, tolower(month.abb))
month_regex <- paste0("(", paste(names(month_map), collapse="|"), ")")  # jan|feb|...

extract_month_nums <- function(s) {
  # returns vector of month numbers found in a string (Jan..Dec)
  m <- str_extract_all(s, regex(month_regex, ignore_case = TRUE))[[1]]
  if (length(m) == 0) integer(0) else unique(month_map[tolower(m)])
}

mid_of_months <- function(v) {
  # choose a representative month if multiple are present
  if (length(v) == 0) return(NA_integer_)
  v[ceiling(length(v)/2)]
}

month_to_season <- function(m) {
  if (is.na(m)) return(NA_character_)
  c("Winter","Winter","Spring","Spring","Spring",
    "Summer","Summer","Summer",
    "Autumn","Autumn","Autumn","Winter")[m]
}

# ---- 2) Clean Month/Season into month + season ----
df_ms <- df_region %>%
  mutate(
    ms_raw   = `Month/Season`,
    ms_clean = ms_raw %>%
      str_replace_all("–", "-") %>%        # en dash -> hyphen
      str_remove_all("\\b(19|20)\\d{2}\\b") %>% # drop years like 2006
      str_remove_all("\\b\\d{1,2}\\b") %>% # drop stray day numbers
      str_squish() %>%
      tolower(),
    # explicit season keywords
    season_word = case_when(
      str_detect(ms_clean, "spring") ~ "Spring",
      str_detect(ms_clean, "summer") ~ "Summer",
      str_detect(ms_clean, "autumn|fall") ~ "Autumn",
      str_detect(ms_clean, "winter") ~ "Winter",
      TRUE ~ NA_character_
    ),
    # extract any month names present; pick a representative month
    month_num = sapply(ms_clean, function(x) mid_of_months(extract_month_nums(x))),
    # derive season from month if no explicit season word
    season = coalesce(season_word, vapply(month_num, month_to_season, character(1))),
    month  = factor(month.abb[month_num], levels = month.abb)
  )

# ---- 3) Monthly summary (use median for robustness) ----
df_month <- df_ms %>%
  filter(!is.na(month)) %>%
  group_by(Region, month) %>%
  summarise(median_flux = median(flux, na.rm = TRUE),
            n = n(), .groups = "drop")

p_month <- ggplot(df_month, aes(x = month, y = median_flux, fill = Region)) +
  geom_col(position = position_dodge(width = 0.8), color = "black", alpha = 0.8) +
  # Optional: show sample size above bars
  geom_text(aes(label = n), position = position_dodge(width = 0.8),
            vjust = -0.3, size = 3, color = "grey20") +
  theme_bw() +
  labs(
    title = "Monthly Variation of N2O Fluxes (median ± counts)",
    x = "Month", y = expression(N[2]*O~flux~' ('*µmol~m^{-2}~d^{-1}*')')
  ) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# ---- 4) Seasonal summary (robust median) ----
df_season <- df_ms %>%
  filter(!is.na(season)) %>%
  mutate(season = factor(season, levels = c("Winter","Spring","Summer","Autumn"))) %>%
  group_by(Region, season) %>%
  summarise(median_flux = median(flux, na.rm = TRUE),
            n = n(), .groups = "drop")

p_season_clean <- ggplot(df_season, aes(x = season, y = median_flux, fill = Region)) +
  geom_col(position = position_dodge(width = 0.75), color = "black", alpha = 0.85) +
  geom_text(aes(label = n), position = position_dodge(width = 0.75),
            vjust = -0.3, size = 3, color = "grey20") +
  theme_bw() +
  labs(
    title = "Seasonal Variation of N2O Fluxes (median ± counts)",
    x = "Season", y = expression(N[2]*O~flux~' ('*µmol~m^{-2}~d^{-1}*')')
  )

# Show
print(p_month)
print(p_season_clean)


# ggsave("month.png", p_month, width = 10, height = 8, dpi = 600)
# ggsave("month_Season.png", p_season_clean, width = 10, height = 8, dpi = 600)


