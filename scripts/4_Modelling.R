# ============================================================
# Script 4: Modelling.R
# Purpose: Modelling N2O Emissions(/Flux)
# ============================================================

# Core libraries
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(sf)
library(raster)
library(ncdf4)
library(randomForest)
library(rnaturalearth)
library(rnaturalearthdata)
library(maps)
library(cowplot)
library(MASS)
library(scales)
library(viridis)
library(caret)
library(xgboost)
library(Matrix)
library(SHAPforxgboost)
library(mgcv)
library(blockCV)
library(factoextra)
library(terra)

# ---- STEP 1: Load and Clean N2O Observational Data ----
# Load Excel file
df <- read_excel("Data_N2O_Guided Research_Shaivan_v1.xlsx")

# Clean decimal format
df <- df %>%
  mutate(`N2O Emission (μmol N₂O m⁻² d⁻¹)` = gsub(",", ".", `N2O Emission (μmol N₂O m⁻² d⁻¹)`))

# Rename and clean up columns
df_clean <- df %>%
  rename(lat = `Latitude (°)`,
         lon = `Longitude (°)`,
         flux_raw = `N2O Emission (μmol N₂O m⁻² d⁻¹)`,
         type = `Point or Region`) %>%
  mutate(
    lat = as.numeric(lat),
    lon = as.numeric(lon),
    flux = flux_raw %>%
      as.character() %>%
      str_replace_all("\u00A0", "") %>%
      str_replace_all("−", "-") %>%
      str_replace_all(",", ".") %>%
      str_trim() %>%
      as.numeric()
  )

# Keep only valid point observations
df_points <- df_clean %>%
  filter(type == "Point", !is.na(lat), !is.na(lon), !is.na(flux))

# Summary
summary(df_points$flux)


# ---- STEP 2: Extract Environmental Predictors from NetCDF ----
# Setup spatial points for extraction
# Convert sampling coordinates to terra SpatVector
coords_vect <- vect(df_points, geom = c("lon", "lat"), crs = "EPSG:4326")

# ---- Helper function: extract surface variable from NetCDF ----
extract_var <- function(nc_file, varname, label, df) {
  r <- rast(nc_file, subds = varname)
  # Assign CRS explicitly to avoid transformation warnings
  crs(r) <- "EPSG:4326"
  surface <- r[[grep("depth=0", names(r))]]
  vals <- terra::extract(surface, coords_vect)[, 2]
  df[[label]] <- vals
  return(df)
}

# ---- Extract WOA environmental variables ----
df_points <- extract_var("woa2023_Nitrate.nc",    "n_an", "nitrate",     df_points)
df_points <- extract_var("woa2023_Phosphate.nc",  "p_an", "phosphate",   df_points)
df_points <- extract_var("woa2023_DisOxygen.nc",  "o_an", "oxygen",      df_points)
df_points <- extract_var("woa2023_Salinity.nc",   "s_an", "salinity",    df_points)
df_points <- extract_var("woa2023_Temperature.nc","t_an","temperature",  df_points)
df_points <- extract_var("woa2018_Density.nc",    "I_an", "density",     df_points)

# ---- Bathymetry (depth) ----
bathy_raster <- rast("GEBCO_2025_sub_ice.nc")
df_points$depth <- terra::extract(bathy_raster, coords_vect)[, 2]

# ---- Chlorophyll (Copernicus, surface only) ----
chl_file1 <- "cmems_mod_glo_bgc_my_0.25deg_P1M-m_1760291720536.nc"      # 1993–2022
chl_file2 <- "cmems_mod_glo_bgc_myint_0.25deg_P1M-m_1760291822771.nc"   # 2023–2025
r1 <- rast(chl_file1)
r2 <- rast(chl_file2)

surface_names_r1 <- grep("chl_depth=0.50576", names(r1), value = TRUE)
surface_names_r2 <- grep("chl_depth=0.50576", names(r2), value = TRUE)

chl_stack1 <- r1[[surface_names_r1]]
chl_stack2 <- r2[[surface_names_r2]]
chl_full_stack <- c(chl_stack1, chl_stack2)
chlorophyll_r <- mean(chl_full_stack, na.rm = TRUE)

df_points$chlorophyll <- terra::extract(chlorophyll_r, coords_vect)[, 2]


# STEP 3: Prepare dataset for modeling
predictor_vars <- c("temperature", "salinity", "oxygen", "nitrate",
                    "phosphate", "chlorophyll", "density", "depth")
target_var <- "flux"

# Keep complete rows only
df_model <- df_points %>%
  dplyr::select(all_of(c(target_var, predictor_vars))) %>%
  filter(complete.cases(.))

# Summary check
print(nrow(df_model))
summary(df_model)

# STEP 4: Split into Training and Test Sets
set.seed(42)
train_indices <- sample(seq_len(nrow(df_model)), size = 0.8 * nrow(df_model))
train_data <- df_model[train_indices, ]
test_data  <- df_model[-train_indices, ]

# STEP 5: Train Random Forest Model
# ================================================================
rf_model <- randomForest(
  flux ~ .,
  data = train_data,
  ntree = 500,
  importance = TRUE
)

print(rf_model)
varImpPlot(rf_model)

# Evaluate performance on test set
pred_test <- predict(rf_model, newdata = test_data)
cat("Test RMSE:", sqrt(mean((pred_test - test_data$flux)^2)), "\n")
cat("Test R²:", cor(pred_test, test_data$flux)^2, "\n")

# Plot observed vs predicted
ggplot(data.frame(obs = test_data$flux, pred = pred_test),
       aes(x = obs, y = pred)) +
  geom_point(alpha = 0.6, color = "darkgreen") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(x = "Observed Flux", y = "Predicted Flux",
       title = "Random Forest N2O Flux Prediction") +
  theme_minimal()

# ---- STEP 6: Evaluate Model Performance ----
# Predict
test_preds <- predict(rf_model, newdata = test_data)

# R² and RMSE
r_squared <- cor(test_preds, test_data$flux)^2
rmse <- sqrt(mean((test_preds - test_data$flux)^2))

cat("Test R²: ", round(r_squared, 3), "\n")
cat("Test RMSE: ", round(rmse, 3), "\n")

# ---- STEP 7: Visualization and Diagnostics ----
# Variable Importance
varImpPlot(rf_model, main = "Variable Importance (Random Forest)")

# Predicted vs Observed
plot(test_data$flux, test_preds,
     xlab = "Observed Flux", ylab = "Predicted Flux",
     main = "Predicted vs Observed N₂O Flux",
     pch = 16, col = "blue")
abline(0, 1, col = "red", lwd = 2)

# Comparing with Linear Model
lm_model <- lm(flux ~ ., data = train_data)
summary(lm_model)

# ------ # ------
# Fit GAM with smooth terms
gam_model <- gam(flux ~ s(temperature) + s(salinity) + s(oxygen) + 
                   s(nitrate) + s(phosphate) + s(chlorophyll) + 
                   s(density) + s(depth), data = train_data)

# Predict and evaluate
gam_preds <- predict(gam_model, newdata = test_data)
gam_r2 <- cor(gam_preds, test_data$flux)^2
gam_rmse <- sqrt(mean((gam_preds - test_data$flux)^2))

cat("GAM R²:", round(gam_r2, 3), "\n")
cat("GAM RMSE:", round(gam_rmse, 3), "\n")

# ---- STEP 8: Add Coordinates as Predictors ----
# Extract lon/lat back from spatial object
coords <- crds(coords_vect)

df_points$lon <- coords[, 1]
df_points$lat <- coords[, 2]

# ---- STEP 9: Add Log-transformed Flux ----
df_points <- df_points %>%
  mutate(log_flux = log(flux + 1e-3))  # Add epsilon to avoid log(0)

# ---- STEP 10: Prepare Enhanced Modelling Dataset ----
predictor_vars_extended <- c("temperature", "salinity", "oxygen", "nitrate", 
                             "phosphate", "chlorophyll", "density", "depth", 
                             "lat", "lon")
target_var_log <- "log_flux"

# Only keep complete cases
model_data_extended <- df_points %>%
  dplyr::select(dplyr::all_of(c(target_var_log, predictor_vars_extended))) %>%
  filter(if_all(everything(), ~ !is.na(.)))

# ---- STEP 11: Split into Training and Test Sets ----
set.seed(123)
train_indices_ext <- sample(seq_len(nrow(model_data_extended)), size = 0.8 * nrow(model_data_extended))
train_data_ext <- model_data_extended[train_indices_ext, ]
test_data_ext  <- model_data_extended[-train_indices_ext, ]

# ---- STEP 12: Train Random Forest on Log-Flux ----
rf_model_ext <- randomForest(
  log_flux ~ ., 
  data = train_data_ext, 
  ntree = 500,
  importance = TRUE
)

# ---- STEP 13: Evaluate Enhanced Model ----
pred_log_flux <- predict(rf_model_ext, newdata = test_data_ext)
pred_flux <- exp(pred_log_flux)  # Back-transform

# Compare to actual (non-log) flux
actual_flux <- exp(test_data_ext$log_flux)
r_squared_ext <- cor(pred_flux, actual_flux)^2
rmse_ext <- sqrt(mean((pred_flux - actual_flux)^2))

cat("\n===== Enhanced Model Evaluation =====\n")
cat("Test R²:", round(r_squared_ext, 3), "\n")
cat("Test RMSE:", round(rmse_ext, 3), "\n")

mae <- mean(abs(pred_flux - actual_flux))
bias <- mean(pred_flux - actual_flux)

# Calibration slope (regression)
cal_model <- lm(actual_flux ~ pred_flux)
slope <- coef(cal_model)[2]

cat("MAE:", round(mae, 3), "\n")
cat("Bias:", round(bias, 3), "\n")
cat("Calibration Slope:", round(slope, 3), "\n")

# ---- STEP 14: Plot Variable Importance ----
varImpPlot(rf_model_ext, main = "Variable Importance (Extended RF Model)")

# ---- STEP 15: Predicted vs Observed Plot (Back-transformed) ----
plot(actual_flux, pred_flux,
     xlab = "Observed Flux (µmol m⁻² d⁻¹)",
     ylab = "Predicted Flux",
     main = "Observed vs Predicted Flux (RF, log-transformed)",
     pch = 16, col = "darkgreen")
abline(0, 1, col = "red", lwd = 2)


# ---- STEP 16: Residual Analysis and Mapping ----
# Calculate residuals (observed - predicted)
test_data_ext$pred_flux <- pred_flux
test_data_ext$actual_flux <- actual_flux
test_data_ext$residual <- test_data_ext$actual_flux - test_data_ext$pred_flux

# inspect residual distribution
# hist(test_data_ext$residual, breaks = 50, main = "Residuals (Observed - Predicted)", xlab = "Residual")

# Plot residuals vs predicted
plot(test_data_ext$pred_flux, test_data_ext$residual,
     xlab = "Predicted Flux", ylab = "Residual",
     main = "Residuals vs Predicted (log-transformed RF)",
     pch = 16, col = "darkred")
abline(h = 0, col = "blue", lwd = 2)


model_data_extended <- df_points %>%
  mutate(log_flux = log(flux)) %>%
  dplyr::select(lat, lon, all_of(c(target_var_log, predictor_vars_extended))) %>%
  filter(if_all(everything(), ~ !is.na(.))) %>%
  mutate(row_id = row_number())  # Store original row index

# Split again using new object
set.seed(123)
train_indices_ext <- sample(seq_len(nrow(model_data_extended)), size = 0.8 * nrow(model_data_extended))
train_data_ext <- model_data_extended[train_indices_ext, ]
test_data_ext  <- model_data_extended[-train_indices_ext, ]

# Predict and calculate residuals again
pred_log_flux <- predict(rf_model_ext, newdata = test_data_ext)
pred_flux <- exp(pred_log_flux)
actual_flux <- exp(test_data_ext$log_flux)
residuals <- actual_flux - pred_flux
test_data_ext$residual <- residuals

# Create a spatial sf object for residuals
residual_map_sf <- st_as_sf(
  test_data_ext,
  coords = c("lon", "lat"),
  crs = 4326
)
world <- map_data("world")

# Set symmetric limits for diverging color scale
max_resid <- max(abs(test_data_ext$residual), na.rm = TRUE)

# Download the land polygons (only needs to run once)
world <- ne_countries(scale = "medium", returnclass = "sf")

# Plot residuals with land background
ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey70") +  # Land
  geom_sf(data = residual_map_sf, aes(color = residual), size = 1) +
  scale_color_gradient2(
    low = "blue", mid = "grey", high = "red",
    midpoint = 0,
    name = "Residual\n(Observed - Predicted)"
  ) +
  coord_sf(xlim = c(-40, 40), ylim = c(20, 75), expand = FALSE) +  # Zoom in to Europe region
  theme_minimal() +
  labs(
    title = "Spatial Distribution of N2O Flux Prediction Residuals",
    subtitle = "Europe + Surrounding Oceans (Blue = Overprediction, Red = Underprediction)"
  )


# ---------------------------------------------------------------------- # 
# ---------------------------------------------------------------------- # 
# ---- STEP A: Comparing to Yang et al., 2020s Predicted N2O Flux
# ---- STEP A1: Loading and Extracting Yang et al.,2020s Predicted Flux ----
yang_file <- "n2oFlux-Yang2020.nc"
yang_flux_raster <- raster(yang_file, varname = "n2oFlux_EnsMean_g-pm2-pyr")  # Ensemble mean

# Convert from g N m⁻² yr⁻¹ to µmol N₂O m⁻² d⁻¹
# Molar mass N₂O = 44 g/mol → 1 g N₂O = (1 / 44) mol = (1 / 44) * 1e6 µmol
yang_flux_raster_umol <- yang_flux_raster * (1e6 / 44) / 365

# Extract Yang’s prediction at your points
coords_sf <- st_as_sf(df_points, coords = c("lon", "lat"), crs = 4326)
df_points$flux_yang <- raster::extract(yang_flux_raster_umol, coords_sf)

summary(df_points$flux_yang)

# ---- STEP A2: Correlate Observed and Yang et al.,2020s Modeled Flux
cor_obs_yang <- cor(df_points$flux, df_points$flux_yang, use = "complete.obs")
cat("Correlation between observed and Yang's flux:", round(cor_obs_yang, 3), "\n")

ggplot(df_points, aes(x = flux, y = flux_yang)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", color = "red") +
  labs(
    title = "Observed N2O Flux vs Yang et al. 2020 Prediction",
    x = "Observed Flux (µmol m⁻² d⁻¹)",
    y = "Yang et al. Predicted Flux (µmol m⁻² d⁻¹)"
  ) +
  theme_minimal()

 # ---- STEP A3: Correlate my RF Predictions vs Yang et al.,2020
# First: copy over lat/lon into test_data_ext
test_data_ext <- test_data_ext %>%
  mutate(
    lat = df_points$lat[as.numeric(rownames(test_data_ext))],
    lon = df_points$lon[as.numeric(rownames(test_data_ext))]
  )

# Then: extract Yang’s predicted flux at test points
test_coords_sf <- st_as_sf(test_data_ext, coords = c("lon", "lat"), crs = 4326)
test_data_ext$flux_yang <- raster::extract(yang_flux_raster_umol, test_coords_sf)

# Summary check
summary(test_data_ext$flux_yang)

# Compute correlation
cor_rf_yang <- cor(pred_flux, test_data_ext$flux_yang, use = "complete.obs")
cat("Correlation between RF predicted flux and Yang's:", round(cor_rf_yang, 3), "\n")

# Plot observed vs Yang
ggplot(test_data_ext, aes(x = pred_flux, y = flux_yang)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "blue") +
  labs(
    title = "My RF Predicted Flux vs Yang et al. 2020 Prediction",
    x = "My Predicted Flux (µmol m⁻² d⁻¹)",
    y = "Yang Predicted Flux (µmol m⁻² d⁻¹)"
  ) +
  theme_minimal()


# ---- STEP B: Use Yang et al.,2020s Flux Predictions as a Predictor in my own Model
# ---- STEP B1: Adding Yang et al.,2020 Flux to my Dataset
# Convert df_points to sf
coords_sf_all <- st_as_sf(df_points, coords = c("lon", "lat"), crs = 4326)

# Extract Yang’s flux for all points
df_points$flux_yang <- raster::extract(yang_flux_raster_umol, coords_sf_all)

# Log-transform again (if needed for model)
df_points$log_flux <- log(df_points$flux)


# ---- STEP B2: Building a New Model that Includes flux_yang
# Define updated predictors
predictors_with_yang <- c("temperature", "salinity", "oxygen", "nitrate", 
                          "phosphate", "chlorophyll", "density", "depth", 
                          "lat", "lon", "flux_yang")

# Only complete cases
df_model_with_yang <- df_points %>%
  dplyr::select(log_flux, dplyr::all_of(predictors_with_yang)) %>%
  dplyr::filter(dplyr::if_all(everything(), ~ !is.na(.)))

# Split again
set.seed(123)
train_idx_yang <- sample(seq_len(nrow(df_model_with_yang)), size = 0.8 * nrow(df_model_with_yang))
train_yang <- df_model_with_yang[train_idx_yang, ]
test_yang  <- df_model_with_yang[-train_idx_yang, ]

# Train new RF model
rf_yang_model <- randomForest(log_flux ~ ., data = train_yang, ntree = 500, importance = TRUE)

# Predict and back-transform
pred_log_flux_yang <- predict(rf_yang_model, newdata = test_yang)
pred_flux_yang <- exp(pred_log_flux_yang)
actual_flux_yang <- exp(test_yang$log_flux)

# Evaluate
r_squared_yang <- cor(pred_flux_yang, actual_flux_yang)^2
rmse_yang <- sqrt(mean((pred_flux_yang - actual_flux_yang)^2))

cat("\n===== Enhanced Model w/ Yang Flux Included =====\n")
cat("Test R²:", round(r_squared_yang, 3), "\n")
cat("Test RMSE:", round(rmse_yang, 3), "\n")


# ---- STEP B3: Visualize Importance ----
# Variable importance
varImpPlot(rf_yang_model, main = "Variable Importance (RF + Yang Flux)")

# Observed vs Predicted
plot(actual_flux_yang, pred_flux_yang,
     xlab = "Observed Flux", ylab = "Predicted Flux",
     main = "Observed vs Predicted (RF + Yang Flux)",
     pch = 16, col = "darkblue")
abline(0, 1, col = "red", lwd = 2)


# ---- STEP C: Add Biome Classification from Yang et al.,2020 ----
# ---- STEP C1: Load Yang’s Biome Raster and Extract Biome Labels
# Load Yang et al. biome raster from NetCDF
yang_file <- "dn2o-mapped-Yang2020.nc"

# Extract biome classification (1 to 6 as per their legend)
biome_raster <- raster(yang_file, varname = "biomes_masks")

# Extract biome category at your observation points
coords_sf_biome <- st_as_sf(df_points, coords = c("lon", "lat"), crs = 4326)
df_points$biome <- raster::extract(biome_raster, coords_sf_biome)

# Convert biome to factor
df_points$biome <- factor(df_points$biome,
                          levels = 1:6,
                          labels = c("Tropical", "Upwelling", "Polar", "Mid-lat", "DeepMix", "SubtropicalGyre"))


# ---- STEP C2: Add Biome to Model and Train Again
# Add biome to your list of predictors
predictors_biome <- c(predictors_with_yang, "biome")  # All previous predictors + biome

# Filter complete cases
df_model_biome <- df_points %>%
  dplyr::select(log_flux, dplyr::all_of(predictors_biome)) %>%
  dplyr::filter(dplyr::if_all(everything(), ~ !is.na(.)))

# Split
set.seed(123)
train_idx_biome <- sample(seq_len(nrow(df_model_biome)), size = 0.8 * nrow(df_model_biome))
train_biome <- df_model_biome[train_idx_biome, ]
test_biome  <- df_model_biome[-train_idx_biome, ]

# Train model
rf_model_biome <- randomForest(log_flux ~ ., data = train_biome, ntree = 500, importance = TRUE)

# Evaluate
pred_biome_log <- predict(rf_model_biome, newdata = test_biome)
pred_biome <- exp(pred_biome_log)
actual_biome <- exp(test_biome$log_flux)

cat("\n===== Model with Biome Included =====\n")
cat("Test R²:", round(cor(pred_biome, actual_biome)^2, 3), "\n")
cat("Test RMSE:", round(sqrt(mean((pred_biome - actual_biome)^2)), 3), "\n")

# Optional: Variable importance
varImpPlot(rf_model_biome, main = "RF Variable Importance (Biome Included)")


# ---- STEP D: Analyze Residuals by Biome or Region
# ---- STEP D1: Compute Residuals and Attach to Biomes
test_biome$residual <- actual_biome - pred_biome
test_biome$biome <- test_biome$biome  # Already present from training split

# ---- STEP D2: Boxplot of Residuals by Biome
ggplot(test_biome, aes(x = biome, y = residual)) +
  geom_boxplot(fill = "skyblue") +
  labs(title = "Residuals by Ocean Biome", x = "Biome", y = "Residual (Observed - Predicted)") +
  theme_minimal()

# ---- STEP D3: Spatial Bias Map # This will be left out for the article/ just something Yang 2020 did
# Start from df_model_biome, which matches test_biome row-for-row
res_coords <- df_model_biome[-train_idx_biome, ] %>%
  dplyr::select(lat, lon) %>%
  mutate(residual = test_biome$residual)

res_coords_europe <- res_coords %>%
  filter(lon >= -30, lon <= 60,
         lat >= 30, lat <= 75)

ggplot(res_coords_europe, aes(x = lon, y = lat, color = residual)) +
  geom_point(size = 1.5) +
  scale_color_gradient2(
    low = "blue", mid = "lightyellow", high = "red", midpoint = 0,
    name = "Residual (µmol m⁻² d⁻¹)"
  ) +
  borders("world", colour = "gray60", fill = "gray90") +
  coord_quickmap(xlim = c(-30, 60), ylim = c(30, 75)) +
  theme_minimal() +
  labs(
    title = "Geographic Residuals over Europe",
    x = "Longitude",
    y = "Latitude"
  )



# ---- STEP E: Making Gridded Predictions ----
# ---- STEP E1: Stack Gridded Predictors into a Raster Brick
# 1. Load raw raster layers (keep only surface = depth 0)
temperature_r <- rast("woa2023_Temperature.nc", subds = "t_an")
temperature_r <- temperature_r[[grep("depth=0", names(temperature_r))]]

salinity_r <- rast("woa2023_Salinity.nc", subds = "s_an")
salinity_r <- salinity_r[[grep("depth=0", names(salinity_r))]]

oxygen_r <- rast("woa2023_DisOxygen.nc", subds = "o_an")
oxygen_r <- oxygen_r[[grep("depth=0", names(oxygen_r))]]

nitrate_r <- rast("woa2023_Nitrate.nc", subds = "n_an")
nitrate_r <- nitrate_r[[grep("depth=0", names(nitrate_r))]]

phosphate_r <- rast("woa2023_Phosphate.nc", subds = "p_an")
phosphate_r <- phosphate_r[[grep("depth=0", names(phosphate_r))]]

density_r <- rast("woa2018_Density.nc", subds = "I_an")
density_r <- density_r[[grep("depth=0", names(density_r))]]

depth_r <- rast("GEBCO_2025_sub_ice.nc")  # already 1 layer

# 2. Load and average chlorophyll (already surface-only)
chl_file1 <- "cmems_mod_glo_bgc_my_0.25deg_P1M-m_1760291720536.nc"
chl_file2 <- "cmems_mod_glo_bgc_myint_0.25deg_P1M-m_1760291822771.nc"

r1 <- rast(chl_file1)
r2 <- rast(chl_file2)

surface_names_r1 <- grep("chl_depth=0.50576", names(r1), value = TRUE)
surface_names_r2 <- grep("chl_depth=0.50576", names(r2), value = TRUE)

chl_stack1 <- r1[[surface_names_r1]]
chl_stack2 <- r2[[surface_names_r2]]

chl_full_stack <- c(chl_stack1, chl_stack2)
chlorophyll_r <- mean(chl_full_stack, na.rm = TRUE)

# 3. Crop to European extent
extent_europe <- ext(-35, 50, 20, 80)
rasters_uncropped <- list(
  temperature_r, salinity_r, oxygen_r, nitrate_r,
  phosphate_r, chlorophyll_r, density_r, depth_r
)
raster_list_cropped <- lapply(rasters_uncropped, function(r) crop(r, extent_europe))

# 4. Align CRS
for (i in seq_along(raster_list_cropped)) {
  crs(raster_list_cropped[[i]]) <- "EPSG:4326"
}

# 5. Align grids
ref_raster <- raster_list_cropped[[1]]
raster_list_aligned <- lapply(raster_list_cropped, function(r)
  resample(r, ref_raster, method = "bilinear")
)

# 6. Stack aligned rasters
predictor_stack <- rast(raster_list_aligned)
names(predictor_stack) <- c("temperature", "salinity", "oxygen", "nitrate",
                            "phosphate", "chlorophyll", "density", "depth")

# 7. Check
predictor_stack
plot(predictor_stack)


# ---- STEP E2: Predict N₂O Flux Across Europe Using RF Model ----
# Add Latitude and Longitude as Raster Layers ---
# --- Rebuild latitude and longitude rasters correctly ---
# Each one should have ONE layer only, aligned to the predictor grid
lat_raster <- rast(predictor_stack[[1]])  # copy geometry
lon_raster <- rast(predictor_stack[[1]])

# Fill with coordinate values
lat_raster[] <- yFromCell(lat_raster, 1:ncell(lat_raster))
lon_raster[] <- xFromCell(lon_raster, 1:ncell(lon_raster))

# Check that they’re single-layer rasters
nlyr(lat_raster)
nlyr(lon_raster)  

# Combine everything
predictor_stack_full <- c(predictor_stack, lat_raster, lon_raster)

# Verify number of layers (should be 10)
nlyr(predictor_stack_full)

names(predictor_stack_full) <- c("temperature", "salinity", "oxygen", "nitrate",
                                 "phosphate", "chlorophyll", "density", "depth",
                                 "lat", "lon")

# Predict flux
predicted_flux <- terra::predict(
  object = predictor_stack_full,
  model  = rf_model,
  na.rm  = TRUE,
  progress = TRUE
)

# Save and plot
names(predicted_flux) <- "Predicted_N2O_Flux"
writeRaster(predicted_flux, "Predicted_N2O_Flux_Europe.tif", overwrite = TRUE)

plot(predicted_flux, col = viridis::viridis(100),
     main = expression("Predicted N"[2]*"O Flux (µmol m"^{-2}*" d"^{-1}*")"))

# Compute total European coastal N2O emission
# ============================================================
# 1. Compute area of each grid cell (in m^2)
#    (terra accounts for varying cell size with latitude)
cell_area_m2 <- terra::cellSize(predicted_flux, unit = "m")

# 2. Convert flux from µmol N2O m^-2 d^-1 to g N2O m^-2 yr^-1
#    µmol -> mol: 1e-6
#    mol -> g N2O: * 44 g/mol
#    per day -> per year: * 365
flux_gN2O_m2_yr <- predicted_flux * 1e-6 * 44 * 365

# 3. Total N2O emission per grid cell (g N2O yr^-1)
cell_emission_gN2O_yr <- flux_gN2O_m2_yr * cell_area_m2

# 4. Sum over all coastal grid cells to get total European emission
g_sum <- terra::global(cell_emission_gN2O_yr,
                       fun = "sum",
                       na.rm = TRUE)

# g_sum is a data.frame with one row and one column (e.g. column "sum")
total_emission_gN2O_yr <- as.numeric(g_sum[1, 1])

# 5. Convert to Tg N2O yr^-1 and also to Tg N yr^-1 
total_emission_TgN2O_yr <- total_emission_gN2O_yr / 1e12  # g -> Tg


# N2O has 28 g N per 44 g N2O
total_emission_TgN_yr <- total_emission_TgN2O_yr * (28 / 44)

cat("\n===== Integrated European Coastal Emissions =====\n")
cat("Total N2O emission:", round(total_emission_TgN2O_yr, 4),
    "Tg N2O yr^-1\n")
cat("Total N2O-N emission:", round(total_emission_TgN_yr, 4),
    "Tg N yr^-1\n")

# Mask land using bathymetry
# ocean_mask <- depth_r < 0
# predicted_flux_ocean <- mask(predicted_flux, ocean_mask)

# Convert to data frame for ggplot
# pred_df <- as.data.frame(predicted_flux_ocean, xy = TRUE)
# names(pred_df)[3] <- "flux"
# pred_df <- pred_df[!is.na(pred_df$flux), ]

# Plot with coastlines
# land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

#ggplot() + geom_raster(data = pred_df, aes(x = x, y = y, fill = flux)) +
# geom_sf(data = land, fill = "gray90", color = "gray50", size = 0.2) +
# scale_fill_viridis(
#   name = expression("Predicted N"[2]*"O Flux (µmol m"^{-2}*" d"^{-1}*")"),
#   option = "C"
# ) +
# coord_sf(xlim = c(-30, 60), ylim = c(30, 75), expand = FALSE) +
# theme_minimal() +
#   labs(
#   title = expression("Predicted N"[2]*"O Flux Across Europe"),
#    x = "Longitude",
#    y = "Latitude"
#  ) +
#  theme(plot.title = element_text(hjust = 0.5, face = "bold")) 


# ---- Compare with Yang et al. (2020) Gridded Flux -----
# Load Yang et al. (2020) flux (g N2O m⁻² yr⁻¹)
yang_flux_raster <- terra::rast("n2oFlux-Yang2020.nc", subds = "n2oFlux_EnsMean_g-pm2-pyr")

# Crop to your study extent
yang_flux_cropped <- terra::crop(yang_flux_raster, predicted_flux)

# Resample to match your model grid
yang_flux_resampled <- terra::resample(yang_flux_cropped, predicted_flux, method = "bilinear")

# Compute mean across 12 layers (if monthly)
yang_flux_mean <- mean(yang_flux_resampled, na.rm = TRUE)

# Convert to µmol N2O m⁻² d⁻¹
# 1 g N₂O = (1 / 44) mol N₂O = (1e6 / 44) µmol N₂O
# Divide by 365 for per day
yang_flux_mean_umol <- yang_flux_mean * (1e6 / (44 * 365))
names(yang_flux_mean_umol) <- "Yang_N2O_Flux_umol_m2_d"

# Plot Yang mean
plot(yang_flux_mean_umol, main = "Yang et al. (2020) Mean N₂O Flux",
     col = viridis::viridis(100))

# ---- Compare with your predicted flux ----
# Extract paired values
vals <- as.data.frame(c(predicted_flux, yang_flux_mean_umol))
names(vals) <- c("your_flux", "yang_flux")
vals <- vals[complete.cases(vals), ]

# Summary and correlation
cat("Number of valid cells:", nrow(vals), "\n")
cat("Correlation between your model and Yang (2020):",
    round(cor(vals$your_flux, vals$yang_flux), 3), "\n")

# ---- Plot comparison ----
par(mfrow = c(1, 3))
plot(predicted_flux, main = "My Predicted N₂O Flux",
     col = viridis::viridis(100))
plot(yang_flux_mean_umol, main = "Yang et al. (2020) N₂O Flux",
     col = viridis::viridis(100))
plot(predicted_flux - yang_flux_mean_umol, main = "Difference (My - Yang)",
     col = viridis::viridis(100))

# (a) Your predicted flux
terra::plot(
  predicted_flux,
  main = "a) This Study's Predicted Flux",
  col  = viridis::viridis(100),
  range = c(0, 0.6)          # << key change
)

# (b) Yang et al. (2020)
terra::plot(
  yang_flux_mean_umol,
  main = "b) Yang et al. (2020) Flux",
  col  = viridis::viridis(100),
  range = c(0, 2)          # << same scale
)

# (c) Difference (can keep its own scale)
terra::plot(
  predicted_flux - yang_flux_mean_umol,
  main = "c) Difference",
  col  = viridis::viridis(100)
)

# ---------------------------------------------------------------
# ---------------------------------------------------------------
# Fix the plor rows again
par(mfrow = c(1,1))

# ---- STEP 17: RandomForest Hyperparameter Tuning (Same model but trying a better configuration) ----
# Define tuning grid
rf_grid <- expand.grid(mtry = c(2, 4, 6, 8, 10))

# Train model using caret
ctrl <- trainControl(method = "cv", number = 5)  # 5-fold CV

rf_tuned <- train(
  log_flux ~ ., 
  data = train_data_ext, 
  method = "rf",
  trControl = ctrl,
  tuneGrid = rf_grid,
  ntree = 500,
  importance = TRUE
)

print(rf_tuned)
plot(rf_tuned)

# ---- STEP 18: Train Basic XGBoost Model (Using predictors available for grid) ----
#  Define matching predictors (for both model and gridded stack)
predictors_basic <- c("lat", "lon", "temperature", "salinity", "oxygen", 
                      "nitrate", "phosphate", "chlorophyll", "density", "depth")

#  Prepare cleaned input data with only these predictors
xgb_data_basic <- df_points %>%
  mutate(log_flux = log(flux + 1e-3)) %>%  # re-define here if needed
  dplyr::select(log_flux, dplyr::all_of(predictors_basic)) %>%
  filter(if_all(everything(), ~ !is.na(.)))

#  Extract predictor matrix and target vector
X_basic <- as.matrix(xgb_data_basic[, predictors_basic])  # predictor matrix
y_basic <- xgb_data_basic$log_flux                        # target vector

#  Add explicit column names to matrix BEFORE DMatrix
colnames(X_basic) <- predictors_basic

#  Train/test split
set.seed(42)
train_idx_basic <- sample(seq_len(nrow(X_basic)), size = 0.8 * nrow(X_basic))
X_train_basic <- X_basic[train_idx_basic, ]
y_train_basic <- y_basic[train_idx_basic]
X_test_basic  <- X_basic[-train_idx_basic, ]
y_test_basic  <- y_basic[-train_idx_basic]

#  DMatrix for XGBoost
dtrain_basic <- xgb.DMatrix(data = X_train_basic, label = y_train_basic)
dtest_basic  <- xgb.DMatrix(data = X_test_basic, label = y_test_basic)

#  Train basic XGBoost model
xgb_model_basic <- xgboost(
  data = dtrain_basic,
  nrounds = 200,
  objective = "reg:squarederror",
  eval_metric = "rmse",
  verbose = 0
)

#  Predict and evaluate
pred_log_basic <- predict(xgb_model_basic, newdata = dtest_basic)
pred_flux_basic <- exp(pred_log_basic)
actual_flux_basic <- exp(y_test_basic)

#  Evaluation metrics
r2_basic <- cor(pred_flux_basic, actual_flux_basic)^2
rmse_basic <- sqrt(mean((pred_flux_basic - actual_flux_basic)^2))
cat("Basic XGBoost Model R²:", round(r2_basic, 3), "\n")
cat("Basic XGBoost Model RMSE:", round(rmse_basic, 3), "\n")



# ---- STEP 19: Gridded Prediction Using Basic XGBoost Model ----
# 1. Convert raster stack to dataframe with coordinates
grid_df_basic <- as.data.frame(predictor_stack_full, xy = TRUE, na.rm = TRUE)

# 2. Ensure proper column names
names(grid_df_basic)[1:2] <- c("lon", "lat")

# 3.Drop missing values
grid_df_basic <- grid_df_basic[complete.cases(grid_df_basic), ]

# 4. Use the exact same predictors as training
predictors_basic <- c("lat", "lon", "temperature", "salinity", "oxygen", 
                      "nitrate", "phosphate", "chlorophyll", "density", "depth")

# 5. Subset and reorder columns
grid_matrix_basic <- as.matrix(grid_df_basic[, predictors_basic])
colnames(grid_matrix_basic) <- predictors_basic  # match training model

# 6. Convert to DMatrix for XGBoost
grid_dmatrix_basic <- xgb.DMatrix(data = grid_matrix_basic)

# 7. Predict and back-transform
grid_log_pred_basic <- predict(xgb_model_basic, newdata = grid_dmatrix_basic)
grid_flux_pred_basic <- exp(grid_log_pred_basic)


# 8. Add predictions back
grid_df_basic$flux_xgb_basic <- grid_flux_pred_basic

# 9. Convert to raster for visualization
raster_flux_basic <- terra::rast(
  terra::ext(predictor_stack_full),
  nrows = nrow(predictor_stack_full),
  ncols = ncol(predictor_stack_full),
  crs   = crs(predictor_stack_full)
)

raster_flux_basic <- terra::rast(
  data.frame(x = grid_df_basic$lon, 
             y = grid_df_basic$lat, 
             flux_xgb_basic = grid_df_basic$flux_xgb_basic),
  crs = crs(predictor_stack_full)
)

names(raster_flux_basic) <- "flux_xgb_basic"

# 10. Plot
plot(raster_flux_basic,
     main = "Gridded N₂O Flux Prediction (Basic XGBoost)",
     col = viridis::viridis(100))

# Formatting
# Remove duplicated columns (e.g. keep only "x" and "y")
grid_df_basic <- grid_df_basic[, !duplicated(names(grid_df_basic))]

# Optionally rename x/y for clarity
# grid_df_basic <- grid_df_basic %>%
#  rename(lon = x, lat = y)
# Remove duplicate column names
grid_df_basic <- grid_df_basic[, !duplicated(names(grid_df_basic))]

# Double-check names
names(grid_df_basic)

# ggplot map
land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  geom_raster(data = grid_df_basic, aes(x = lon, y = lat, fill = flux_xgb_basic)) +
  geom_sf(data = land, fill = "gray90", color = "gray50") +
  scale_fill_viridis(
    name = expression("Predicted N"[2]*"O Flux (µmol m"^{-2}*" d"^{-1}*")"),
    option = "C"
  ) +
  coord_sf(xlim = c(-30, 60), ylim = c(30, 75), expand = FALSE) +
  theme_minimal() +
  labs(
    title = "Gridded N2O Flux Prediction (Basic XGBoost)",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

# ---- STEP 20: Train Tuned XGBoost Model (Better Hyperparameters) ----
# 1. Define same predictors used earlier
predictors_basic <- c("lat", "lon", "temperature", "salinity", "oxygen", 
                      "nitrate", "phosphate", "chlorophyll", "density", "depth")

# 2. Prepare cleaned dataset 
xgb_data <- df_points %>%
  mutate(log_flux = log(flux + 1e-3)) %>%
  dplyr::select(log_flux, dplyr::all_of(predictors_basic)) %>%
  filter(if_all(everything(), ~ !is.na(.)))

# 3. Train/test split
set.seed(123)
train_idx_xgb <- sample(seq_len(nrow(xgb_data)), size = 0.8 * nrow(xgb_data))
train_xgb <- xgb_data[train_idx_xgb, ]
test_xgb  <- xgb_data[-train_idx_xgb, ]

# 4. Convert to matrix format
dtrain <- xgb.DMatrix(data = as.matrix(train_xgb[, predictors_basic]), label = train_xgb$log_flux)
dtest  <- xgb.DMatrix(data = as.matrix(test_xgb[, predictors_basic]), label = test_xgb$log_flux)

# 5. Define tuned parameters 
params <- list(
  objective = "reg:squarederror",
  eval_metric = "rmse",
  eta = 0.05,            # Learning rate (lower = slower, but more stable)
  max_depth = 6,         # Max depth of trees
  subsample = 0.8,       # Subsampling of rows
  colsample_bytree = 0.8 # Subsampling of features
)

# 6. Train tuned XGBoost model
xgb_model_tuned <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 300,   # Slightly higher than basic model
  verbose = 0
)

# 7. Predict and evaluate
pred_log_tuned <- predict(xgb_model_tuned, newdata = dtest)
pred_flux_tuned <- exp(pred_log_tuned)
actual_flux_tuned <- exp(test_xgb$log_flux)

# 8. Evaluate model
r2_tuned <- cor(pred_flux_tuned, actual_flux_tuned)^2
rmse_tuned <- sqrt(mean((pred_flux_tuned - actual_flux_tuned)^2))

cat("Tuned XGBoost Model R²:", round(r2_tuned, 3), "\n")
cat("Tuned XGBoost Model RMSE:", round(rmse_tuned, 3), "\n")

# 9. Plot observed vs predicted
plot(actual_flux_tuned, pred_flux_tuned,
     xlab = "Observed Flux", ylab = "Predicted Flux",
     main = "Observed vs Predicted (Tuned XGBoost)",
     pch = 16, col = "darkorange")
abline(0, 1, col = "red", lwd = 2)

# 10. Variable Importance
importance_tuned <- xgb.importance(model = xgb_model_tuned)
xgb.plot.importance(importance_matrix = importance_tuned, top_n = 10)


# ---- STEP 21: Gridded Prediction Using Tuned XGBoost Model ----
# Convert predictor stack (with lat/lon included) to dataframe
grid_df_tuned <- as.data.frame(predictor_stack_full, xy = TRUE)

# Remove duplicate names if they exist
grid_df_tuned <- grid_df_tuned[, !duplicated(names(grid_df_tuned))]

#    Rename x/y to lon/lat if needed
names(grid_df_tuned)[names(grid_df_tuned) == "x"] <- "lon"
names(grid_df_tuned)[names(grid_df_tuned) == "y"] <- "lat"

# Drop rows with missing values
grid_df_tuned <- grid_df_tuned[complete.cases(grid_df_tuned), ]

# Define predictor list (to match training model)
predictors_basic <- c("lat", "lon", "temperature", "salinity", "oxygen", 
                      "nitrate", "phosphate", "chlorophyll", "density", "depth")

# Verify that all predictors exist before selecting
missing_cols <- setdiff(predictors_basic, names(grid_df_tuned))
if (length(missing_cols) > 0) {
  stop(paste("Missing columns in grid_df_tuned:", paste(missing_cols, collapse = ", ")))
}

# Extract matrix for prediction
grid_matrix_tuned <- as.matrix(grid_df_tuned[, predictors_basic])
colnames(grid_matrix_tuned) <- predictors_basic

# Convert to DMatrix and predict
grid_dmatrix_tuned <- xgb.DMatrix(data = grid_matrix_tuned)
grid_log_flux_tuned <- predict(xgb_model_tuned, grid_dmatrix_tuned)
grid_flux_tuned <- exp(grid_log_flux_tuned)

# Add predictions to dataframe
grid_df_tuned$flux_xgb_tuned <- grid_flux_tuned

# Convert to raster
raster_flux_tuned <- rasterFromXYZ(grid_df_tuned[, c("lon", "lat", "flux_xgb_tuned")])

# Plot raster (basic)
plot(raster_flux_tuned,
     main = "Gridded N₂O Flux Prediction (Tuned XGBoost)",
     col = viridis::viridis(100))

# ---- CLEANING DUPLICATE COLUMNS ----
# Remove duplicate column names
grid_df_tuned <- grid_df_tuned[, !duplicated(names(grid_df_tuned))]

# Rename coordinate columns to match expectation (x → lon, y → lat)
if ("x" %in% names(grid_df_tuned)) names(grid_df_tuned)[names(grid_df_tuned) == "x"] <- "lon"
if ("y" %in% names(grid_df_tuned)) names(grid_df_tuned)[names(grid_df_tuned) == "y"] <- "lat"

# Confirm no more duplicates
anyDuplicated(names(grid_df_tuned))
names(grid_df_tuned)

# ggplot version (for clean map)
land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  geom_raster(data = grid_df_tuned, aes(x = lon, y = lat, fill = flux_xgb_tuned)) +
  geom_sf(data = land, fill = "gray90", color = "gray50") +
  scale_fill_viridis(
    name = expression("Predicted N"[2]*"O Flux (µmol m"^{-2}*" d"^{-1}*")"),
    option = "C"
  ) +
  coord_sf(xlim = c(-30, 60), ylim = c(30, 75), expand = FALSE) +
  theme_minimal() +
  labs(
    title = "Gridded N₂O Flux Prediction (Tuned XGBoost)",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

# -------------------------------------------
# Scenario Analysis: N₂O Flux Sensitivity
# Based on predictor_stack and XGBoost tuned model
# -------------------------------------------
# Clone base stack
stack_base <- predictor_stack

# Define scenario modifications
stack_temp2C <- stack_base
stack_temp2C[["temperature"]] <- stack_temp2C[["temperature"]] + 2

stack_oxygen20 <- stack_base
stack_oxygen20[["oxygen"]] <- stack_oxygen20[["oxygen"]] * 0.8

stack_nutrient20 <- stack_base
stack_nutrient20[["nitrate"]] <- stack_nutrient20[["nitrate"]] * 1.2
stack_nutrient20[["phosphate"]] <- stack_nutrient20[["phosphate"]] * 1.2

stack_multi <- stack_base
stack_multi[["temperature"]] <- stack_multi[["temperature"]] + 2
stack_multi[["oxygen"]] <- stack_multi[["oxygen"]] * 0.8
stack_multi[["nitrate"]] <- stack_multi[["nitrate"]] * 1.2
stack_multi[["phosphate"]] <- stack_multi[["phosphate"]] * 1.2

# Convert to prediction-ready format
to_grid_df <- function(stack) {
  df <- as.data.frame(stack, xy = TRUE)
  df <- df[complete.cases(df), ]
  df
}

df_base <- to_grid_df(stack_base)
df_temp <- to_grid_df(stack_temp2C)
df_o2   <- to_grid_df(stack_oxygen20)
df_nut  <- to_grid_df(stack_nutrient20)
df_multi <- to_grid_df(stack_multi)

# --- CLEAN COLUMN NAMES ---
fix_coords <- function(df) {
  # Remove duplicate names
  df <- df[, !duplicated(names(df))]
  
  # Rename coordinate columns
  if ("x" %in% names(df)) names(df)[names(df) == "x"] <- "lon"
  if ("y" %in% names(df)) names(df)[names(df) == "y"] <- "lat"
  
  return(df)
}

df_base  <- fix_coords(df_base)
df_temp  <- fix_coords(df_temp)
df_o2    <- fix_coords(df_o2)
df_nut   <- fix_coords(df_nut)
df_multi <- fix_coords(df_multi)

predict_flux <- function(df, predictors, model) {
  mat <- as.matrix(df[, predictors])
  colnames(mat) <- predictors
  pred <- predict(model, newdata = mat)
  exp(pred)
}


predictors <- predictors_basic  
df_base$flux <- predict_flux(df_base, predictors, xgb_model_tuned)
df_temp$flux <- predict_flux(df_temp, predictors, xgb_model_tuned)
df_o2$flux   <- predict_flux(df_o2, predictors, xgb_model_tuned)
df_nut$flux  <- predict_flux(df_nut, predictors, xgb_model_tuned)
df_multi$flux <- predict_flux(df_multi, predictors, xgb_model_tuned)

# Calculate deltas (flux difference from baseline)
df_temp$delta <- df_temp$flux - df_base$flux
df_o2$delta   <- df_o2$flux - df_base$flux
df_nut$delta  <- df_nut$flux - df_base$flux
df_multi$delta <- df_multi$flux - df_base$flux

# Visualize difference
plot_scenario <- function(df, title) {
  ggplot(df, aes(x = lon, y = lat, fill = delta)) +
    geom_raster() +
    geom_sf(data = land, inherit.aes = FALSE,
            fill = "gray90", color = "gray50", size = 0.1) +
    scale_fill_gradient2(
      low = "blue", mid = "grey", high = "red", midpoint = 0,
      name = expression(Delta * "N"[2]*"O Flux (µmol m"^{-2}*" d"^{-1}*")")
    ) +
    coord_sf(xlim = c(-30, 60), ylim = c(30, 75), expand = FALSE) +
    theme_minimal() +
    labs(title = title, x = "Longitude", y = "Latitude") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

plot_scenario(df_temp,  "Scenario: +2°C Temperature")
plot_scenario(df_o2,    "Scenario: -20% Oxygen")
plot_scenario(df_nut,   "Scenario: +20% Nutrients")
plot_scenario(df_multi, "Scenario: Combined Drivers")


plot_scenario <- function(df, title) {
  ggplot(df, aes(x = lon, y = lat, fill = delta)) +
    geom_raster() +
    geom_sf(data = land, inherit.aes = FALSE,
            fill = "gray90", color = "gray50", size = 0.1) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      limits = c(-0.5, 0.5),        # ADD A LIMIT
      oob = scales::squish,     
      midpoint = 0,
      name = expression(Delta * "N"[2]*"O Flux (µmol m"^{-2}*" d"^{-1}*")")
    ) +
    coord_sf(xlim = c(-30, 60), ylim = c(30, 75), expand = FALSE) +
    theme_minimal() +
    labs(title = title, x = "Longitude", y = "Latitude") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
}


# ----------------------------------------------------------------
# ---- Model Uncertainty: Bootstrap Random Forests ----
# =============================================================
#️ Parameters
n_bootstrap <- 30  # can increase to 100 if need be but it takes too long on my laptop
set.seed(42)

# ---- Ensure lat/lon layers exist ----
lat_raster <- rast(predictor_stack[[1]])
lon_raster <- rast(predictor_stack[[1]])
lat_raster[] <- yFromCell(lat_raster, 1:ncell(lat_raster))
lon_raster[] <- xFromCell(lon_raster, 1:ncell(lon_raster))

predictor_stack_full <- c(predictor_stack, lat_raster, lon_raster)
names(predictor_stack_full) <- c(
  "temperature", "salinity", "oxygen", "nitrate", "phosphate",
  "chlorophyll", "density", "depth", "lat", "lon"
)

# ---- Initialize storage ----
bootstrap_rasters <- list()

# ----Loop through bootstraps ----
for (i in 1:n_bootstrap) {
  cat("Training bootstrap model:", i, "\n")
  # Sample with replacement from training data
  boot_idx <- sample(nrow(train_data_ext), replace = TRUE)
  boot_data <- train_data_ext[boot_idx, ]
  # Remove row_id if it exists
  if ("row_id" %in% names(boot_data)) {
    boot_data <- dplyr::select(boot_data, -row_id)
  }
  # Train bootstrap Random Forest
  rf_boot <- randomForest(log_flux ~ ., data = boot_data, ntree = 500)
  # Predict log-flux using terra
  pred_log_boot <- terra::predict(
    predictor_stack_full,
    rf_boot,
    na.rm = TRUE,
    progress = FALSE
  )
  # Back-transform (log → original)
  flux_boot <- exp(pred_log_boot)
  # Store each bootstrap raster
  bootstrap_rasters[[i]] <- flux_boot
}

# ---- Combine bootstrap rasters ----
bootstrap_stack <- rast(bootstrap_rasters)

# ----Compute mean and standard deviation ----
flux_mean <- mean(bootstrap_stack, na.rm = TRUE)
flux_sd   <- app(bootstrap_stack, fun = sd, na.rm = TRUE)

# ---- Plot results ----
plot(flux_mean, main = "Bootstrap Mean N₂O Flux", col = viridis::viridis(100))
plot(flux_sd, main = "Bootstrap Uncertainty (SD)", col = viridis::magma(100))


# Convert to dataframe for ggplot
uncertainty_df <- as.data.frame(flux_sd, xy = TRUE)
names(uncertainty_df)[3] <- "uncertainty"
uncertainty_df <- na.omit(uncertainty_df)

# Land shapefile
land <- ne_countries(scale = "medium", returnclass = "sf")

# Plot
ggplot() +
  geom_raster(data = uncertainty_df, aes(x = x, y = y, fill = uncertainty)) +
  geom_sf(data = land, fill = "gray90", color = "gray50", size = 0.2) +
  scale_fill_viridis(
    name = expression("Prediction SD (µmol m"^{-2}*" d"^{-1}*")"),
    option = "C"
  ) +
  coord_sf(xlim = c(-30, 60), ylim = c(30, 75), expand = FALSE) +
  labs(
    title = "Prediction Uncertainty (Bootstrap SD)",
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# ---- Plot Prediction Uncertainty (from flux_sd) ----
# Convert to dataframe for plotting
uncertainty_df <- as.data.frame(flux_sd, xy = TRUE)
names(uncertainty_df)[3] <- "uncertainty"

# Remove NAs
uncertainty_df <- na.omit(uncertainty_df)

# Load land polygons for context
land <- ne_countries(scale = "medium", returnclass = "sf")

# Plot the uncertainty (standard deviation)
ggplot() +
  geom_raster(data = uncertainty_df, aes(x = x, y = y, fill = uncertainty)) +
  geom_sf(data = land, fill = "gray90", color = "gray50", size = 0.2) +
  scale_fill_viridis(
    name = expression("Prediction SD (µmol m"^{-2}*" d"^{-1}*")"),
    option = "C"
  ) +
  coord_sf(xlim = c(-30, 60), ylim = c(30, 75), expand = FALSE) +
  labs(
    title = "Prediction Uncertainty (Bootstrap SD)",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


# ---- Relative predictive uncertainty: SD / mean (in %) ----
eps <- 1e-6  # small number to avoid division by ~0
flux_mean_safe <- flux_mean
flux_mean_safe[abs(flux_mean_safe) < eps] <- NA  # ignore cells with ~0 mean flux

flux_rel_uncertainty     <- flux_sd / flux_mean_safe        # unitless
flux_rel_uncertainty_pct <- flux_rel_uncertainty * 100      # %

names(flux_rel_uncertainty_pct) <- "rel_uncertainty_pct"

summary(values(flux_rel_uncertainty_pct))

# Convert relative uncertainty to dataframe for ggplot
uncertainty_df <- as.data.frame(flux_rel_uncertainty_pct, xy = TRUE)
names(uncertainty_df)[3] <- "rel_uncertainty_pct"
uncertainty_df <- na.omit(uncertainty_df)

ggplot() +
  geom_raster(data = uncertainty_df,
              aes(x = x, y = y, fill = rel_uncertainty_pct)) +
  geom_sf(data = land, fill = "gray90", color = "gray50", size = 0.2) +
  scale_fill_viridis(
    name   = "Relative uncertainty (%)",
    option = "C",
    # you can tweak limits after seeing summary(), e.g.:
    # limits = c(0, 80),
    oob    = scales::squish
  ) +
  coord_sf(xlim = c(-30, 60), ylim = c(30, 75), expand = FALSE) +
  labs(
    title = "Prediction Uncertainty (SD / Mean, %)",
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# =================================================
# ---- MODEL INTERPRETABILITY ----
# =================================================
# ----- RANDOM FOREST -----
# If you trained rf_model (not rf_model_ext)
varImpPlot(rf_model, main = "Variable Importance (Random Forest)")

# ----- XGBOOST SHAP INTERPRETATION -----
# Prepare SHAP input data (same predictors used during training)
X_shap <- as.matrix(train_xgb[, predictors_basic])
colnames(X_shap) <- predictors_basic

# Compute SHAP values
shap_values <- shap.values(xgb_model_tuned, X_shap)

# Convert to long format and visualize summary
shap_long <- shap.prep(shap_contrib = shap_values$shap_score, X_train = X_shap)
shap.plot.summary(shap_long)

# ---- CLUSTERING OF OBSERVATION ENVIRONMENTS ----
df_cluster <- df_model %>%
  dplyr::select(predictor_vars) %>%
  scale() %>%
  as.data.frame()

# ---- K-MEANS CLUSTERING ----
set.seed(42)
k_clusters <- kmeans(df_cluster, centers = 4)  # Try 3–6 clusters

# Add cluster membership back to df_model
df_model$cluster <- as.factor(k_clusters$cluster)

# Visualize cluster centroids
fviz_cluster(k_clusters, data = df_cluster)

# Add spatial coordinates
df_model$lat <- df_points$lat[match(rownames(df_model), rownames(df_points))]
df_model$lon <- df_points$lon[match(rownames(df_model), rownames(df_points))]

# Plot clusters on map
ggplot(df_model, aes(x = lon, y = lat, color = cluster)) +
  geom_point(size = 2, alpha = 0.8) +
  borders("world", colour = "gray60", fill = "gray90") +
  coord_quickmap(xlim = c(-30, 60), ylim = c(30, 75)) +
  labs(
    title = "Clustered Environmental Profiles of Coastal Systems",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


library(pdp)

pdp_sal <- partial(rf_model, pred.var = "salinity")
plot(pdp_sal)

# Partial dependence for salinity
p_sal <- partial(
  rf_model,
  pred.var = "salinity",
  grid.resolution = 50,      # smooth curve
  train = train_data_ext     # your training data
)

ggplot(p_sal, aes(x = salinity, y = yhat)) +
  geom_line(size = 1.2, colour = "#1b9e77") +
  labs(
    x = "Salinity",
    y = "Partial dependence (log-flux)",
    title = "Partial Dependence of N₂O Flux on Salinity"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


vars <- c("salinity", "depth", "temperature ", "oxygen")

p_list <- lapply(vars, function(v) {
  pd <- partial(rf_model, pred.var = v, grid.resolution = 50)
  ggplot(pd, aes_string(x = v, y = "yhat")) +
    geom_line(size = 1.2, colour = "#1b9e77") +
    labs(x = v, y = "Partial dependence") +
    theme_minimal(base_size = 12)
})

library(patchwork)
wrap_plots(p_list, ncol = 2) +
  plot_annotation(
    title = "Partial Dependence Plots for Top Predictors"
  ) &
  theme(
    plot.title = element_text(hjust = 0.5)
  )

vars <- c("salinity", "depth", "temperature", "oxygen")

# Custom axis labels WITH UNITS
var_labels <- c(
  salinity    = "Salinity (psu)",
  depth       = "Depth (m)",
  temperature = "Temperature (°C)",
  oxygen      = "Oxygen (µmol/kg)"
)

p_list <- lapply(vars, function(v) {
  pd <- partial(rf_model, pred.var = v, grid.resolution = 50)
  
  ggplot(pd, aes_string(x = v, y = "yhat")) +
    geom_line(size = 1.2, colour = "#1b9e77") +
    labs(
      x = var_labels[[v]],      # ← uses the units here
      y = "Partial dependence"
    ) +
    theme_minimal(base_size = 12)
})

wrap_plots(p_list, ncol = 2) +
  plot_annotation(title = "Partial Dependence Plots for Top Predictors") &
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
