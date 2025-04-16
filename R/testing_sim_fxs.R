# Simulation of animal movement and midge distribution for disease risk assessment
# Using sf and terra packages as requested

# Load required packages
library(sf)
library(terra)
library(tidyterra)
library(ggplot2)
library(dplyr)
library(spsurvey)

set.seed(123)
# 5. Run the simulation -----------------

# Create the simulated landscape

# Test the landscape function:
study_area <- create_study_area()
str(study_area) # this is an sf object, just a polygon
plot(study_area)
st_crs(study_area)
st_area(study_area) # using that crs it has units in meters

# create the non-overlapping water bodies:
# values should be between 2-5
water_bodies <- create_water_bodies(study_area, n = 4)
ggplot() +
   geom_sf(data = study_area) +
   geom_sf(data = water_bodies)

# Create the non-overlapping feeders outside water bodies
feeders <- create_feeders(study_area, 10, water_bodies, 100, 100)
ggplot() +
   geom_sf(data = study_area) +
   geom_sf(data = water_bodies) +
   geom_sf(data = feeders)

# Create the environmental layers
env_rasters <- create_env_rasters(study_area, water_bodies)
plot(env_rasters)


# Simulate animal movement
animal_tracks <- simulate_animal_movement(
   study_area = study_area,
   water_bodies = water_bodies,
   feeders = feeders,
   n_animals = 5,
   n_steps = 200
)

animal_tracks |>
   # filter(animal_id == 1) |>
   ggplot() +
   geom_sf(data = study_area) +
   geom_sf(aes(color = factor(animal_id))) +
   geom_sf(data = water_bodies) +
   geom_sf(data = feeders)

# check UD
# Create animal utilization distribution
animal_ud <- create_animal_ud(animal_tracks, study_area, resolution = 10, smoothing_factor = 9)
plot(animal_ud)

ggplot() +
   geom_sf(data = study_area) +
   tidyterra::geom_spatraster(data = animal_ud) +
   scale_fill_viridis_c(name = "Utlization\nDistribution", option = "plasma") +
   # geom_sf(data = feeders, color = "white") +
   # geom_sf(data = water_bodies, fill = "grey80", alpha = 0.5) +
   theme_void() +
   labs(title = "Animal Utilization Distribution")

# Simulate midge data
midge_sim <- simulate_midge_data(env_rasters, water_bodies)
head(midge_sim)

midge_data <- midge_sim$midge_data
midge_data |>
   ggplot() +
   geom_sf(aes(color = factor(presence))) +
   geom_sf(data = water_bodies, fill = NA, color = "blue")

# Get a sample of 20-100 traps with the data for the sdm
midge_data_sample <- slice_sample(midge_data, n = 20)
midge_data_sample |>
   ggplot() +
   geom_sf(aes(color = factor(presence))) +
   geom_sf(data = water_bodies, fill = NA, color = "blue")

# Use a more complex sampling: GRTS
library(spsurvey)
# check that the data is sf object:
midge_data
grts_sample <- spsurvey::grts(
   midge_data,
   100
)

midge_data_grts <- grts_sample$sites_base |>
   dplyr::select(presence, elevation, water_dist, veg_index, temperature, geometry)

midge_data_grts|>
   ggplot() +
   geom_sf(aes(color = factor(presence))) +
   geom_sf(data = water_bodies, fill = NA, color = "blue")

# Fit midge distribution model
# midge_sdm <- fit_midge_sdm(midge_data_sample, env_rasters)
midge_sdm <- fit_midge_sdm(midge_data_grts, env_rasters)
midge_prediction <- midge_sdm$prediction
plot(midge_prediction)

# Plot midge distribution model
ggplot() +
   geom_sf(data = study_area) +
   tidyterra::geom_spatraster(data = midge_prediction) +
   geom_sf(data = water_bodies, alpha = 0.5) +
   scale_fill_distiller(palette = "BrBG") +
   # geom_sf(data = feeders, color = "white") +
   # geom_sf(data = water_bodies, fill = "grey80", alpha = 0.5) +
   theme_minimal() +
   # labs(title = "Midge probability of occurrence")
   # scale_fill_viridis_c(name = "Midge\nProbability", option = "plasma") +
   scale_shape_manual(values = c(4, 16), name = "Midge Presence") +
   theme_void() +
   labs(title = "Midge Species Distribution Model",
        x = "Easting", y = "Northing")

# Create final risk map
risk_map <- create_risk_map(animal_ud, midge_prediction)
plot(risk_map)

# Plot final risk map
ggplot() +
   geom_sf(data = study_area) +
   tidyterra::geom_spatraster(data = risk_map) +
   geom_sf(data = water_bodies, alpha = 0.5) +
   geom_sf(data = feeders, color = "red") +
   scale_fill_distiller(palette = "BrBG") +
   # geom_sf(data = feeders, color = "white") +
   # geom_sf(data = water_bodies, fill = "grey80", alpha = 0.5) +
   # scale_fill_viridis_c(name = "Disease Risk", option = "magma") +
   theme_void() +
   labs(title = "Normalized Disease Risk Map",
        subtitle = "Animal Use × Midge Probability",
        x = "Easting", y = "Northing")


# Calculate the risk for each feeder
# - risk_map: output from create_risk_map()
# - feeders_sf: an sf object with feeder locations

# Calculate basic point risk
feeder_risks <- calculate_feeder_risk(
   risk_surface = risk_map,
   feeders = feeders
)
feeder_risks |> arrange(desc(risk_point))

# Calculate more comprehensive metrics with a 100m buffer
feeder_risks_detailed <- calculate_feeder_risk(
   risk_surface = risk_map,
   feeders = feeders,
   buffer_radius = 30,
   metrics = c("point", "mean", "max", "quantile"),
   quantile_probs = c(0.5, 0.75, 0.9, 0.95)
)

feeder_risks_detailed
