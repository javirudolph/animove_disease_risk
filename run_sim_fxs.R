# Source the existing simulation functions
source("simulation_functions.R")

# 5. Run the simulation -----------------

# Create the simulated landscape
set.seed(508)
study_area <- create_study_area()
plot(study_area)
water_bodies <- create_water_bodies(study_area, n_lakes = 3)
# plot(water_bodies)
ggplot() +
   geom_sf(data = study_area) +
   geom_sf(data = water_bodies)
feeders <- create_feeders(study_area)
plot(feeders)
env_rasters <- create_env_rasters(study_area)
plot(env_rasters)

# Simulate animal movement
animal_tracks <- simulate_animal_movement(
   study_area = study_area,
   water_bodies = water_bodies,
   feeders = feeders,
   n_animals = 3,
   n_steps = 1000
)

animal_tracks |>
   # filter(animal_id == 1) |>
   ggplot() +
   geom_sf(aes(color = factor(animal_id)))


# Create animal utilization distribution
animal_ud <- create_animal_ud(animal_tracks, env_rasters)
plot(animal_ud)

# Simulate midge data
midge_sim <- simulate_midge_data(env_rasters)
head(midge_sim)

midge_data <- midge_sim$midge_data
midge_data |>
   ggplot() +
   geom_sf(aes(color = factor(presence)))

# Fit midge distribution model
midge_sdm <- fit_midge_sdm(midge_data, env_rasters)
midge_prediction <- midge_sdm$prediction
plot(midge_prediction)

# Create final risk map
risk_map <- create_risk_map(animal_ud, midge_prediction)
plot(risk_map)


# 6. Plot results -----------------

# Convert rasters to data frames for plotting
animal_ud_df <- as.data.frame(animal_ud, xy = TRUE)
midge_prob_df <- as.data.frame(midge_prediction, xy = TRUE)
risk_map_df <- as.data.frame(risk_map, xy = TRUE)

# Convert sf objects to data frames for plotting
study_area_df <- st_coordinates(st_cast(study_area, "MULTILINESTRING"))
water_bodies_df <- st_coordinates(st_cast(water_bodies, "MULTILINESTRING"))
feeders_df <- st_coordinates(feeders)
animal_points_df <- st_coordinates(animal_tracks)

# Plot animal movement and utilization distribution
animal_tracks |>
   mutate(x = st_coordinates(geometry)[,1],
          y = st_coordinates(geometry)[,2]) |>
   ggplot() +
   facet_wrap(~ factor(animal_id)) +
   geom_path(aes(x = x, y = y, group = factor(animal_id), color = factor(animal_id))) +
   theme_minimal() +
   geom_sf(data = feeders) +
   geom_sf(data = water_bodies)

# Plot study area
ggplot() +
   geom_sf(data = study_area, fill = "NA", color = "black", linewidth = 1) +
   geom_sf(data = feeders, aes(color = "Feeders")) +
   geom_sf(data = water_bodies, aes(fill = "Water\nBodies"), alpha = 0.5) +
   scale_color_manual(values = "red") +
   scale_fill_manual(values = "grey80") +
   theme_void() +
   labs(title = "Study Area",
        x = "Easting", y = "Northing") +
   theme(legend.title = element_blank()) -> plot_study_area

# Plot utilization distribution
ggplot() +
   geom_sf(data = study_area) +
   tidyterra::geom_spatraster(data = animal_ud) +
   scale_fill_viridis_c(name = "Utlization\nDistribution", option = "plasma") +
   # geom_sf(data = feeders, color = "white") +
   # geom_sf(data = water_bodies, fill = "grey80", alpha = 0.5) +
   theme_void() +
   labs(title = "Animal Utilization Distribution",
        x = "Easting", y = "Northing") -> plot_animal_ud

# Plot midge distribution model
ggplot() +
   geom_sf(data = study_area) +
   tidyterra::geom_spatraster(data = midge_prediction) +
   scale_fill_distiller(palette = "BrBG") +
   # geom_sf(data = feeders, color = "white") +
   # geom_sf(data = water_bodies, fill = "grey80", alpha = 0.5) +
   theme_minimal() +
   # labs(title = "Midge probability of occurrence")
   scale_fill_viridis_c(name = "Midge\nProbability", option = "plasma") +
   scale_shape_manual(values = c(4, 16), name = "Midge Presence") +
   theme_void() +
   labs(title = "Midge Species Distribution Model",
        x = "Easting", y = "Northing") -> plot_midge_prob

# Plot final risk map
ggplot() +
   geom_sf(data = study_area) +
   tidyterra::geom_spatraster(data = risk_map) +
   scale_fill_distiller(palette = "BrBG") +
   # geom_sf(data = feeders, color = "white") +
   # geom_sf(data = water_bodies, fill = "grey80", alpha = 0.5) +
   scale_fill_viridis_c(name = "Disease Risk", option = "magma") +
   theme_void() +
   labs(title = "Normalized Disease Risk Map",
        subtitle = "Animal Use × Midge Probability",
        x = "Easting", y = "Northing") -> plot_disease_risk


cowplot::plot_grid(
   plot_study_area,
   plot_animal_ud,
   plot_midge_prob,
   plot_disease_risk,
   ncol = 2
)
