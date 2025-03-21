# First steps in simulation

# javirudolph 28-Mar-2025

# I would like this to be a complete simulation
# so, no external packages needed

# Libraries ---------------------------------------------------------------
library(tidyverse)
library(sf)
library(terra)

source("R/sim_landscape_fxs.R")

# Create spatial data -----------------------------------------------------


# Main simulation
# Create study area
boundary <- create_study_area(width = 5000, height = 5000)
plot(boundary)

# Create water bodies
water <- create_water_bodies(boundary, n_lakes = 5)
plot(water)

# Create feeders
feeders <- create_feeders(boundary, n_feeders = 8)
plot(feeders)
st_buffer(feeders, 50) |> ggplot() + geom_sf() + theme_minimal()

r <- rast(ext = st_bbox(boundary), resolution = 10)  # Adjust resolution as needed
r <- rasterize(boundary, r, field = 1)
feeder_raster <- rasterize(feeders, r, field = 1)

# Distance to feeders (for attraction)
feeder_dist <- distance(feeder_raster)

# Create environmental rasters
env_rasters <- create_env_rasters(boundary, resolution = 50)
temperature <- env_rasters[[1]]
humidity <- env_rasters[[2]]
vegetation <- env_rasters[[3]]

# Create output directory if it doesn't exist
# dir.create("simulated_data", showWarnings = FALSE)

# # Save shapefiles using sf
# st_write(boundary, "simulated_data/boundary.shp", delete_layer = TRUE)
# st_write(water_bodies, "simulated_data/water.shp", delete_layer = TRUE)
# st_write(feeders, "simulated_data/feeders.shp", delete_layer = TRUE)
#
# # Save rasters using terra
# writeRaster(env_rasters$temperature, "simulated_data/temperature.tif", overwrite = TRUE)
# writeRaster(env_rasters$humidity, "simulated_data/humidity.tif", overwrite = TRUE)
# writeRaster(env_rasters$vegetation, "simulated_data/vegetation.tif", overwrite = TRUE)

# Display a quick visualization
par(mfrow = c(2, 2))
plot(boundary$geometry, main = "Study Area with Features")
plot(water_bodies$geometry, col = "blue", add = TRUE)
plot(feeders$geometry, col = "red", pch = 19, add = TRUE)

# Plot rasters with terra
plot(env_rasters$temperature, main = "Temperature")
plot(env_rasters$humidity, main = "Humidity")
plot(env_rasters$vegetation, main = "Vegetation")



# Animal movement ---------------------------------------------------------

# Create movement simulation function
simulate_animal_movement <- function(start_point, n_steps,
                                     step_length_mean, step_length_sd,
                                     turning_angle_kappa) {

   # Initialize trajectory
   trajectory <- matrix(NA, nrow = n_steps + 1, ncol = 2)
   trajectory[1, ] <- start_point

   # Current direction (random start)
   current_direction <- runif(1, 0, 2*pi)

   for(i in 1:n_steps) {
      # Generate step length from gamma distribution
      step <- rgamma(1, shape = step_length_mean^2/step_length_sd^2,
                     scale = step_length_sd^2/step_length_mean)

      # Calculate potential new positions by sampling multiple angles
      trial_angles <- seq(0, 2*pi, length.out = 36)
      valid_positions <- matrix(NA, length(trial_angles), 2)
      attraction_scores <- numeric(length(trial_angles))

      for(j in 1:length(trial_angles)) {
         # Calculate potential position
         new_direction <- rwrappedcauchy(1, current_direction, turning_angle_kappa)

         # Add bias towards feeders
         feeder_bias <- trial_angles[j]

         # Combine direction with bias
         trial_direction <- circular::mean.circular(c(new_direction, feeder_bias),
                                                    c(0.7, 0.3))  # Weights

         pot_x <- trajectory[i, 1] + step * cos(trial_direction)
         pot_y <- trajectory[i, 2] + step * sin(trial_direction)

         valid_positions[j, ] <- c(pot_x, pot_y)

         # Check if position is within boundary and not in water
         test_point <- c(pot_x, pot_y)
         cell_vals <- extract(c(r, water_raster), matrix(test_point, ncol=2))

         if(!is.na(cell_vals[1,1])) {
            if(cell_vals[1,1] == 1 && (is.na(cell_vals[1,2]) || cell_vals[1,2] == 0)) {
               # Calculate attraction score based on distance to feeders
               attraction_score <- extract(feeder_dist, matrix(test_point, ncol=2))
               attraction_scores[j] <- attraction_score[1,1]
            } else {
               attraction_scores[j] <- NA
            }
         } else {
            attraction_scores[j] <- NA
         }
      }

      # Select position with best score (lower distance to feeders is better)
      valid_idx <- which(!is.na(attraction_scores))

      if(length(valid_idx) > 0) {
         best_idx <- valid_idx[which.min(attraction_scores[valid_idx])]
         trajectory[i+1, ] <- valid_positions[best_idx, ]

         # Update current direction
         current_direction <- atan2(trajectory[i+1, 2] - trajectory[i, 2],
                                    trajectory[i+1, 1] - trajectory[i, 1])
      } else {
         # If no valid position, stay in place
         trajectory[i+1, ] <- trajectory[i, ]
      }
   }

   return(trajectory)
}


# Disease vector ----------------------------------------------------------


