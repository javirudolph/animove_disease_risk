# Create a simulated study area boundary (polygon)
create_study_area <- function(center_x = 0, center_y = 0, width = 1000, height = 1000) {

   # center_x = 0; center_y = 0; width = 1000; height = 1000
   # Create a rectangular boundary
   boundary_coords <- rbind(
      c(center_x - width/2, center_y - height/2),
      c(center_x + width/2, center_y - height/2),
      c(center_x + width/2, center_y + height/2),
      c(center_x - width/2, center_y + height/2),
      c(center_x - width/2, center_y - height/2)  # Close the polygon
   )

   # Create sf polygon
   boundary_sf <- st_polygon(list(boundary_coords)) |>
      st_sfc() |>
      st_sf(geometry = _, id = 1)

   # Add projection (UTM zone 17N as an example)
   st_crs(boundary_sf) <- "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"

   return(boundary_sf)
}

# Create simulated water bodies (several small lakes)
create_water_bodies <- function(boundary_sf, n_lakes = 3) {
   # Get the extent of the boundary
   ext <- st_bbox(boundary_sf)
   water_polys <- list()

   for(i in 1:n_lakes) {
      # Random center within boundary
      center_x <- runif(1, ext["xmin"] + (ext["xmax"] - ext["xmin"]) * 0.2,
                        ext["xmax"] - (ext["xmax"] - ext["xmin"]) * 0.2)
      center_y <- runif(1, ext["ymin"] + (ext["ymax"] - ext["ymin"]) * 0.2,
                        ext["ymax"] - (ext["ymax"] - ext["ymin"]) * 0.2)

      # Random radius (smaller lakes)
      radius <- runif(1, (ext["xmax"] - ext["xmin"]) * 0.05,
                      (ext["xmax"] - ext["xmin"]) * 0.1)

      # Create a circle approximation with 20 points
      angles <- seq(0, 2*pi, length.out = 20)
      circle_x <- center_x + radius * cos(angles)
      circle_y <- center_y + radius * sin(angles)
      lake_coords <- cbind(circle_x, circle_y)

      # Ensure the polygon is properly closed - first and last points must be identical
      if(!identical(lake_coords[1,], lake_coords[nrow(lake_coords),])) {
         lake_coords <- rbind(lake_coords, lake_coords[1,])
      }

      # Create polygon
      water_polys[[i]] <- st_polygon(list(lake_coords))
   }

   # Create sf object
   water_sf <- st_sf(geometry = st_sfc(water_polys), id = 1:n_lakes)
   st_crs(water_sf) <- st_crs(boundary_sf)

   return(water_sf)
}

# Create feeder locations (points)
create_feeders <- function(boundary_sf, n_feeders = 5) {
   # Get the extent of the boundary
   ext <- st_bbox(boundary_sf)

   # Generate random points
   x_coords <- runif(n_feeders, ext["xmin"] + (ext["xmax"] - ext["xmin"]) * 0.1,
                     ext["xmax"] - (ext["xmax"] - ext["xmin"]) * 0.1)
   y_coords <- runif(n_feeders, ext["ymin"] + (ext["ymax"] - ext["ymin"]) * 0.1,
                     ext["ymax"] - (ext["ymax"] - ext["ymin"]) * 0.1)

   # Create sf object
   feeders_sf <- st_as_sf(data.frame(x = x_coords, y = y_coords, id = 1:n_feeders),
                          coords = c("x", "y"))
   st_crs(feeders_sf) <- st_crs(boundary_sf)

   return(feeders_sf)
}

# Create simulated environmental rasters
# Create simulated environmental rasters using terra's focal operations
create_env_rasters <- function(boundary_sf, resolution = 10) {
   # Create a template raster with terra
   r <- rast(ext = st_bbox(boundary_sf), resolution = resolution)
   crs(r) <- st_crs(boundary_sf)$wkt

   # Temperature raster (gradient from south to north with some noise)
   temperature <- r
   # Create coordinates matrix
   coords <- crds(temperature)

   # Normalize y-coordinates
   y_coords <- coords[, 2]
   y_range <- range(y_coords)
   y_normalized <- (y_coords - y_range[1]) / (y_range[2] - y_range[1])

   # Warmer in south, cooler in north with noise
   temp_values <- 30 - 15 * y_normalized + rnorm(length(y_normalized), 0, 2)

   # Fill temperature raster
   temperature[] <- temp_values

   # Humidity raster with spatial autocorrelation
   # Start with random noise
   humidity <- r
   humidity[] <- runif(ncell(humidity), 40, 90)  # Random values between 40-90%

   # Apply gaussian smoothing to create spatial autocorrelation
   # We use several passes of focal operation with different window sizes
   humidity <- focal(humidity, w=matrix(1/9, nrow=3, ncol=3), fun="mean")
   humidity <- focal(humidity, w=matrix(1/25, nrow=5, ncol=5), fun="mean")
   humidity <- focal(humidity, w=matrix(1/49, nrow=7, ncol=7), fun="mean")

   # Add some relationship to temperature (higher humidity in cooler areas)
   humidity <- humidity - 0.2 * scale(temperature)

   # Rescale humidity to 40-90% range
   humidity <- stretch(humidity, 40, 90)

   # Vegetation raster (more patchy spatial pattern)
   # Start with random noise
   vegetation <- r
   vegetation[] <- runif(ncell(vegetation), 0, 100)

   # Create patches using focal smoothing
   # First a smaller window
   w1 <- matrix(1/9, nrow=3, ncol=3)
   vegetation <- focal(vegetation, w=w1, fun="mean")

   # Then add some variability back
   rand_overlay <- r
   rand_overlay[] <- runif(ncell(vegetation), -20, 20)
   vegetation <- vegetation + rand_overlay

   # # Clamp values to prevent extreme outliers
   vegetation <- clamp(vegetation, 0, 100)
   #
   # Now smooth again with a different window size
   w2 <- matrix(1/25, nrow=5, ncol=5)
   vegetation <- focal(vegetation, w=w2, fun="mean")

   # Create more distinct edges between patches
   # The minmax function in your code needs to be replaced with proper min/max calls
   veg_min <- min(values(vegetation), na.rm=TRUE)
   veg_max <- max(values(vegetation), na.rm=TRUE)
   vegetation <- (vegetation - veg_min) / (veg_max - veg_min) * 100

   # Apply a threshold effect to create more distinct patches
   # Need to extract values first, then apply cut, then put back into the raster
   veg_values <- values(vegetation)
   threshold <- quantile(veg_values, probs=seq(0, 1, 0.2), na.rm=TRUE)
   vegetation_classes <- cut(veg_values, breaks=threshold)

   # Convert the factor to numeric for the raster
   vegetation_classified <- r
   vegetation_classified[] <- as.numeric(vegetation_classes)

   # Make sure all rasters have the same extent
   temperature <- crop(temperature, ext(r))
   humidity <- crop(humidity, ext(r))
   vegetation <- crop(vegetation_classified, ext(r))

   return(list(temperature = temperature, humidity = humidity, vegetation = vegetation))
}


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
         cell_x <- cellFromXY(r, c(pot_x, pot_y))[1]

         if(!is.na(cell_x)) {
            if(r[cell_x] == 1 && (is.na(water_raster[cell_x]) || water_raster[cell_x] == 0)) {
               # Calculate attraction score based on distance to feeders
               attraction_scores[j] <- extract(feeder_dist, c(pot_x, pot_y))
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
