create_water_bodies_old <- function(boundary_sf, n_lakes = 2) {
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


create_feeders_old <- function(study_area, n = 5) {
   # Extract bounding box
   bbox <- st_bbox(study_area)

   # Create n random points for feeders
   x <- runif(n, bbox["xmin"] + 100, bbox["xmax"] - 100)
   y <- runif(n, bbox["ymin"] + 100, bbox["ymax"] - 100)

   # Create points
   feeders <- st_sfc(lapply(1:n, function(i) st_point(c(x[i], y[i]))), crs = st_crs(study_area))
   feeders <- st_sf(id = 1:n, geometry = feeders)

   return(feeders)
}

simulate_animal_movement_old <- function(
      study_area,
      water_bodies,
      feeders,
      n_animals = 5,
      n_steps = 1000,
      step_scale = 20,
      k = 3,  # concentration parameter for von Mises distribution
      feeder_attraction = 0.2  # strength of attraction to feeders
) {

   # Extract study area bounding box
   bbox <- st_bbox(study_area)

   # Initialize animal positions (random within study area but not in water)
   init_positions <- list()
   for (i in 1:n_animals) {
      valid_position <- FALSE
      while (!valid_position) {
         x <- runif(1, bbox["xmin"] + 50, bbox["xmax"] - 50)
         y <- runif(1, bbox["ymin"] + 50, bbox["ymax"] - 50)
         pt <- st_point(c(x, y))
         pt_sf <- st_sf(geometry = st_sfc(pt, crs = st_crs(study_area)))

         # Check if inside study area and not in water
         if (st_within(pt_sf, study_area, sparse = FALSE)[1] &&
             !any(st_within(pt_sf, water_bodies, sparse = FALSE))) {
            init_positions[[i]] <- c(x, y)
            valid_position <- TRUE
         }
      }
   }

   # Function to calculate angle towards nearest feeder
   angle_to_nearest_feeder <- function(pos, feeders) {
      dists <- apply(st_coordinates(feeders)[, 1:2], 1, function(f_pos) {
         sqrt(sum((pos - f_pos)^2))
      })
      nearest_idx <- which.min(dists)
      nearest_pos <- st_coordinates(feeders)[nearest_idx, 1:2]
      return(atan2(nearest_pos[2] - pos[2], nearest_pos[1] - pos[1]))
   }

   # Simulate movement for each animal
   animal_tracks <- list()

   for (a in 1:n_animals) {
      # Initialize
      positions <- matrix(NA, nrow = n_steps + 1, ncol = 2)
      positions[1, ] <- init_positions[[a]]

      # Initial direction (random)
      current_angle <- runif(1, 0, 2 * pi)

      for (t in 1:n_steps) {
         current_pos <- positions[t, ]

         # Calculate angle towards nearest feeder
         feeder_angle <- angle_to_nearest_feeder(current_pos, feeders)

         # Generate turning angle from von Mises distribution (wrapped normal)
         # Mix random direction with attraction to feeders
         if (runif(1) < feeder_attraction) {
            # Move towards feeder with some random variation
            turning_angle <- rnorm(1, 0, pi/8)
            proposed_angle <- feeder_angle + turning_angle
         } else {
            # Correlated random walk - direction correlated with previous
            # k controls concentration - higher k means less variable turning angles
            turning_angle <- rnorm(1, 0, pi/k)
            proposed_angle <- current_angle + turning_angle
         }

         # Generate step length from gamma distribution
         step_length <- rgamma(1, shape = 2, scale = step_scale)

         # Calculate proposed new position
         proposed_x <- current_pos[1] + step_length * cos(proposed_angle)
         proposed_y <- current_pos[2] + step_length * sin(proposed_angle)
         proposed_pos <- c(proposed_x, proposed_y)

         # Check if proposed position is valid (within study area and not in water)
         pt <- st_point(proposed_pos)
         pt_sf <- st_sf(geometry = st_sfc(pt, crs = st_crs(study_area)))

         if (st_within(pt_sf, study_area, sparse = FALSE)[1] &&
             !any(st_within(pt_sf, water_bodies, sparse = FALSE))) {
            # Accept move
            positions[t + 1, ] <- proposed_pos
            current_angle <- proposed_angle
         } else {
            # Reject move, stay in place
            positions[t + 1, ] <- current_pos
            # Change angle more drastically for next attempt
            current_angle <- runif(1, 0, 2 * pi)
         }
      }

      # Convert to sf object
      track_sf <- st_sf(
         animal_id = a,
         step = 0:n_steps,
         geometry = st_sfc(
            lapply(1:(n_steps + 1), function(i) st_point(positions[i, ])),
            crs = st_crs(study_area)
         )
      )

      animal_tracks[[a]] <- track_sf
   }

   # Combine all tracks
   all_tracks <- do.call(rbind, animal_tracks)

   return(all_tracks)
}


create_animal_ud <- function(animal_tracks, env_rasters) {
   # Extract points for KDE
   points <- st_coordinates(animal_tracks)

   # Create a KDE raster based on the extent and resolution of environmental rasters
   r_template <- env_rasters[[1]]

   # Create a raster of counts/density
   kde_rast <- rast(r_template)
   values(kde_rast) <- 0

   # Count points in each cell
   animal_vect <- vect(points)
   utilization <- rasterize(animal_vect, kde_rast, fun = "count", background = 0)
   utilization <- utilization + 0.001
   # plot(utilization)

   # total_points <- nrow(points)
   # utilization_prob <- utilization / total_points
   # plot(utilization_prob)

   # Smooth the density using a focal operation (moving window)
   smoothed_kde <- focal(utilization, w = matrix(1, 5, 5), fun = mean, na.rm = TRUE)
   # plot(smoothed_kde)

   # Normalize to 0-1 scale
   ud_rast <- (smoothed_kde - minmax(smoothed_kde)[1]) / (minmax(smoothed_kde)[2] - minmax(smoothed_kde)[1]) + 0.001
   names(ud_rast) <- "animal_ud"

   return(ud_rast)
}
