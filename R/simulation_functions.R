# Simulation of animal movement and midge distribution for disease risk assessment
# Using sf and terra packages as requested

# Load required packages
# library(sf)
# library(terra)
# library(tidyterra)
# library(ggplot2)
# library(dplyr)

# set.seed(123) # For reproducibility

# 1. Create simulated landscape -----------------

# Create a simulated boundary (study area)
create_study_area <- function(xmin = 0, xmax = 1000, ymin = 0, ymax = 1000, crs_code = 32616) {
   # Create polygon for the study area
   study_area <- st_polygon(list(rbind(
      c(xmin, ymin),
      c(xmax, ymin),
      c(xmax, ymax),
      c(xmin, ymax),
      c(xmin, ymin)
   )))

   # Convert to sf object
   study_area_sf <- st_sfc(study_area, crs = crs_code) # UTM 16N as an example
   study_area_sf <- st_sf(geometry = study_area_sf)

   return(study_area_sf)
}


# Create simulated water bodies (several small lakes)
create_water_bodies <- function(boundary_sf, n_lakes = 2, max_attempts = 100) {
   # Get the extent of the boundary
   ext <- st_bbox(boundary_sf)
   water_polys <- list()

   for(i in 1:n_lakes) {
      attempts <- 0
      valid_lake <- FALSE

      while(!valid_lake && attempts < max_attempts) {
         attempts <- attempts + 1

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
         candidate_poly <- st_polygon(list(lake_coords))

         # Check if the new lake overlaps with any existing lakes
         if(length(water_polys) == 0) {
            valid_lake <- TRUE  # First lake is always valid
         } else {
            # Create temporary sf object for existing lakes
            existing_lakes <- st_sfc(water_polys)
            candidate_lake <- st_sfc(candidate_poly)

            # Set CRS to match boundary
            st_crs(existing_lakes) <- st_crs(boundary_sf)
            st_crs(candidate_lake) <- st_crs(boundary_sf)

            # Check for intersection
            intersects <- any(st_intersects(existing_lakes, candidate_lake, sparse = FALSE))

            if(!intersects) {
               valid_lake <- TRUE
            }
         }

         if(valid_lake) {
            water_polys[[i]] <- candidate_poly
         }
      }

      # If we couldn't place a non-overlapping lake after max attempts, warn the user
      if(!valid_lake) {
         warning(paste("Could not place lake", i, "without overlap after", max_attempts, "attempts"))
         # Create a very small lake as a placeholder
         tiny_radius <- (ext["xmax"] - ext["xmin"]) * 0.01
         angles <- seq(0, 2*pi, length.out = 20)
         circle_x <- ext["xmin"] + (ext["xmax"] - ext["xmin"]) * 0.1 * i
         circle_y <- ext["ymin"] + (ext["ymax"] - ext["ymin"]) * 0.1
         lake_coords <- cbind(circle_x + tiny_radius * cos(angles),
                              circle_y + tiny_radius * sin(angles))
         lake_coords <- rbind(lake_coords, lake_coords[1,])  # Close the polygon
         water_polys[[i]] <- st_polygon(list(lake_coords))
      }
   }

   # Create sf object
   water_sf <- st_sf(geometry = st_sfc(water_polys), id = 1:length(water_polys))
   st_crs(water_sf) <- st_crs(boundary_sf)

   return(water_sf)
}

# Create feeder locations that attract animals
# New feeders function that doesn't overlap water
create_feeders <- function(study_area, n = 5, water_bodies = NULL, buffer_dist = 50, max_attempts = 100) {
   # Extract bounding box
   ext <- st_bbox(study_area)

   # Set a minimum distance between feeders (using buffer_dist)
   feeders_list <- list()
   placed_count <- 0

   # If water_bodies is provided, ensure it has the same CRS as study_area
   if(!is.null(water_bodies)) {
      if(st_crs(water_bodies) != st_crs(study_area)) {
         water_bodies <- st_transform(water_bodies, st_crs(study_area))
      }
      # Create a unified water polygon for faster intersection checks
      water_union <- st_union(water_bodies)
   }

   # Try to place each feeder
   for(i in 1:n) {
      attempts <- 0
      valid_feeder <- FALSE

      while(!valid_feeder && attempts < max_attempts) {
         attempts <- attempts + 1

         # Generate random location
         # x <- runif(1, bbox["xmin"] + 100, bbox["xmax"] - 100)
         # y <- runif(1, bbox["ymin"] + 100, bbox["ymax"] - 100)
         x <- runif(1, ext["xmin"] + (ext["xmax"] - ext["xmin"]) * 0.1,
                           ext["xmax"] - (ext["xmax"] - ext["xmin"]) * 0.1)
         y <- runif(1, ext["ymin"] + (ext["ymax"] - ext["ymin"]) * 0.1,
                           ext["ymax"] - (ext["ymax"] - ext["ymin"]) * 0.1)

         # Create candidate point
         candidate_point <- st_point(c(x, y))
         candidate_sfc <- st_sfc(candidate_point, crs = st_crs(study_area))

         # Create buffer around point for checking distance to other feeders
         candidate_buffer <- st_buffer(candidate_sfc, buffer_dist)

         # Check if point is valid:
         # 1. Doesn't overlap with existing feeders
         # 2. Not in a water body (if specified)

         # First check: Is it in the study area?
         in_study_area <- st_intersects(study_area, candidate_sfc, sparse = FALSE)[1]

         if(in_study_area) {
            # Second check: Is it away from other feeders?
            if(length(feeders_list) == 0) {
               away_from_feeders <- TRUE  # First feeder is always valid
            } else {
               # Create temporary sf object for existing feeders with buffer
               existing_buffers <- st_sfc(lapply(feeders_list, function(p) st_buffer(p, buffer_dist)),
                                          crs = st_crs(study_area))

               # Check if new feeder's buffer intersects with any existing feeder's buffer
               intersects_feeders <- any(st_intersects(existing_buffers, candidate_buffer, sparse = FALSE))
               away_from_feeders <- !intersects_feeders
            }

            # Third check (optional): Is it outside water bodies?
            if(!is.null(water_bodies)) {
               outside_water <- !any(st_intersects(water_union, candidate_sfc, sparse = FALSE))
            } else {
               outside_water <- TRUE  # No water bodies to check
            }

            # If all conditions are met, accept this feeder
            if(away_from_feeders && outside_water) {
               valid_feeder <- TRUE
            }
         }

         if(valid_feeder) {
            feeders_list[[i]] <- candidate_point
            placed_count <- placed_count + 1
         }
      }

      # If we couldn't place this feeder after max attempts, warn the user
      if(!valid_feeder) {
         warning(paste("Could not place feeder", i, "with valid constraints after", max_attempts, "attempts"))
      }
   }

   if(placed_count == 0) {
      stop("Could not place any feeders with the given constraints")
   }

   # Create sf object from the successfully placed feeders
   feeders_sfc <- st_sfc(feeders_list, crs = st_crs(study_area))
   feeders_sf <- st_sf(id = 1:length(feeders_list), geometry = feeders_sfc)

   return(feeders_sf)
}

# Generate environmental rasters for midge habitat modeling
create_env_rasters <- function(study_area, water_bodies) {
   # Convert study area to a SpatRaster for extent
   r <- rast(ext(st_bbox(study_area)), resolution = 10, crs = st_crs(study_area)$wkt)

   # Create environmental layers
   env_layers <- list()

   # Layer 1: Elevation (higher in some areas, lower in others)
   elevation <- rast(r)
   m <- matrix(ncol = ncol(elevation), nrow = nrow(elevation))
   for (i in 1:nrow(m)) {
      for (j in 1:ncol(m)) {
         # Create a surface with some hills and valleys
         x <- j / ncol(m) * 10
         y <- i / nrow(m) * 10
         m[i, j] <- 100 + 50 * sin(x) * cos(y) + 20 * sin(3 * x) + 10 * cos(5 * y)
      }
   }
   values(elevation) <- as.vector(m)
   names(elevation) <- "elevation"
   env_layers$elevation <- elevation

   # Layer 2: Distance to water (midges often breed near water)
   water_dist <- rast(r)
   water_dist_vals <- distance(rasterize(vect(water_bodies), r))
   values(water_dist) <- values(water_dist_vals)
   names(water_dist) <- "water_dist"
   env_layers$water_dist <- water_dist

   # Layer 3: Vegetation index (random but with spatial autocorrelation)
   veg_index <- rast(r)
   nr <- nrow(veg_index)
   nc <- ncol(veg_index)
   m <- matrix(ncol = nc, nrow = nr)

   # Create a surface with spatial autocorrelation
   base_vals <- matrix(rnorm(25), nrow = 5)
   # Expand to full size using bilinear interpolation
   for (i in 1:nr) {
      for (j in 1:nc) {
         # Map to positions in base_vals
         bi <- 1 + (i - 1) / (nr - 1) * (nrow(base_vals) - 1)
         bj <- 1 + (j - 1) / (nc - 1) * (ncol(base_vals) - 1)

         # Interpolate
         i1 <- floor(bi)
         i2 <- ceiling(bi)
         j1 <- floor(bj)
         j2 <- ceiling(bj)

         # Handle edge cases
         if (i2 > nrow(base_vals)) i2 <- nrow(base_vals)
         if (j2 > ncol(base_vals)) j2 <- ncol(base_vals)

         # Calculate weights
         wi <- bi - i1
         wj <- bj - j1

         # Bilinear interpolation
         if (i1 == i2 && j1 == j2) {
            m[i, j] <- base_vals[i1, j1]
         } else if (i1 == i2) {
            m[i, j] <- (1 - wj) * base_vals[i1, j1] + wj * base_vals[i1, j2]
         } else if (j1 == j2) {
            m[i, j] <- (1 - wi) * base_vals[i1, j1] + wi * base_vals[i2, j1]
         } else {
            m[i, j] <- (1 - wi) * (1 - wj) * base_vals[i1, j1] +
               wi * (1 - wj) * base_vals[i2, j1] +
               (1 - wi) * wj * base_vals[i1, j2] +
               wi * wj * base_vals[i2, j2]
         }
      }
   }

   # Scale to 0-1 range
   m <- (m - min(m)) / (max(m) - min(m))
   values(veg_index) <- as.vector(m)
   names(veg_index) <- "veg_index"
   env_layers$veg_index <- veg_index


   # Layer 4:Temperature raster (gradient from south to north with some noise)
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
   env_layers$temperature <- temperature

   # Combine all layers into a SpatRaster
   env_stack <- rast(env_layers)

   return(env_stack)
}

# 2. Animal movement simulation -----------------

# Simulate movement using a correlated random walk
simulate_animal_movement <- function(
      study_area,
      water_bodies,
      feeders,
      n_animals = 5,
      n_steps = 1000,
      step_scale = 20,
      k = 3,  # concentration parameter for von Mises distribution
      feeder_attraction = 0.2,  # strength of attraction to feeders
      max_feeder_influence_dist = 500  # max distance at which feeders influence movement
) {

   # Extract study area bounding box
   bbox <- st_bbox(study_area)

   # Initialize animal positions (random within study area but not in water)
   init_positions <- list()

   # Calculate buffer as a percentage of study area dimensions
   # Using 5% of the smallest dimension as buffer
   width <- bbox["xmax"] - bbox["xmin"]
   height <- bbox["ymax"] - bbox["ymin"]
   buffer_size <- 0.05 * min(width, height)

   for (i in 1:n_animals) {
      valid_position <- FALSE
      while (!valid_position) {
         x <- runif(1, bbox["xmin"] + buffer_size, bbox["xmax"] - buffer_size)
         y <- runif(1, bbox["ymin"] + buffer_size, bbox["ymax"] - buffer_size)
         pt <- st_point(c(x, y))
         pt_sf <- st_sf(geometry = st_sfc(pt, crs = st_crs(study_area)))

         # Check if inside study area and not in water
         if (st_within(pt_sf, study_area, sparse = FALSE)[1] &&
             !any(st_intersects(pt_sf, water_bodies, sparse = FALSE))) {
            init_positions[[i]] <- c(x, y)
            valid_position <- TRUE
         }
      }
   }

   # Function to calculate angle towards nearest feeder
   angle_to_nearest_feeder <- function(pos, feeders) {
      feeder_coords <- st_coordinates(feeders)[, 1:2]
      dists <- apply(feeder_coords, 1, function(f_pos) {
         sqrt(sum((pos - f_pos)^2))
      })
      nearest_idx <- which.min(dists)
      nearest_pos <- feeder_coords[nearest_idx, ]
      angle <- atan2(nearest_pos[2] - pos[2], nearest_pos[1] - pos[1])
      min_dist <- dists[nearest_idx]

      return(list(angle = angle, distance = min_dist))
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

         # Calculate angle towards nearest feeder and distance
         feeder_info <- angle_to_nearest_feeder(current_pos, feeders)
         feeder_angle <- feeder_info$angle
         feeder_dist <- feeder_info$distance

         # Adjust feeder attraction based on distance
         # Stronger attraction when closer to feeder, weaker when far away
         distance_factor <- max(0, 1 - feeder_dist / max_feeder_influence_dist)
         effective_attraction <- feeder_attraction * distance_factor

         # Determine if this step will be influenced by feeders
         feeder_influenced <- runif(1) < effective_attraction

         if (feeder_influenced) {
            # Move towards feeder with some random variation
            turning_angle <- rnorm(1, 0, pi/4)  # More variability
            proposed_angle <- feeder_angle + turning_angle
         } else {
            # Correlated random walk - direction correlated with previous
            turning_angle <- rnorm(1, 0, pi/k)
            proposed_angle <- current_angle + turning_angle
         }

         # Generate step length from gamma distribution
         # Shorter steps when near feeder, longer steps when away
         local_scale <- step_scale * (1 - 0.5 * distance_factor)  # Adjust step length based on feeder proximity
         step_length <- rgamma(1, shape = 2, scale = local_scale)

         # Calculate proposed new position
         proposed_x <- current_pos[1] + step_length * cos(proposed_angle)
         proposed_y <- current_pos[2] + step_length * sin(proposed_angle)
         proposed_pos <- c(proposed_x, proposed_y)

         # Check if proposed position is valid (within study area and not in water)
         pt <- st_point(proposed_pos)
         pt_sf <- st_sf(geometry = st_sfc(pt, crs = st_crs(study_area)))

         if (st_within(pt_sf, study_area, sparse = FALSE)[1] &&
             !any(st_intersects(pt_sf, water_bodies, sparse = FALSE))) {
            # Accept move
            positions[t + 1, ] <- proposed_pos
            current_angle <- proposed_angle
         } else {
            # Water avoidance behavior
            # Try multiple alternative angles to find a valid move
            found_valid_move <- FALSE

            for (attempt in 1:8) {  # Try 8 different directions
               # Turn away from proposed direction
               avoidance_angle <- proposed_angle + attempt * pi/4

               # Shorter step when avoiding obstacles
               avoidance_step <- step_length * 0.6

               # Calculate new proposed position
               avoid_x <- current_pos[1] + avoidance_step * cos(avoidance_angle)
               avoid_y <- current_pos[2] + avoidance_step * sin(avoidance_angle)
               avoid_pos <- c(avoid_x, avoid_y)

               # Check if valid
               avoid_pt <- st_point(avoid_pos)
               avoid_pt_sf <- st_sf(geometry = st_sfc(avoid_pt, crs = st_crs(study_area)))

               if (st_within(avoid_pt_sf, study_area, sparse = FALSE)[1] &&
                   !any(st_intersects(avoid_pt_sf, water_bodies, sparse = FALSE))) {
                  # Found valid move
                  positions[t + 1, ] <- avoid_pos
                  current_angle <- avoidance_angle
                  found_valid_move <- TRUE
                  break
               }
            }

            if (!found_valid_move) {
               # If no valid move found, stay in place
               positions[t + 1, ] <- current_pos
               # Change angle more drastically for next attempt
               current_angle <- runif(1, 0, 2 * pi)
            }
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

# Create kernel density estimation of animal space use

create_animal_ud <- function(
      animal_tracks,
      study_area,
      resolution = 10,  # Default resolution in map units (e.g., meters)
      smoothing_factor = 5  # Size of smoothing window
) {
   # Extract points for KDE
   points <- st_coordinates(animal_tracks)

   # Create a template raster based on the study area
   # Get the extent from the study area
   sa_bbox <- st_bbox(study_area)

   # Create template raster with specified resolution
   template_rast <- rast(ext(st_bbox(study_area)), resolution = resolution, crs = st_crs(study_area)$wkt)

   # Initialize raster with zeros
   kde_rast <- template_rast
   values(kde_rast) <- 0

   # Mask the raster to the study area (optional but recommended)
   study_area_vect <- vect(study_area)
   kde_rast <- mask(kde_rast, study_area_vect)

   # Convert animal track points to vector for rasterization
   animal_vect <- vect(points)

   # Count points in each cell
   utilization <- rasterize(animal_vect, kde_rast, fun = "count", background = 0)

   # Add small value to avoid zeros (for log transformations if needed later)
   utilization <- utilization + 0.001

   # Smooth the density using a focal operation (moving window)
   # Create window matrix based on smoothing factor
   w <- matrix(1, smoothing_factor, smoothing_factor)
   smoothed_kde <- focal(utilization, w = w, fun = mean, na.rm = TRUE)

   # Normalize to 0-1 scale
   ud_min <- minmax(smoothed_kde)[1]
   ud_max <- minmax(smoothed_kde)[2]

   # Check if min equals max (all values the same)
   if (ud_min == ud_max) {
      ud_rast <- smoothed_kde
      ud_rast[] <- 1  # Set all to 1 if no variation
   } else {
      ud_rast <- (smoothed_kde - ud_min) / (ud_max - ud_min)
   }

   # Add small value to avoid zeros
   ud_rast <- ud_rast + 0.0001

   # Set name
   names(ud_rast) <- "animal_ud"

   return(ud_rast)
}


# 3. Midge distribution simulation -----------------

# Simulate midge presence/absence data based on environmental variables
simulate_midge_data <- function(env_rasters, water_bodies, n_samples = 500) {
   # Create random sampling points throughout the study area
   r_template <- env_rasters[[1]]

   # Generate more points than needed because some will be filtered out
   n_buffer <- round(n_samples * 1.5)  # Generate 50% more points as buffer
   random_cells <- sample(1:ncell(r_template), n_buffer)
   xy <- xyFromCell(r_template, random_cells)

   # Convert the xy matrix to a data frame
   xy_df <- as.data.frame(xy)
   names(xy_df) <- c("x", "y")

   # Convert to sf object
   points_sf <- st_as_sf(xy_df,
                    coords = c("x", "y"),
                    crs = st_crs(crs(r_template)))

   # Filter out points that fall within water bodies
   if (st_crs(points_sf) != st_crs(water_bodies)) {
      points_sf <- st_transform(points_sf, st_crs(water_bodies))
   }

   intersection <- st_intersects(points_sf, water_bodies)
   no_intersection <- lengths(intersection) == 0
   points_outside_water <- points_sf[no_intersection, ]

   # If we have too few points after filtering, generate more
   if (nrow(points_outside_water) < n_samples) {
      warning("Not enough points outside water bodies. Regenerating more points...")
      return(simulate_midge_data(env_rasters, water_bodies, n_samples))
   }

   # Take only the required number of points
   points_outside_water <- slice_sample(points_outside_water, n = n_samples)

   # Extract coordinates from sf object
   xy_filtered <- st_coordinates(points_outside_water)

   # Extract environmental values
   env_values <- extract(env_rasters, xy_filtered)

   # Define relationship between environmental variables and midge presence
   # Water proximity increases probability
   # Vegetation has moderate positive effect
   # Elevation has negative effect (midges prefer lower areas)
   # Normalize environmental variables to 0-1 scale for easier coefficient interpretation
   env_values$elevation <- scale01(env_values$elevation)
   env_values$water_dist <- scale01(env_values$water_dist)
   env_values$temperature <- scale01(env_values$temperature)

   # Inverse of water distance (proximity)
   env_values$water_prox <- 1 - env_values$water_dist

   # Calculate linear predictor
   linear_pred <- -1 +
      3 * env_values$water_prox + # Strong positive effect of water proximity
      0.5 * env_values$veg_index + # Moderate positive effect of vegetation
      -2 * env_values$elevation + # Negative effect of elevation
      -3 * env_values$temperature # Strong negative effect of temperature

   # Add some random noise
   linear_pred <- linear_pred + rnorm(n_samples, 0, 0.2)

   # Convert to probability using logistic function
   prob_presence <- 1 / (1 + exp(-linear_pred))

   # Generate presence/absence data
   presence <- rbinom(n_samples, 1, prob_presence)

   # Create sf object with the data
   midge_data <- st_sf(
      presence = presence,
      elevation = env_values$elevation,
      water_dist = env_values$water_dist,
      veg_index = env_values$veg_index,
      temperature = env_values$temperature,
      geometry = points_outside_water$geometry
   )

   return(list(midge_data = midge_data, true_prob = prob_presence))
}

# Fit a species distribution model for midges
fit_midge_sdm <- function(midge_data, env_rasters) {
   # Fit a logistic regression model
   sdm_model <- glm(
      presence ~ elevation + water_dist + veg_index + temperature,
      family = binomial(link = "logit"),
      data = midge_data
   )

   # Print model summary
   print(summary(sdm_model))

   # Check if layer names match model variables
   print("Raster layer names:")
   print(names(env_rasters))
   print("Model variables:")
   print(names(coef(sdm_model)))

   env_rasters_scaled <- env_rasters
   for (var in c("elevation", "water_dist", "veg_index", "temperature")) {
      if (var %in% names(env_rasters)) {
         # Scale each raster layer using the same parameters as the training data
         env_rasters_scaled[[var]] <- app(env_rasters[[var]], scale01)
      } else {
         warning(paste("Variable", var, "not found in raster stack"))
      }
   }


   # Predict to the entire study area
   sdm_prediction <- terra::predict(env_rasters_scaled,
                                    sdm_model, type = "response")
   names(sdm_prediction) <- "midge_prob"

   return(list(model = sdm_model, prediction = sdm_prediction))
}


# 4. Combine animal use and midge risk to create risk maps -----------------

create_risk_map <- function(animal_ud, midge_sdm_prediction, normalize = TRUE) {

   # Ensure rasters are aligned
   if (crs(animal_ud) != crs(midge_sdm_prediction) ||
       ext(animal_ud) != ext(midge_sdm_prediction) ||
       ncol(animal_ud) != ncol(midge_sdm_prediction) ||
       nrow(animal_ud) != nrow(midge_sdm_prediction)) {
      midge_sdm_prediction <- resample(midge_sdm_prediction, animal_ud)
   }

   # Calculate combined risk as the product of animal use and midge probability
   risk_map <- animal_ud * midge_sdm_prediction
   names(risk_map) <- "disease_risk"

   if(normalize) {
      # Normalize to 0-1 scale
      risk_map <- (risk_map - minmax(risk_map)[1]) / (minmax(risk_map)[2] - minmax(risk_map)[1]) + 0.001
      names(risk_map) <- "disease_risk"
   }

   return(risk_map)
}

# want to calculate specific risk for each feeder
calculate_feeder_risk <- function(
      risk_surface,       # Risk surface raster (output from create_risk_map)
      feeders,            # SF object containing feeder locations
      buffer_radius = 0,  # Optional buffer radius around feeders (in map units)
      metrics = c("point", "mean", "max", "quantile"), # Risk metrics to calculate
      quantile_probs = c(0.5, 0.75, 0.9)  # Quantiles to calculate if "quantile" in metrics
) {
   # Ensure feeders is an sf object
   if (!inherits(feeders, "sf")) {
      stop("feeders must be an sf object")
   }

   # Ensure CRS match
   if (sf::st_crs(feeders) != terra::crs(risk_surface)) {
      feeders <- sf::st_transform(feeders, terra::crs(risk_surface))
   }

   # Initialize results dataframe with feeder coordinates
   feeder_coords <- sf::st_coordinates(feeders)
   results <- data.frame(
      feeder_id = seq_len(nrow(feeders)),
      x = feeder_coords[, 1],
      y = feeder_coords[, 2]
   )

   # Add any attributes from the original feeders object
   if (ncol(sf::st_drop_geometry(feeders)) > 0) {
      results <- cbind(results, sf::st_drop_geometry(feeders))
   }

   # Extract point risk values at feeder locations
   if ("point" %in% metrics) {
      point_values <- terra::extract(risk_surface, feeder_coords)
      results$risk_point <- point_values$disease_risk
   }

   # If buffer radius > 0, calculate additional metrics within buffer
   if (buffer_radius > 0) {
      # Create buffered areas around feeders
      feeder_buffers <- sf::st_buffer(feeders, dist = buffer_radius)

      # Loop through each feeder to calculate metrics within its buffer
      for (i in 1:nrow(feeders)) {
         # Extract single buffer
         buffer <- feeder_buffers[i, ]

         # Crop and mask risk surface to buffer
         buffer_risk <- terra::crop(risk_surface, terra::ext(buffer))
         buffer_risk <- terra::mask(buffer_risk, terra::vect(buffer))

         # Extract values within buffer
         buffer_values <- terra::values(buffer_risk)
         buffer_values <- buffer_values[!is.na(buffer_values)]

         # Skip calculations if no valid values in buffer
         if (length(buffer_values) == 0) {
            if ("mean" %in% metrics) results$risk_mean[i] <- NA
            if ("max" %in% metrics) results$risk_max[i] <- NA
            if ("quantile" %in% metrics) {
               for (prob in quantile_probs) {
                  col_name <- paste0("risk_q", prob * 100)
                  results[[col_name]][i] <- NA
               }
            }
            next
         }

         # Calculate requested metrics
         if ("mean" %in% metrics) {
            results$risk_mean[i] <- mean(buffer_values, na.rm = TRUE)
         }

         if ("max" %in% metrics) {
            results$risk_max[i] <- max(buffer_values, na.rm = TRUE)
         }

         if ("quantile" %in% metrics) {
            quant_values <- stats::quantile(buffer_values, probs = quantile_probs, na.rm = TRUE)
            for (j in 1:length(quantile_probs)) {
               col_name <- paste0("risk_q", quantile_probs[j] * 100)
               results[[col_name]][i] <- quant_values[j]
            }
         }
      }
   }

   return(results)
}



## Accessory functions
# Helper function to scale values to 0-1
scale01 <- function(x) {
   return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

