# Simulation of animal movement and midge distribution for disease risk assessment
# Using sf and terra packages as requested

# Load required packages
library(sf)
library(terra)
library(tidyterra)
library(ggplot2)
library(dplyr)

set.seed(123) # For reproducibility

# 1. Create simulated landscape -----------------

# Create a simulated boundary (study area)
create_study_area <- function(xmin = 0, xmax = 1000, ymin = 0, ymax = 1000) {
   # Create polygon for the study area
   study_area <- st_polygon(list(rbind(
      c(xmin, ymin),
      c(xmax, ymin),
      c(xmax, ymax),
      c(xmin, ymax),
      c(xmin, ymin)
   )))

   # Convert to sf object
   study_area_sf <- st_sfc(study_area, crs = 32616) # UTM 16N as an example
   study_area_sf <- st_sf(geometry = study_area_sf)

   return(study_area_sf)
}

# Create water bodies to be avoided by animals
# create_water_bodies <- function(study_area) {
#    # Extract bounding box
#    bbox <- st_bbox(study_area)
#
#    # Create two lakes
#    lake1 <- st_polygon(list(rbind(
#       c(bbox["xmin"] + 200, bbox["ymin"] + 200),
#       c(bbox["xmin"] + 300, bbox["ymin"] + 200),
#       c(bbox["xmin"] + 300, bbox["ymin"] + 300),
#       c(bbox["xmin"] + 200, bbox["ymin"] + 300),
#       c(bbox["xmin"] + 200, bbox["ymin"] + 200)
#    )))
#
#    lake2 <- st_polygon(list(rbind(
#       c(bbox["xmin"] + 700, bbox["ymin"] + 600),
#       c(bbox["xmin"] + 800, bbox["ymin"] + 600),
#       c(bbox["xmin"] + 850, bbox["ymin"] + 700),
#       c(bbox["xmin"] + 750, bbox["ymin"] + 750),
#       c(bbox["xmin"] + 700, bbox["ymin"] + 650),
#       c(bbox["xmin"] + 700, bbox["ymin"] + 600)
#    )))
#
#    # Convert to sf object
#    water_bodies <- st_sfc(lake1, lake2, crs = st_crs(study_area))
#    water_bodies <- st_sf(id = c(1, 2), geometry = water_bodies)
#
#    return(water_bodies)
# }

# Create simulated water bodies (several small lakes)
create_water_bodies <- function(boundary_sf, n_lakes = 2) {
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

# Create feeder locations that attract animals
create_feeders <- function(study_area, n = 5) {
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

# Generate environmental rasters for midge habitat modeling
create_env_rasters <- function(study_area, n_layers = 3) {
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
      n_animals = 10,
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

# Create kernel density estimation of animal space use
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


# 3. Midge distribution simulation -----------------

# Simulate midge presence/absence data based on environmental variables
simulate_midge_data <- function(env_rasters, n_samples = 500) {
   # Create random sampling points throughout the study area
   r_template <- env_rasters[[1]]
   random_cells <- sample(1:ncell(r_template), n_samples)
   xy <- xyFromCell(r_template, random_cells)

   # Extract environmental values
   env_values <- extract(env_rasters, xy)

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
      geometry = st_sfc(
         lapply(1:nrow(xy), function(i) st_point(xy[i, ])),
         crs = crs(r_template)
      )
   )

   return(list(midge_data = midge_data, true_prob = prob_presence))
}

# Helper function to scale values to 0-1
scale01 <- function(x) {
   return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
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

   midge_sdm_prediction <- midge_sdm$prediction
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
