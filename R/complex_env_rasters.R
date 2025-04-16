# Updated Temperature Range:
#
#    Changed the temperature scale to Fahrenheit as requested
# Open areas now range from 80-100°F (south to north gradient)
# Areas with dense canopy can be as cool as 65°F (35°F cooling effect)
# Added some random variation to create natural temperature patterns
# Named the layer "temperature_f" to indicate the units
#
#
# Added Humidity Layer:
#
#    Created a relative humidity layer (0-100%)
# Humidity is heavily influenced by canopy cover:
#
#    Sparse/no canopy areas: around 40% humidity baseline
# Dense forest areas: up to 90% humidity baseline
#
#
# Proximity to water bodies increases humidity (up to 20% boost very close to water)
# Added random variation to create natural humidity patterns
# Combined effects ensure realistic correlations between humidity, forest cover, and water proximity
#
# Maintained Other Features:
#
#    Distance to water bodies layer works the same as before
# Clustered canopy cover with user-specified number of patches
# All layers maintain the same resolution and spatial extent
#
#
#
# This implementation creates a realistic ecological gradient where:
#
#    Forested areas are cooler and more humid
# Open areas are hotter and drier
# Areas near water have increased humidity
# Temperature varies along a north-south gradient
# All environmental variables have appropriate spatial correlation and realistic values


create_env_rasters <- function(study_area, water_bodies, canopy_clusters = 3, resolution = 10) {
   # Convert study area to a SpatRaster for extent
   r <- rast(ext(st_bbox(study_area)), resolution = resolution, crs = st_crs(study_area)$wkt)

   # Create environmental layers list
   env_layers <- list()

   # Layer 1: Distance to water bodies
   if (length(st_geometry(water_bodies)) > 0) {
      # Make sure water bodies have the same CRS as study area
      if (st_crs(water_bodies) != st_crs(study_area)) {
         water_bodies <- st_transform(water_bodies, st_crs(study_area))
      }

      # Rasterize water bodies
      water_rast <- rasterize(vect(water_bodies), r, field = 1, background = 0)

      # Calculate distance to water
      water_dist <- distance(water_rast)
      names(water_dist) <- "water_dist"
      env_layers$water_dist <- water_dist
   } else {
      # If no water bodies, create uniform "far from water" layer
      water_dist <- rast(r)
      values(water_dist) <- rep(1000, ncell(water_dist))  # Arbitrary large distance
      names(water_dist) <- "water_dist"
      env_layers$water_dist <- water_dist
   }

   # Layer 2: Canopy cover with specified number of clusters
   canopy_cover <- rast(r)

   # Create a base for the clusters
   nr <- nrow(canopy_cover)
   nc <- ncol(canopy_cover)

   # Generate random cluster centers
   cluster_centers <- data.frame(
      x = runif(canopy_clusters, 1, nc),
      y = runif(canopy_clusters, 1, nr),
      strength = runif(canopy_clusters, 0.5, 1)  # How dense the canopy is in this cluster
   )

   # Fill the canopy cover raster based on distance to cluster centers
   m <- matrix(0, nrow = nr, ncol = nc)

   for (i in 1:nr) {
      for (j in 1:nc) {
         # Calculate influence from each cluster center
         for (k in 1:canopy_clusters) {
            # Calculate distance to this cluster center
            dist_to_center <- sqrt((j - cluster_centers$x[k])^2 + (i - cluster_centers$y[k])^2)

            # Calculate influence based on distance and strength
            # Further from center = less canopy cover
            # Use Gaussian-like falloff
            radius <- min(nc, nr) / 4  # Control the size of forest patches
            influence <- cluster_centers$strength[k] * exp(-(dist_to_center/radius)^2)

            # Add influence to this cell
            m[i, j] <- max(m[i, j], influence)  # Take max influence from any cluster
         }
      }
   }

   # Scale to 0-1 range
   m <- (m - min(m)) / (max(m) - min(m))

   # Convert to canopy percentage (0-100%)
   m <- m * 100

   values(canopy_cover) <- as.vector(m)
   names(canopy_cover) <- "canopy_cover"
   env_layers$canopy_cover <- canopy_cover

   # Layer 3: Temperature gradient (south to north with some noise)
   temperature <- rast(r)

   # Create coordinates matrix
   coords <- crds(temperature)

   # Normalize y-coordinates (north-south)
   y_coords <- coords[, 2]
   y_range <- range(y_coords)
   y_normalized <- (y_coords - y_range[1]) / (y_range[2] - y_range[1])

   # Get canopy effect for each cell
   canopy_effect <- m[cbind(
      pmin(nr, pmax(1, round((coords[, 2] - y_range[1]) / (y_range[2] - y_range[1]) * nr))),
      pmin(nc, pmax(1, round((coords[, 1] - ext(r)[1]) / (ext(r)[2] - ext(r)[1]) * nc)))
   )] / 100  # Convert percentage to proportion

   # Base temperature gradient: 100°F in south to 80°F in north in open areas
   base_temp <- 100 - 20 * y_normalized

   # Apply cooling effect of canopy (up to 35°F cooler in dense forest)
   # This ensures the minimum is around 65°F in dense forest areas
   canopy_cooling <- 35 * canopy_effect

   # Add some random noise
   noise <- rnorm(length(y_normalized), 0, 2)

   # Calculate final temperature in Fahrenheit
   temp_values <- base_temp - canopy_cooling + noise

   # Fill temperature raster
   temperature[] <- temp_values
   names(temperature) <- "temperature_f"
   env_layers$temperature_f <- temperature

   # Layer 4: Humidity (correlated with canopy cover and proximity to water)
   humidity <- rast(r)

   # Base humidity from canopy cover (forests have higher humidity)
   # Forest effect: 40-90% humidity depending on canopy density
   forest_humidity <- 40 + 50 * canopy_effect

   # Water proximity effect (higher humidity near water)
   if ("water_dist" %in% names(env_layers)) {
      # Get normalized water distance values (0-1)
      water_dist_values <- values(env_layers$water_dist)
      max_dist <- max(water_dist_values, na.rm = TRUE)
      if (max_dist > 0) {
         water_dist_norm <- water_dist_values / max_dist
      } else {
         water_dist_norm <- water_dist_values  # All zeros
      }

      # Water proximity increases humidity (up to 20% boost very close to water)
      water_humidity <- 20 * (1 - pmin(1, water_dist_norm / 0.3))
   } else {
      water_humidity <- 0
   }

   # Add some random variation
   humidity_noise <- rnorm(ncell(humidity), 0, 5)

   # Calculate final humidity (ensure it stays in 0-100% range)
   humidity_values <- pmin(100, pmax(0, forest_humidity + water_humidity + humidity_noise))

   # Fill humidity raster
   humidity[] <- humidity_values
   names(humidity) <- "relative_humidity"
   env_layers$relative_humidity <- humidity

   # Combine all layers into a SpatRaster
   env_stack <- rast(env_layers)
   return(env_stack)
}

create_env_rasters <- function(study_area, water_bodies, canopy_clusters = 3, resolution = 10) {
   # Convert study area to a SpatRaster for extent
   r <- rast(ext(st_bbox(study_area)), resolution = resolution, crs = st_crs(study_area)$wkt)

   # Create environmental layers list
   env_layers <- list()

   # Layer 1: Distance to water bodies
   if (length(st_geometry(water_bodies)) > 0) {
      # Make sure water bodies have the same CRS as study area
      if (st_crs(water_bodies) != st_crs(study_area)) {
         water_bodies <- st_transform(water_bodies, st_crs(study_area))
      }

      # Rasterize water bodies
      water_rast <- rasterize(vect(water_bodies), r)

      # Calculate distance to water
      water_dist <- distance(water_rast)
      names(water_dist) <- "water_dist"
      env_layers$water_dist <- water_dist
   } else {
      # If no water bodies, create uniform "far from water" layer
      water_dist <- rast(r)
      values(water_dist) <- rep(1000, ncell(water_dist))  # Arbitrary large distance
      names(water_dist) <- "water_dist"
      env_layers$water_dist <- water_dist
   }

   # Layer 2: Canopy cover with specified number of clusters
   canopy_cover <- rast(r)

   # Create a base for the clusters
   nr <- nrow(canopy_cover)
   nc <- ncol(canopy_cover)

   # Generate random cluster centers
   cluster_centers <- data.frame(
      x = runif(canopy_clusters, 1, nc),
      y = runif(canopy_clusters, 1, nr),
      strength = runif(canopy_clusters, 0.5, 1)  # How dense the canopy is in this cluster
   )

   # Fill the canopy cover raster based on distance to cluster centers
   m <- matrix(0, nrow = nr, ncol = nc)

   for (i in 1:nr) {
      for (j in 1:nc) {
         # Calculate influence from each cluster center
         for (k in 1:canopy_clusters) {
            # Calculate distance to this cluster center
            dist_to_center <- sqrt((j - cluster_centers$x[k])^2 + (i - cluster_centers$y[k])^2)

            # Calculate influence based on distance and strength
            # Further from center = less canopy cover
            # Use Gaussian-like falloff
            radius <- min(nc, nr) / 4  # Control the size of forest patches
            influence <- cluster_centers$strength[k] * exp(-(dist_to_center/radius)^2)

            # Add influence to this cell
            m[i, j] <- max(m[i, j], influence)  # Take max influence from any cluster
         }
      }
   }

   # Scale to 0-1 range
   m <- (m - min(m)) / (max(m) - min(m))

   # Convert to canopy percentage (0-100%)
   m <- m * 100

   values(canopy_cover) <- as.vector(m)
   names(canopy_cover) <- "canopy_cover"
   env_layers$canopy_cover <- canopy_cover

   # Layer 3: Temperature gradient (south to north with some noise)
   temperature <- rast(r)

   # Create coordinates matrix
   coords <- crds(temperature)

   # Normalize y-coordinates (north-south)
   y_coords <- coords[, 2]
   y_range <- range(y_coords)
   y_normalized <- (y_coords - y_range[1]) / (y_range[2] - y_range[1])

   # Get canopy effect for each cell
   canopy_effect <- m[cbind(
      pmin(nr, pmax(1, round((coords[, 2] - y_range[1]) / (y_range[2] - y_range[1]) * nr))),
      pmin(nc, pmax(1, round((coords[, 1] - ext(r)[1]) / (ext(r)[2] - ext(r)[1]) * nc)))
   )] / 100  # Convert percentage to proportion

   # Base temperature gradient: 100°F in south to 80°F in north in open areas
   base_temp <- 100 - 20 * y_normalized

   # Apply cooling effect of canopy (up to 35°F cooler in dense forest)
   # This ensures the minimum is around 65°F in dense forest areas
   canopy_cooling <- 35 * canopy_effect

   # Add some random noise (reduced noise for smoother pattern)
   noise <- rnorm(length(y_normalized), 0, 1.5)

   # Calculate final temperature in Fahrenheit
   temp_values <- base_temp - canopy_cooling + noise

   # Fill temperature raster
   temperature[] <- temp_values
   names(temperature) <- "temperature_f"
   env_layers$temperature_f <- temperature

   # Layer 4: Humidity (correlated with canopy cover and proximity to water)
   humidity <- rast(r)

   # Base humidity level throughout the landscape
   base_humidity <- 50  # 50% base humidity

   # Forest effect: additional 0-35% humidity depending on canopy density
   forest_humidity <- 35 * canopy_effect

   # Water proximity effect (more subtle)
   if ("water_dist" %in% names(env_layers)) {
      water_dist_values <- values(env_layers$water_dist)
      max_dist <- max(water_dist_values, na.rm = TRUE)

      if (max_dist > 0) {
         # Normalize distances
         water_dist_norm <- water_dist_values / max_dist

         # More gradual effect from water proximity (up to 10% boost, gradual falloff)
         # Using logistic decay for smoother transition
         water_effect_scale <- 0.2  # Controls how quickly effect drops with distance
         water_humidity <- 10 * (1 / (1 + exp(water_dist_norm / water_effect_scale - 4)))
      } else {
         water_humidity <- 0
      }
   } else {
      water_humidity <- 0
   }

   # Add subtle spatial autocorrelation instead of pure noise
   # Using a smoother pattern based on position
   humidity_smooth <- matrix(0, nrow = nr, ncol = nc)
   base_freqs <- matrix(rnorm(16), nrow = 4)  # Smaller base matrix for smoother variation

   # Interpolate to create smoother pattern
   for (i in 1:nr) {
      for (j in 1:nc) {
         # Map to positions in base_freqs
         bi <- 1 + (i - 1) / (nr - 1) * (nrow(base_freqs) - 1)
         bj <- 1 + (j - 1) / (nc - 1) * (ncol(base_freqs) - 1)

         # Bilinear interpolation for smooth pattern
         i1 <- floor(bi)
         i2 <- ceiling(bi)
         j1 <- floor(bj)
         j2 <- ceiling(bj)

         # Handle edge cases
         if (i2 > nrow(base_freqs)) i2 <- nrow(base_freqs)
         if (j2 > ncol(base_freqs)) j2 <- ncol(base_freqs)

         # Calculate weights
         wi <- bi - i1
         wj <- bj - j1

         # Bilinear interpolation
         if (i1 == i2 && j1 == j2) {
            humidity_smooth[i, j] <- base_freqs[i1, j1]
         } else if (i1 == i2) {
            humidity_smooth[i, j] <- (1 - wj) * base_freqs[i1, j1] + wj * base_freqs[i1, j2]
         } else if (j1 == j2) {
            humidity_smooth[i, j] <- (1 - wi) * base_freqs[i1, j1] + wi * base_freqs[i2, j1]
         } else {
            humidity_smooth[i, j] <- (1 - wi) * (1 - wj) * base_freqs[i1, j1] +
               wi * (1 - wj) * base_freqs[i2, j1] +
               (1 - wi) * wj * base_freqs[i1, j2] +
               wi * wj * base_freqs[i2, j2]
         }
      }
   }

   # Scale the smooth variation to have less impact
   humidity_smooth <- (humidity_smooth - min(humidity_smooth)) / (max(humidity_smooth) - min(humidity_smooth))
   humidity_variation <- 5 * (humidity_smooth - 0.5)  # Scale to +/- 2.5%

   # Calculate final humidity (ensure it stays in 0-100% range)
   humidity_values <- pmin(100, pmax(0, base_humidity + forest_humidity + water_humidity + humidity_variation))

   # Fill humidity raster
   humidity[] <- humidity_values
   names(humidity) <- "relative_humidity"
   env_layers$relative_humidity <- humidity

   # Combine all layers into a SpatRaster
   env_stack <- rast(env_layers)
   return(env_stack)
}

envs <- create_env_rasters(study_area, water_bodies)
plot(envs)
