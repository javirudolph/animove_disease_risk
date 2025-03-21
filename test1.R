# Required packages
library(sf)
library(terra)
library(CircStats)  # For circular statistics (turning angles)
library(adehabitatHR)  # For kernel density estimation

# Load your simulated spatial data using sf
boundary <- st_read("simulated_data/boundary.shp", quiet = TRUE)
water <- st_read("simulated_data/water.shp", quiet = TRUE)
feeders <- st_read("simulated_data/feeders.shp", quiet = TRUE)

# Load your simulated environmental data using terra
temperature <- rast("simulated_data/temperature.tif")
humidity <- rast("simulated_data/humidity.tif")
vegetation <- rast("simulated_data/vegetation.tif")

# Convert to raster for easier calculations
r <- rast(ext = st_bbox(boundary), resolution = 10)  # Adjust resolution as needed
r <- rasterize(boundary, r, field = 1)
water_raster <- rasterize(water, r, field = 1)
feeder_raster <- rasterize(feeders, r, field = 1)

# Distance to feeders (for attraction)
feeder_dist <- distance(feeder_raster)

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

# Run simulation for multiple animals
n_animals <- 10
n_steps <- 1000
animal_trajectories <- list()

for(a in 1:n_animals) {
   # Random starting location within boundary
   valid_start <- FALSE
   while(!valid_start) {
      ext <- st_bbox(boundary)
      start_x <- runif(1, ext["xmin"], ext["xmax"])
      start_y <- runif(1, ext["ymin"], ext["ymax"])

      # Check if starting point is valid using point-in-polygon test
      point_sf <- st_sfc(st_point(c(start_x, start_y)), crs = st_crs(boundary))
      in_boundary <- st_intersects(point_sf, boundary, sparse = FALSE)[1,1]

      if(in_boundary) {
         # Check if not in water
         in_water <- any(st_intersects(point_sf, water, sparse = FALSE)[1,])
         if(!in_water) {
            valid_start <- TRUE
         }
      }
   }

   animal_trajectories[[a]] <- simulate_animal_movement(
      c(start_x, start_y),
      n_steps,
      step_length_mean = 50,  # Adjust for your animal
      step_length_sd = 20,
      turning_angle_kappa = 0.7  # Higher = more directional persistence
   )
}

# Create utilization distribution from trajectories
all_points <- do.call(rbind, animal_trajectories)
all_points_sf <- st_as_sf(data.frame(x = all_points[,1], y = all_points[,2]),
                          coords = c("x", "y"), crs = st_crs(boundary))

# Calculate KDE using adehabitatHR (works with sf objects)
all_points_sp <- as(all_points_sf, "Spatial")  # temporary conversion for adehabitatHR
ud <- kernelUD(all_points_sp, grid = 100)
ud_raster <- raster(ud)  # temporary raster object

# Convert to terra raster
ud_terra <- rast(ud_raster)

# Now, simulate the midge abundance and create the species distribution model
# Standardize variables
temp_mean <- global(temperature, fun = "mean", na.rm = TRUE)[1,1]
temp_sd <- global(temperature, fun = "sd", na.rm = TRUE)[1,1]
temp_std <- (temperature - temp_mean) / temp_sd

humid_mean <- global(humidity, fun = "mean", na.rm = TRUE)[1,1]
humid_sd <- global(humidity, fun = "sd", na.rm = TRUE)[1,1]
humid_std <- (humidity - humid_mean) / humid_sd

veg_mean <- global(vegetation, fun = "mean", na.rm = TRUE)[1,1]
veg_sd <- global(vegetation, fun = "sd", na.rm = TRUE)[1,1]
veg_std <- (vegetation - veg_mean) / veg_sd

# Create a "true" species distribution model
# Midges prefer warm, humid areas with moderate vegetation
midge_probability <- 0.7 * temp_std + 0.8 * humid_std - 0.3 * (veg_std^2) +
   0.3 * (temp_std * humid_std)

# Add spatial autocorrelation
library(fields)
# Create a spatial correlation field
r_template <- rast(ext = ext(midge_probability), resolution = res(midge_probability)[1])
pts <- as.data.frame(xyFromCell(r_template, 1:ncell(r_template)))

# Generate a spatial field with Matern correlation
set.seed(123)
sim_field <- sim.rf(pts, setup = TRUE, theta = 500, cov.args = list(Covariance = "Matern", range = 1000, smoothness = 1))

# Create a raster of the field
spatial_effect <- r_template
spatial_effect[] <- sim_field
spatial_effect <- scale(spatial_effect)  # standardize

# Combine effects
midge_probability <- scale(midge_probability + spatial_effect)
# Convert to probability scale using logistic transformation
midge_probability <- 1 / (1 + exp(-midge_probability))

# Simulate sampling locations
boundary_buffer <- st_buffer(boundary, dist = -100)  # Buffer to ensure points are within
sample_points <- st_sample(boundary_buffer, size = 100, type = "random")
sample_points_sf <- st_sf(geometry = sample_points)

# Extract environmental data
sample_data <- data.frame(
   x = st_coordinates(sample_points_sf)[,1],
   y = st_coordinates(sample_points_sf)[,2]
)

# Extract values at sample points
sample_data$temperature <- extract(temperature, sample_points_sf)[,1]
sample_data$humidity <- extract(humidity, sample_points_sf)[,1]
sample_data$vegetation <- extract(vegetation, sample_points_sf)[,1]

# Add midge counts based on probability
true_prob <- extract(midge_probability, sample_points_sf)[,1]
sample_data$count <- rpois(nrow(sample_data), lambda = true_prob * 50)  # Scale factor of 50

# Build SDM with the sampled data
library(mgcv)
sdm_model <- gam(count ~ s(temperature, k = 4) + s(humidity, k = 4) +
                    s(vegetation, k = 4) + s(x, y, k = 15),
                 family = poisson, data = sample_data)

# Create prediction raster stack
pred_stack <- c(temperature, humidity, vegetation)
names(pred_stack) <- c("temperature", "humidity", "vegetation")

# Create coordinates for prediction
coords <- crds(pred_stack)
pred_df <- as.data.frame(coords)
names(pred_df) <- c("x", "y")

# Extract environmental values
pred_df$temperature <- extract(temperature, coords)[,1]
pred_df$humidity <- extract(humidity, coords)[,1]
pred_df$vegetation <- extract(vegetation, coords)[,1]

# Make predictions with NA handling
pred_df <- na.omit(pred_df)
predicted_values <- predict(sdm_model, pred_df, type = "response")

# Create prediction raster
predicted_raster <- rast(temperature)  # Use as template
predicted_raster[] <- NA  # Reset values

# Fill in predicted values
cell_indices <- cellFromXY(predicted_raster, pred_df[, c("x", "y")])
predicted_raster[cell_indices] <- predicted_values

# Scale to 0-1 probability
min_val <- minmax(predicted_raster)[1]
max_val <- minmax(predicted_raster)[2]
predicted_prob <- (predicted_raster - min_val) / (max_val - min_val)

# Create risk map by multiplying animal utilization and midge probability
# Ensure both rasters have the same extent and resolution
ud_terra_resampled <- resample(ud_terra, predicted_prob)
risk_map <- ud_terra_resampled * predicted_prob

# Plot the results
par(mfrow = c(1, 3))
plot(ud_terra_resampled, main = "Animal Space Use")
plot(predicted_prob, main = "Midge Distribution")
plot(risk_map, main = "Disease Risk")

# Save results
writeRaster(risk_map, "simulated_data/disease_risk.tif", overwrite = TRUE)
writeRaster(ud_terra_resampled, "simulated_data/animal_utilization.tif", overwrite = TRUE)
writeRaster(predicted_prob, "simulated_data/midge_probability.tif", overwrite = TRUE)
