
create_water_bodies_irreg <- function(boundary_sf, n_lakes = 2,
                                complexity = 8,
                                irregularity = 0.4,
                                elongation = c(0.6, 1.4)) {
   # Get the extent of the boundary
   ext <- st_bbox(boundary_sf)
   water_polys <- list()

   for(i in 1:n_lakes) {
      # Random center within boundary
      # The 0.2 value creates a 20% buffer around the edges
      # center_x <- runif(1, ext["xmin"] + (ext["xmax"] - ext["xmin"]) * 0.2,
      #                   ext["xmax"] - (ext["xmax"] - ext["xmin"]) * 0.2)
      # center_y <- runif(1, ext["ymin"] + (ext["ymax"] - ext["ymin"]) * 0.2,
      #                   ext["ymax"] - (ext["ymax"] - ext["ymin"]) * 0.2)
      # Sample points within the boundary with a buffer
      buffer_distance <- -1 * (ext["xmax"] - ext["xmin"]) * 0.2
      smaller_boundary <- st_buffer(boundary_sf, buffer_distance)
      random_point <- st_sample(smaller_boundary, 1)
      point_coords <- st_coordinates(random_point)
      center_x <- point_coords[1]
      center_y <- point_coords[2]

      # Base radius for this lake
      # Generate a random radius between 5% and 10% of the study area
      base_radius <- runif(1, (ext["xmax"] - ext["xmin"]) * 0.05,
                           (ext["xmax"] - ext["xmin"]) * 0.1)

      # Random elongation factor and rotation
      x_scale <- runif(1, elongation[1], elongation[2])
      y_scale <- runif(1, elongation[1], elongation[2])
      rotation_angle <- runif(1, 0, 2*pi)

      # Create more points for smoother outline
      n_points <- 30 + complexity
      angles <- seq(0, 2*pi, length.out = n_points)

      # Add irregular variations to the radius
      # Generate noise, with larger frequencies for more natural shapes
      noise_scale <- base_radius * irregularity
      radius_noise <- rnorm(n_points, 0, noise_scale)

      # Add some smoothness by applying moving average to the noise
      window_size <- 3
      smoothed_noise <- stats::filter(radius_noise, rep(1/window_size, window_size), circular = TRUE)

      # Calculate variable radius with the noise
      radius <- base_radius + smoothed_noise

      # Create coordinates with elongation and rotation
      basic_x <- radius * cos(angles)
      basic_y <- radius * sin(angles)

      # Apply elongation
      stretched_x <- basic_x * x_scale
      stretched_y <- basic_y * y_scale

      # Apply rotation
      rotated_x <- stretched_x * cos(rotation_angle) - stretched_y * sin(rotation_angle)
      rotated_y <- stretched_x * sin(rotation_angle) + stretched_y * cos(rotation_angle)

      # Translate to center
      final_x <- center_x + rotated_x
      final_y <- center_y + rotated_y

      # Create points
      lake_coords <- cbind(final_x, final_y)

      # Ensure the polygon is properly closed
      if(!identical(lake_coords[1,], lake_coords[nrow(lake_coords),])) {
         lake_coords <- rbind(lake_coords, lake_coords[1,])
      }

      # Create polygon
      water_poly <- tryCatch({
         st_polygon(list(lake_coords))
      }, error = function(e) {
         # Fallback to simpler shape if there's an error
         message("Error in lake ", i, ": ", e$message, ". Creating simpler shape.")
         angles <- seq(0, 2*pi, length.out = 20)
         circle_x <- center_x + base_radius * cos(angles)
         circle_y <- center_y + base_radius * sin(angles)
         lake_coords <- cbind(circle_x, circle_y)
         lake_coords <- rbind(lake_coords, lake_coords[1,])
         st_polygon(list(lake_coords))
      })

      # Check if water body is within boundary
      if(st_is_valid(water_poly) && st_intersects(st_sfc(water_poly, crs = st_crs(boundary_sf)),
                                                  boundary_sf, sparse = FALSE)[1,1]) {
         water_polys[[i]] <- water_poly
      } else {
         # If not valid or outside boundary, create simpler shape
         angles <- seq(0, 2*pi, length.out = 20)
         circle_x <- center_x + base_radius * 0.7 * cos(angles)
         circle_y <- center_y + base_radius * 0.7 * sin(angles)
         lake_coords <- cbind(circle_x, circle_y)
         lake_coords <- rbind(lake_coords, lake_coords[1,])
         water_polys[[i]] <- st_polygon(list(lake_coords))
      }
   }

   # Create sf object
   water_sf <- st_sf(geometry = st_sfc(water_polys), id = 1:length(water_polys))
   st_crs(water_sf) <- st_crs(boundary_sf)

   return(water_sf)
}
