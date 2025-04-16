create_irregular_study_area <- function(center_x = 0, center_y = 0,
                                        base_size = 100,
                                        irregularity = 0.2,
                                        complexity = 12) {
   # Generate points around a circle, but with varying distance from center
   n_points <- 30 + complexity
   angles <- seq(0, 2*pi, length.out = n_points)

   # Add irregular variations to the radius
   noise_scale <- base_size * irregularity
   radius_noise <- rnorm(n_points, 0, noise_scale)

   # Apply smoothing for more natural edges
   window_size <- 5
   smoothed_noise <- stats::filter(radius_noise, rep(1/window_size, window_size), circular = TRUE)

   # Calculate variable radius with the noise
   radius <- base_size + smoothed_noise

   # Create polygon coordinates
   poly_x <- center_x + radius * cos(angles)
   poly_y <- center_y + radius * sin(angles)

   # Create closed polygon
   poly_coords <- cbind(poly_x, poly_y)
   poly_coords <- rbind(poly_coords, poly_coords[1,])

   # Create polygon
   boundary_poly <- st_polygon(list(poly_coords))

   # Create sf object
   boundary_sf <- st_sf(geometry = st_sfc(boundary_poly))
   st_crs(boundary_sf) <- 4326  # Set to WGS84, change as needed

   return(boundary_sf)
}

# Alternative version that creates more diverse shapes with "bulges" and "inlets"
create_diverse_study_area <- function(center_x = 0, center_y = 0,
                                      base_size = 100,
                                      n_vertices = 20,
                                      roughness = 0.3,
                                      bulge_factor = 0.4) {
   # Create a base shape with random points
   angles <- sort(runif(n_vertices, 0, 2*pi))

   # Add variations to create more organic shapes
   # Some points will be pushed out (bulges), others pulled in (inlets)
   radius_variations <- numeric(n_vertices)

   # Create several bulges/inlets
   n_features <- round(n_vertices / 4)
   feature_centers <- sample(1:n_vertices, n_features)
   feature_widths <- rpois(n_features, 3) + 2
   feature_types <- sample(c(-1, 1), n_features, replace = TRUE, prob = c(0.4, 0.6))

   for(i in 1:n_features) {
      center <- feature_centers[i]
      width <- feature_widths[i]
      feature_type <- feature_types[i] * bulge_factor * base_size

      # Apply the feature to nearby points
      for(j in 1:n_vertices) {
         # Calculate distance along the circle (in vertices)
         dist <- min(abs(j - center), n_vertices - abs(j - center))
         if(dist <= width) {
            # Apply smoothed effect based on distance
            effect <- feature_type * cos(dist/width * pi/2)^2
            radius_variations[j] <- radius_variations[j] + effect
         }
      }
   }

   # Add small random noise
   fine_noise <- rnorm(n_vertices, 0, roughness * base_size)

   # Calculate final radius for each point
   radii <- base_size + radius_variations + fine_noise

   # Convert to coordinates
   poly_x <- center_x + radii * cos(angles)
   poly_y <- center_y + radii * sin(angles)

   # Create closed polygon
   poly_coords <- cbind(poly_x, poly_y)
   poly_coords <- rbind(poly_coords, poly_coords[1,])

   # Create polygon
   boundary_poly <- st_polygon(list(poly_coords))

   # Create sf object
   boundary_sf <- st_sf(geometry = st_sfc(boundary_poly))
   st_crs(boundary_sf) <- 4326  # Set to WGS84, change as needed

   return(boundary_sf)
}

# Example usage
# For a simple irregular area:
study_area <- create_irregular_study_area(base_size = 10000, irregularity = 0.15)
plot(study_area)

# For a more organic shape with distinct features:
complex_area <- create_diverse_study_area(base_size = 10000,
                                         n_vertices = 24,
                                         roughness = 0.2,
                                         bulge_factor = 0.5)
plot(complex_area)
