create_risk_map <- function(
      animal_ud,                 # Animal utilization distribution
      vector_distribution,       # Disease vector distribution (e.g., midges)
      risk_method = "weighted",  # Method to calculate risk: "product", "weighted", "threshold", "exposure_time", "fuzzy"
      weights = c(animal = 0.5, vector = 0.5), # Weights for animal and vector components
      animal_threshold = 0.2,    # Threshold for animal presence (for threshold method)
      vector_threshold = 0.3,    # Threshold for vector presence (for threshold method)
      time_weight = 0.7,         # Weight for time component in exposure time method
      uncertain_buffer = 0.2,    # Buffer for fuzzy logic classification
      normalize = TRUE           # Whether to normalize final result
) {
   # Ensure rasters are aligned
   if (crs(animal_ud) != crs(vector_distribution) ||
       ext(animal_ud) != ext(vector_distribution) ||
       ncol(animal_ud) != ncol(vector_distribution) ||
       nrow(animal_ud) != nrow(vector_distribution)) {
      vector_distribution <- resample(vector_distribution, animal_ud)
   }

   # Calculate risk based on selected method
   if (risk_method == "product") {
      # Simple product (original method)
      risk_map <- animal_ud * vector_distribution

   } else if (risk_method == "weighted") {
      # Weighted combination - allows different importance for animal vs vector
      # Normalize weights if they don't sum to 1
      weights <- weights / sum(weights)
      risk_map <- (weights["animal"] * animal_ud) + (weights["vector"] * vector_distribution)

   } else if (risk_method == "threshold") {
      # Threshold-based - risk only where both exceed minimum thresholds
      animal_presence <- animal_ud > animal_threshold
      vector_presence <- vector_distribution > vector_threshold

      # Calculate conditional risk - high risk where both are present
      # Only consider vector density where animals are present
      risk_map <- vector_distribution * animal_presence

      # Optional: can further modify by animal density where vectors are present
      risk_map <- risk_map * (animal_ud / max(values(animal_ud), na.rm = TRUE)) * vector_presence

   } else if (risk_method == "exposure_time") {
      # Risk as function of exposure time (animal UD represents time spent)
      # Higher weight on time spent (animal UD) vs. vector probability

      # First normalize both inputs to 0-1 if they aren't already
      animal_norm <- (animal_ud - minmax(animal_ud)[1]) / (minmax(animal_ud)[2] - minmax(animal_ud)[1])
      vector_norm <- (vector_distribution - minmax(vector_distribution)[1]) /
         (minmax(vector_distribution)[2] - minmax(vector_distribution)[1])

      # Exposure time weighted more heavily than simple presence
      risk_map <- (time_weight * animal_norm) * ((1 - time_weight) + time_weight * vector_norm)

   } else if (risk_method == "fuzzy") {
      # Fuzzy logic approach - accounts for uncertainty in both distributions

      # Define fuzzy membership functions
      # High animal use areas (core)
      animal_core <- animal_ud >= (max(values(animal_ud), na.rm = TRUE) * (1 - uncertain_buffer))
      # Moderate animal use (uncertain)
      animal_uncertain <- animal_ud >= (max(values(animal_ud), na.rm = TRUE) * 0.3) &
         animal_ud < (max(values(animal_ud), na.rm = TRUE) * (1 - uncertain_buffer))

      # High vector probability
      vector_high <- vector_distribution >= (max(values(vector_distribution), na.rm = TRUE) * (1 - uncertain_buffer))
      # Moderate vector probability
      vector_uncertain <- vector_distribution >= (max(values(vector_distribution), na.rm = TRUE) * 0.3) &
         vector_distribution < (max(values(vector_distribution), na.rm = TRUE) * (1 - uncertain_buffer))

      # Create risk categories with fuzzy logic
      risk_map <- rast(animal_ud)
      values(risk_map) <- 0

      # Very high risk - core animal areas AND high vector probability
      very_high_risk <- animal_core & vector_high
      # High risk - core animal areas AND uncertain vector OR uncertain animal AND high vector
      high_risk <- (animal_core & vector_uncertain) | (animal_uncertain & vector_high)
      # Moderate risk - uncertain in both or high in one but low in other
      moderate_risk <- (animal_uncertain & vector_uncertain) |
         (animal_core & !vector_high & !vector_uncertain) |
         (vector_high & !animal_core & !animal_uncertain)

      # Assign risk values
      risk_map[very_high_risk] <- 1.0
      risk_map[high_risk] <- 0.7
      risk_map[moderate_risk] <- 0.4

      # Optional: can multiply by continuous values for more granularity
      risk_map <- risk_map * (0.5 + 0.5 * (animal_ud / max(values(animal_ud), na.rm = TRUE))) *
         (0.5 + 0.5 * (vector_distribution / max(values(vector_distribution), na.rm = TRUE)))

   } else {
      stop("Invalid risk_method. Choose from: 'product', 'weighted', 'threshold', 'exposure_time', or 'fuzzy'")
   }

   # Normalize final risk map if requested
   if (normalize) {
      risk_min <- minmax(risk_map)[1]
      risk_max <- minmax(risk_map)[2]

      # Check if all values are the same
      if (risk_min == risk_max) {
         risk_map <- risk_map / risk_max  # Set to 1 if no variation
      } else {
         risk_map <- (risk_map - risk_min) / (risk_max - risk_min)
      }

      # Add small constant to avoid zeros
      risk_map <- risk_map + 0.001
   }

   names(risk_map) <- "disease_risk"
   return(risk_map)
}

