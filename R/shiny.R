# Disease Risk Simulation Shiny App - Standard Theme
# A simplified app focusing on key parameters with standard Shiny styling

# Load required packages
library(shiny)
library(sf)
library(terra)
library(tidyterra)
library(ggplot2)
library(dplyr)
library(spsurvey)
library(shinyjs)

# Source the simulation functions
source("R/simulation_functions.R")

# UI
ui <- fluidPage(
   useShinyjs(),

   # Title
   titlePanel("Disease Risk Simulation"),

   # Layout
   sidebarLayout(
      sidebarPanel(
         # Key parameters to modify
         sliderInput("water_bodies", "Number of Water Bodies:",
                     min = 1, max = 8, value = 4, step = 1),

         sliderInput("n_feeders", "Number of Feeders:",
                     min = 5, max = 15, value = 10, step = 1),

         sliderInput("n_animals", "Number of Animals:",
                     min = 1, max = 10, value = 5, step = 1),

         sliderInput("n_steps", "Number of Movement Steps:",
                     min = 50, max = 300, value = 200, step = 50),

         actionButton("run_simulation", "Run Simulation",
                      class = "btn-primary",
                      style = "width: 100%; margin-top: 20px;"),

         hr(),
         helpText("Note: Other parameters are set to default values")
      ),

      mainPanel(
         tabsetPanel(
            tabPanel("Landscape", plotOutput("landscape_plot")),
            tabPanel("Animal Movement", plotOutput("movement_plot")),
            tabPanel("Risk Map", plotOutput("risk_plot"))
         )
      )
   )
)

# Server
server <- function(input, output, session) {

   # Reactive values to store simulation results
   sim_results <- reactiveValues(
      study_area = NULL,
      water_bodies = NULL,
      feeders = NULL,
      env_rasters = NULL,
      animal_tracks = NULL,
      animal_ud = NULL,
      midge_prediction = NULL,
      risk_map = NULL
   )

   # Run simulation when button is clicked
   observeEvent(input$run_simulation, {

      # Disable run button during simulation
      disable("run_simulation")

      # Show progress notification
      withProgress(message = 'Running simulation', value = 0, {

         # Step 1: Create study area (fixed parameters)
         incProgress(0.1, detail = "Creating landscape")
         sim_results$study_area <- create_study_area()

         # Step 2: Create water bodies
         incProgress(0.2, detail = "Adding water bodies")
         sim_results$water_bodies <- create_water_bodies(sim_results$study_area, n = input$water_bodies)

         # Step 3: Create feeders
         incProgress(0.3, detail = "Adding feeders")
         # Using only the parameters that exist in your function
         sim_results$feeders <- create_feeders(
            sim_results$study_area,
            input$n_feeders,
            sim_results$water_bodies
         )

         # Step 4: Create environmental rasters
         incProgress(0.4, detail = "Generating environment")
         sim_results$env_rasters <- create_env_rasters(sim_results$study_area, sim_results$water_bodies)

         # Step 5: Simulate animal movement
         incProgress(0.5, detail = "Simulating animal movement")
         sim_results$animal_tracks <- simulate_animal_movement(
            study_area = sim_results$study_area,
            water_bodies = sim_results$water_bodies,
            feeders = sim_results$feeders,
            n_animals = input$n_animals,
            n_steps = input$n_steps
         )

         # Step 6: Create animal utilization distribution
         incProgress(0.6, detail = "Calculating animal utilization")
         sim_results$animal_ud <- create_animal_ud(
            sim_results$animal_tracks,
            sim_results$study_area,
            resolution = 10,
            smoothing_factor = 9
         )

         # Step 7: Simulate midge data
         incProgress(0.7, detail = "Simulating midges")
         midge_sim <- simulate_midge_data(sim_results$env_rasters)
         sim_results$midge_data <- midge_sim$midge_data

         # Step 8: Get midge samples (using GRTS with default parameters)
         incProgress(0.8, detail = "Sampling midges")
         grts_sample <- spsurvey::grts(sim_results$midge_data, 100)
         sim_results$midge_sample <- grts_sample$sites_base %>%
            dplyr::select(presence, elevation, water_dist, veg_index, temperature, geometry)

         # Step 9: Fit midge distribution model
         incProgress(0.9, detail = "Modeling midge distribution")
         midge_sdm <- fit_midge_sdm(sim_results$midge_sample, sim_results$env_rasters)
         sim_results$midge_prediction <- midge_sdm$prediction

         # Step 10: Create risk map
         incProgress(1.0, detail = "Creating risk map")
         sim_results$risk_map <- create_risk_map(sim_results$animal_ud, sim_results$midge_prediction)

      })

      # Re-enable run button
      enable("run_simulation")
   })

   # Landscape plot
   output$landscape_plot <- renderPlot({
      req(sim_results$study_area, sim_results$water_bodies, sim_results$feeders)

      ggplot() +
         geom_sf(data = sim_results$study_area, fill = "white", color = "black") +
         geom_sf(data = sim_results$water_bodies, fill = "lightblue", color = "blue") +
         geom_sf(data = sim_results$feeders, color = "red", size = 3) +
         labs(title = "Simulated Landscape",
              subtitle = paste(input$water_bodies, "water bodies and", input$n_feeders, "feeders")) +
         theme_minimal()
   })

   # Animal movement plot
   output$movement_plot <- renderPlot({
      req(sim_results$animal_tracks)

      ggplot() +
         geom_sf(data = sim_results$study_area, fill = "white", color = "black") +
         geom_sf(data = sim_results$water_bodies, fill = "lightblue", color = "blue") +
         geom_sf(data = sim_results$feeders, color = "red", size = 3) +
         geom_sf(data = sim_results$animal_tracks, aes(color = factor(animal_id)), size = 0.5) +
         labs(title = "Animal Movement Tracks",
              subtitle = paste(input$n_animals, "animals with", input$n_steps, "steps")) +
         scale_color_viridis_d(name = "Animal ID") +
         theme_minimal()
   })

   # Risk map plot
   output$risk_plot <- renderPlot({
      req(sim_results$risk_map)

      terra::plot(sim_results$risk_map, main = "Disease Risk Map", col = hcl.colors(50, "Spectral", rev = TRUE))
   })
}

# Run the application
shinyApp(ui = ui, server = server)
