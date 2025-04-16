# Enhanced Disease Risk Simulation Dashboard App
# Including environmental rasters and distribution visualizations

# Load required packages
library(shiny)
library(shinydashboard)
library(sf)
library(terra)
library(tidyterra)
library(ggplot2)
library(dplyr)
library(spsurvey)
library(shinyjs)
library(shinycssloaders)

# Source the simulation functions
source("R/simulation_functions.R")

# UI
ui <- dashboardPage(
   dashboardHeader(title = "Disease Risk Simulation"),
   dashboardSidebar(
      sidebarMenu(
         menuItem("Simulation Parameters", tabName = "params", icon = icon("sliders")),
         menuItem("Landscape", tabName = "landscape", icon = icon("map")),
         menuItem("Environment", tabName = "environment", icon = icon("tree")),
         menuItem("Animal Data", tabName = "animals", icon = icon("paw")),
         menuItem("Midge Data", tabName = "midges", icon = icon("bug")),
         menuItem("Risk Assessment", tabName = "risk", icon = icon("exclamation-triangle")),
         menuItem("About", tabName = "about", icon = icon("info-circle"))
      )
   ),
   dashboardBody(
      useShinyjs(),
      tabItems(
         # Parameters tab
         tabItem(tabName = "params",
                 fluidRow(
                    box(
                       title = "Landscape Parameters", width = 6, status = "primary",
                       sliderInput("water_bodies", "Number of Water Bodies:",
                                   min = 1, max = 8, value = 4, step = 1),
                       sliderInput("n_feeders", "Number of Feeders:",
                                   min = 5, max = 15, value = 10, step = 1)
                    ),
                    box(
                       title = "Animal Movement Parameters", width = 6, status = "primary",
                       sliderInput("n_animals", "Number of Animals:",
                                   min = 1, max = 10, value = 5, step = 1),
                       sliderInput("n_steps", "Number of Movement Steps:",
                                   min = 50, max = 300, value = 200, step = 50)
                    )
                 ),
                 fluidRow(
                    box(
                       width = 12, status = "success",
                       actionButton("run_simulation", "Run Simulation",
                                    icon = icon("play"),
                                    style = "color: #fff; background-color: #28a745; border-color: #28a745; width: 100%; height: 60px; font-size: 24px")
                    )
                 ),
                 fluidRow(
                    box(
                       title = "Note", width = 12, status = "info",
                       "This simulation uses default values for parameters not shown above. The simulation process may take a moment to complete."
                    )
                 )
         ),

         # Landscape tab
         tabItem(tabName = "landscape",
                 fluidRow(
                    box(
                       title = "Simulated Landscape", width = 12, status = "primary",
                       withSpinner(plotOutput("landscape_plot", height = "600px")),
                       downloadButton("download_landscape", "Download Plot")
                    )
                 )
         ),

         # Environment tab
         tabItem(tabName = "environment",
                 fluidRow(
                    box(
                       title = "Environmental Variables", width = 12, status = "primary",
                       selectInput("env_layer", "Select Environmental Layer:",
                                   choices = c("Elevation" = "elevation",
                                               "Water Distance" = "water_dist",
                                               "Vegetation Index" = "veg_index",
                                               "Temperature" = "temperature"),
                                   selected = "elevation"),
                       withSpinner(plotOutput("env_plot", height = "600px")),
                       downloadButton("download_env", "Download Plot")
                    )
                 ),
                 fluidRow(
                    box(
                       title = "Environmental Variables Description", width = 12, status = "info",
                       HTML("<ul>
                         <li><strong>Elevation:</strong> Simulated terrain elevation across the study area.</li>
                         <li><strong>Water Distance:</strong> Distance from each point to the nearest water body.</li>
                         <li><strong>Vegetation Index:</strong> Simulated vegetation density (similar to NDVI).</li>
                         <li><strong>Temperature:</strong> Simulated temperature patterns across the landscape.</li>
                       </ul>")
                    )
                 )
         ),

         # Animal Data tab
         tabItem(tabName = "animals",
                 fluidRow(
                    box(
                       title = "Animal Movement Tracks", width = 6, status = "primary",
                       withSpinner(plotOutput("movement_plot", height = "400px")),
                       downloadButton("download_movement", "Download Plot")
                    ),
                    box(
                       title = "Animal Utilization Distribution", width = 6, status = "primary",
                       withSpinner(plotOutput("ud_plot", height = "400px")),
                       downloadButton("download_ud", "Download Plot")
                    )
                 ),
                 fluidRow(
                    box(
                       title = "Animal Data Interpretation", width = 12, status = "info",
                       HTML("<p>The left plot shows the movement tracks of individual animals in the simulation. Each color represents a different animal.</p>
                       <p>The right plot displays the animal utilization distribution, which shows the probability of finding animals in different areas of the landscape.
                          Warmer colors (red, orange) indicate areas with higher animal activity.</p>")
                    )
                 )
         ),

         # Midge Data tab
         tabItem(tabName = "midges",
                 fluidRow(
                    box(
                       title = "Midge Sample Points", width = 6, status = "primary",
                       withSpinner(plotOutput("midge_sample_plot", height = "400px")),
                       downloadButton("download_midge_sample", "Download Plot")
                    ),
                    box(
                       title = "Predicted Midge Distribution", width = 6, status = "primary",
                       withSpinner(plotOutput("midge_prediction_plot", height = "400px")),
                       downloadButton("download_midge_prediction", "Download Plot")
                    )
                 ),
                 fluidRow(
                    box(
                       title = "Midge Data Interpretation", width = 12, status = "info",
                       HTML("<p>The left plot shows the midge sampling points used to build the distribution model. Points are colored by presence (1) or absence (0) of midges.</p>
                       <p>The right plot displays the predicted probability of midge presence across the landscape based on environmental variables.
                          Blues indicate higher probability of midge presence.</p>")
                    )
                 )
         ),

         # Risk tab
         tabItem(tabName = "risk",
                 fluidRow(
                    box(
                       title = "Disease Risk Map", width = 12, status = "warning",
                       withSpinner(plotOutput("risk_plot", height = "600px")),
                       downloadButton("download_risk", "Download Plot")
                    )
                 ),
                 fluidRow(
                    box(
                       title = "Risk Interpretation", width = 12, status = "info",
                       HTML("<p>The disease risk map combines animal utilization distribution and midge distribution to identify areas of potential disease transmission.</p>
                       <p>Areas with both high animal activity and high midge presence represent higher risk for disease transmission (shown in red).</p>
                       <p>This information can be used to target surveillance and control efforts in high-risk areas.</p>")
                    )
                 )
         ),

         # About tab
         tabItem(tabName = "about",
                 box(
                    title = "About This Simulator", width = 12,
                    HTML("<p>This application simulates animal movement and midge distribution to calculate disease risk maps.</p>
                     <p>The simulation follows these steps:</p>
                     <ol>
                       <li>Create a simulated landscape with water bodies and feeding stations</li>
                       <li>Generate environmental layers (elevation, water distance, vegetation, temperature)</li>
                       <li>Simulate animal movement within the landscape</li>
                       <li>Generate animal utilization distribution (where animals spend their time)</li>
                       <li>Simulate midge distribution based on environmental factors</li>
                       <li>Combine animal and midge distributions to create a disease risk map</li>
                     </ol>
                     <p>Adjust parameters in the 'Simulation Parameters' tab and click 'Run Simulation' to see results.</p>")
                 )
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
      midge_data = NULL,
      midge_sample = NULL,
      midge_prediction = NULL,
      risk_map = NULL
   )

   # Run simulation when button is clicked
   observeEvent(input$run_simulation, {

      # Disable run button during simulation
      disable("run_simulation")

      # Show progress notification
      withProgress(message = 'Running simulation...', value = 0, {

         # Step 1: Create study area
         incProgress(0.1, detail = "Creating landscape")
         sim_results$study_area <- create_study_area()

         # Step 2: Create water bodies
         incProgress(0.1, detail = "Adding water bodies")
         sim_results$water_bodies <- create_water_bodies(sim_results$study_area, n = input$water_bodies)

         # Step 3: Create feeders - with corrected parameters
         incProgress(0.1, detail = "Adding feeders")
         sim_results$feeders <- create_feeders(
            sim_results$study_area,
            input$n_feeders,
            sim_results$water_bodies
         )

         # Step 4: Create environmental rasters
         incProgress(0.1, detail = "Generating environmental layers")
         sim_results$env_rasters <- create_env_rasters(sim_results$study_area, sim_results$water_bodies)

         # Step 5: Simulate animal movement
         incProgress(0.2, detail = "Simulating animal movement")
         sim_results$animal_tracks <- simulate_animal_movement(
            study_area = sim_results$study_area,
            water_bodies = sim_results$water_bodies,
            feeders = sim_results$feeders,
            n_animals = input$n_animals,
            n_steps = input$n_steps
         )

         # Step 6: Create animal utilization distribution
         incProgress(0.1, detail = "Calculating animal utilization")
         sim_results$animal_ud <- create_animal_ud(
            sim_results$animal_tracks,
            sim_results$study_area,
            resolution = 10,
            smoothing_factor = 9
         )

         # Step 7: Simulate midge data
         incProgress(0.1, detail = "Simulating midge distribution")
         midge_sim <- simulate_midge_data(sim_results$env_rasters)
         sim_results$midge_data <- midge_sim$midge_data

         # Step 8: Get midge samples
         incProgress(0.1, detail = "Sampling midges")
         grts_sample <- spsurvey::grts(sim_results$midge_data, 100)
         sim_results$midge_sample <- grts_sample$sites_base %>%
            dplyr::select(presence, elevation, water_dist, veg_index, temperature, geometry)

         # Step 9: Fit midge distribution model
         incProgress(0.1, detail = "Modeling midge distribution")
         midge_sdm <- fit_midge_sdm(sim_results$midge_sample, sim_results$env_rasters)
         sim_results$midge_prediction <- midge_sdm$prediction

         # Step 10: Create risk map
         incProgress(0.1, detail = "Creating risk map")
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
         geom_sf(data = sim_results$feeders, fill = "orange", color = "red", size = 3) +
         labs(title = "Simulated Landscape",
              subtitle = paste(input$water_bodies, "water bodies and", input$n_feeders, "feeders")) +
         theme_minimal()
   })

   # Environmental raster plot
   output$env_plot <- renderPlot({
      req(sim_results$env_rasters, input$env_layer)

      selected_raster <- sim_results$env_rasters[[input$env_layer]]

      # Create a nicer title based on selection
      titles <- list(
         elevation = "Elevation (m)",
         water_dist = "Distance to Water Bodies (m)",
         veg_index = "Vegetation Index",
         temperature = "Temperature (°C)"
      )

      # Color schemes for different variables
      color_schemes <- list(
         elevation = terrain.colors(100),
         water_dist = colorRampPalette(c("darkblue", "lightblue", "white"))(100),
         veg_index = colorRampPalette(c("brown", "yellow", "darkgreen"))(100),
         temperature = heat.colors(100)
      )

      plot(selected_raster,
           main = titles[[input$env_layer]],
           col = color_schemes[[input$env_layer]])
   })

   # Animal movement plot
   output$movement_plot <- renderPlot({
      req(sim_results$animal_tracks)

      ggplot() +
         geom_sf(data = sim_results$study_area, fill = "white", color = "black") +
         geom_sf(data = sim_results$water_bodies, fill = "lightblue", color = "blue") +
         geom_sf(data = sim_results$feeders, fill = "orange", color = "red", size = 3) +
         geom_sf(data = sim_results$animal_tracks, aes(color = factor(animal_id)), size = 0.5) +
         labs(title = "Animal Movement Tracks",
              subtitle = paste(input$n_animals, "animals with", input$n_steps, "steps")) +
         scale_color_viridis_d(name = "Animal ID") +
         theme_minimal()
   })

   # Animal utilization distribution plot
   output$ud_plot <- renderPlot({
      req(sim_results$animal_ud)

      plot(sim_results$animal_ud,
           main = "Animal Utilization Distribution",
           col = colorRampPalette(c("white", "yellow", "orange", "red"))(100))
   })

   # Midge sample points plot
   output$midge_sample_plot <- renderPlot({
      req(sim_results$midge_sample)

      ggplot() +
         geom_sf(data = sim_results$study_area, fill = "white", color = "black") +
         geom_sf(data = sim_results$water_bodies, fill = "lightblue", color = "blue") +
         geom_sf(data = sim_results$midge_sample, aes(color = factor(presence)), size = 2) +
         scale_color_manual(values = c("0" = "gray", "1" = "darkblue"),
                            name = "Midge Presence",
                            labels = c("0" = "Absent", "1" = "Present")) +
         labs(title = "Midge Sampling Points",
              subtitle = "GRTS sampling design with presence/absence data") +
         theme_minimal()
   })

   # Midge prediction plot
   output$midge_prediction_plot <- renderPlot({
      req(sim_results$midge_prediction)

      plot(sim_results$midge_prediction,
           main = "Predicted Midge Distribution",
           col = colorRampPalette(c("white", "lightblue", "blue", "darkblue"))(100))
   })

   # Risk map plot
   output$risk_plot <- renderPlot({
      req(sim_results$risk_map)

      plot(sim_results$risk_map,
           main = "Disease Risk Map",
           col = colorRampPalette(c("blue", "green", "yellow", "orange", "red"))(100))
   })

   # Download handlers for plots
   output$download_landscape <- downloadHandler(
      filename = function() {
         paste("landscape-plot-", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
         ggsave(file, plot = {
            ggplot() +
               geom_sf(data = sim_results$study_area, fill = "white", color = "black") +
               geom_sf(data = sim_results$water_bodies, fill = "lightblue", color = "blue") +
               geom_sf(data = sim_results$feeders, fill = "orange", color = "red", size = 3) +
               labs(title = "Simulated Landscape",
                    subtitle = paste(input$water_bodies, "water bodies and", input$n_feeders, "feeders")) +
               theme_minimal()
         }, width = 10, height = 8, dpi = 300)
      }
   )

   output$download_env <- downloadHandler(
      filename = function() {
         paste("environment-", input$env_layer, "-", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
         png(file, width = 3000, height = 2400, res = 300)
         selected_raster <- sim_results$env_rasters[[input$env_layer]]

         titles <- list(
            elevation = "Elevation (m)",
            water_dist = "Distance to Water Bodies (m)",
            veg_index = "Vegetation Index",
            temperature = "Temperature (°C)"
         )

         color_schemes <- list(
            elevation = terrain.colors(100),
            water_dist = colorRampPalette(c("darkblue", "lightblue", "white"))(100),
            veg_index = colorRampPalette(c("brown", "yellow", "darkgreen"))(100),
            temperature = heat.colors(100)
         )

         plot(selected_raster,
              main = titles[[input$env_layer]],
              col = color_schemes[[input$env_layer]])
         dev.off()
      }
   )

   output$download_movement <- downloadHandler(
      filename = function() {
         paste("animal-movement-", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
         ggsave(file, plot = {
            ggplot() +
               geom_sf(data = sim_results$study_area, fill = "white", color = "black") +
               geom_sf(data = sim_results$water_bodies, fill = "lightblue", color = "blue") +
               geom_sf(data = sim_results$feeders, fill = "orange", color = "red", size = 3) +
               geom_sf(data = sim_results$animal_tracks, aes(color = factor(animal_id)), size = 0.5) +
               labs(title = "Animal Movement Tracks",
                    subtitle = paste(input$n_animals, "animals with", input$n_steps, "steps")) +
               scale_color_viridis_d(name = "Animal ID") +
               theme_minimal()
         }, width = 10, height = 8, dpi = 300)
      }
   )

   output$download_ud <- downloadHandler(
      filename = function() {
         paste("animal-utilization-", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
         png(file, width = 3000, height = 2400, res = 300)
         plot(sim_results$animal_ud,
              main = "Animal Utilization Distribution",
              col = colorRampPalette(c("white", "yellow", "orange", "red"))(100))
         dev.off()
      }
   )

   output$download_midge_sample <- downloadHandler(
      filename = function() {
         paste("midge-samples-", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
         ggsave(file, plot = {
            ggplot() +
               geom_sf(data = sim_results$study_area, fill = "white", color = "black") +
               geom_sf(data = sim_results$water_bodies, fill = "lightblue", color = "blue") +
               geom_sf(data = sim_results$midge_sample, aes(color = factor(presence)), size = 2) +
               scale_color_manual(values = c("0" = "gray", "1" = "darkblue"),
                                  name = "Midge Presence",
                                  labels = c("0" = "Absent", "1" = "Present")) +
               labs(title = "Midge Sampling Points",
                    subtitle = "GRTS sampling design with presence/absence data") +
               theme_minimal()
         }, width = 10, height = 8, dpi = 300)
      }
   )

   output$download_midge_prediction <- downloadHandler(
      filename = function() {
         paste("midge-prediction-", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
         png(file, width = 3000, height = 2400, res = 300)
         plot(sim_results$midge_prediction,
              main = "Predicted Midge Distribution",
              col = colorRampPalette(c("white", "lightblue", "blue", "darkblue"))(100))
         dev.off()
      }
   )

   output$download_risk <- downloadHandler(
      filename = function() {
         paste("risk-map-", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
         png(file, width = 3000, height = 2400, res = 300)
         plot(sim_results$risk_map,
              main = "Disease Risk Map",
              col = colorRampPalette(c("blue", "green", "yellow", "orange", "red"))(100))
         dev.off()
      }
   )
}

# Run the application
shinyApp(ui = ui, server = server)
