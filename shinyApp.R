# app.R - Disease Risk Simulation Shiny App
# Assumes all simulation functions are in a separate script called "simulation_functions.R"

# Load required packages
library(shiny)
library(shinydashboard)
library(sf)
library(terra)
library(tidyterra)
library(ggplot2)
library(dplyr)
library(spsurvey)
library(plotly)
library(shinycssloaders)

# Source the simulation functions
source("R/simulation_functions.R")

# UI
ui <- dashboardPage(
   dashboardHeader(title = "Disease Risk Simulation"),
   dashboardSidebar(
      sidebarMenu(
         menuItem("Simulation Parameters", tabName = "params", icon = icon("sliders")),
         menuItem("Simulation Results", tabName = "results", icon = icon("map")),
         menuItem("About", tabName = "about", icon = icon("info-circle"))
      )
   ),
   dashboardBody(
      tabItems(
         # Parameters tab
         tabItem(tabName = "params",
                 fluidRow(
                    box(
                       title = "Landscape Parameters", width = 6, status = "primary",
                       sliderInput("study_area_width", "Study Area Width (m):",
                                   min = 500, max = 2000, value = 1000, step = 100),
                       sliderInput("study_area_height", "Study Area Height (m):",
                                   min = 500, max = 2000, value = 1000, step = 100),
                       sliderInput("water_bodies", "Number of Water Bodies:",
                                   min = 1, max = 10, value = 4, step = 1),
                       sliderInput("n_feeders", "Number of Feeders:",
                                   min = 5, max = 20, value = 10, step = 1),
                       sliderInput("min_dist_to_water", "Min Distance from Feeders to Water (m):",
                                   min = 10, max = 200, value = 100, step = 10),
                       sliderInput("min_dist_between_feeders", "Min Distance Between Feeders (m):",
                                   min = 10, max = 200, value = 100, step = 10)
                    ),
                    box(
                       title = "Animal Movement Parameters", width = 6, status = "primary",
                       sliderInput("n_animals", "Number of Animals:",
                                   min = 1, max = 20, value = 5, step = 1),
                       sliderInput("n_steps", "Number of Movement Steps:",
                                   min = 50, max = 500, value = 200, step = 50),
                       sliderInput("resolution", "UD Resolution (m):",
                                   min = 5, max = 50, value = 10, step = 5),
                       sliderInput("smoothing_factor", "UD Smoothing Factor:",
                                   min = 1, max = 20, value = 9, step = 1)
                    )
                 ),
                 fluidRow(
                    box(
                       title = "Midge Sampling Parameters", width = 6, status = "primary",
                       radioButtons("sampling_method", "Sampling Method:",
                                    choices = c("Random" = "random", "GRTS" = "grts"),
                                    selected = "grts"),
                       sliderInput("n_samples", "Number of Midge Sampling Points:",
                                   min = 20, max = 200, value = 100, step = 10)
                    ),
                    box(
                       width = 6, status = "success",
                       actionButton("run_simulation", "Run Simulation",
                                    icon = icon("play"),
                                    style = "color: #fff; background-color: #28a745; border-color: #28a745; width: 100%; height: 60px; font-size: 24px")
                    )
                 )
         ),

         # Results tab
         tabItem(tabName = "results",
                 fluidRow(
                    tabBox(
                       title = "Simulation Results", width = 12,
                       tabPanel("Landscape",
                                withSpinner(plotOutput("landscape_plot", height = "600px")),
                                downloadButton("download_landscape", "Download Plot")),
                       tabPanel("Animal Movement",
                                withSpinner(plotOutput("movement_plot", height = "600px")),
                                downloadButton("download_movement", "Download Plot")),
                       tabPanel("Animal Utilization",
                                withSpinner(plotOutput("ud_plot", height = "600px")),
                                downloadButton("download_ud", "Download Plot")),
                       tabPanel("Midge Distribution",
                                withSpinner(plotOutput("midge_plot", height = "600px")),
                                downloadButton("download_midge", "Download Plot")),
                       tabPanel("Risk Map",
                                withSpinner(plotOutput("risk_plot", height = "600px")),
                                downloadButton("download_risk", "Download Plot"))
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

# check out the ui first
shinyApp(ui, function(...){})

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

      withProgress(message = 'Running simulation...', value = 0, {

         # Step 1: Create study area
         incProgress(0.1, detail = "Creating landscape")
         sim_results$study_area <- create_study_area(width = input$study_area_width, height = input$study_area_height)

         # Step 2: Create water bodies
         incProgress(0.1, detail = "Adding water bodies")
         sim_results$water_bodies <- create_water_bodies(sim_results$study_area, n = input$water_bodies)

         # Step 3: Create feeders
         incProgress(0.1, detail = "Adding feeders")
         sim_results$feeders <- create_feeders(
            sim_results$study_area,
            n = input$n_feeders,
            water_bodies = sim_results$water_bodies,
            min_dist_to_water = input$min_dist_to_water,
            min_dist_between = input$min_dist_between_feeders
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
            resolution = input$resolution,
            smoothing_factor = input$smoothing_factor
         )

         # Step 7: Simulate midge data
         incProgress(0.1, detail = "Simulating midge distribution")
         midge_sim <- simulate_midge_data(sim_results$env_rasters)
         sim_results$midge_data <- midge_sim$midge_data

         # Step 8: Get midge samples
         incProgress(0.1, detail = "Sampling midges")
         if (input$sampling_method == "random") {
            sim_results$midge_sample <- dplyr::slice_sample(sim_results$midge_data, n = input$n_samples)
         } else { # GRTS
            grts_sample <- spsurvey::grts(sim_results$midge_data, input$n_samples)
            sim_results$midge_sample <- grts_sample$sites_base %>%
               dplyr::select(presence, elevation, water_dist, veg_index, temperature, geometry)
         }

         # Step 9: Fit midge distribution model
         incProgress(0.1, detail = "Modeling midge distribution")
         midge_sdm <- fit_midge_sdm(sim_results$midge_sample, sim_results$env_rasters)
         sim_results$midge_prediction <- midge_sdm$prediction

         # Step 10: Create risk map
         incProgress(0.1, detail = "Creating risk map")
         sim_results$risk_map <- create_risk_map(sim_results$animal_ud, sim_results$midge_prediction)
      })

   })

   # Landscape plot
   output$landscape_plot <- renderPlot({
      req(sim_results$study_area, sim_results$water_bodies, sim_results$feeders)

      ggplot() +
         geom_sf(data = sim_results$study_area, fill = "white", color = "black") +
         geom_sf(data = sim_results$water_bodies, fill = "lightblue", color = "blue") +
         geom_sf(data = sim_results$feeders, fill = "orange", color = "red", size = 3) +
         labs(title = "Simulated Landscape",
              subtitle = "Study area with water bodies and feeders") +
         theme_minimal()
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
              subtitle = paste("Tracks for", input$n_animals, "animals with", input$n_steps, "steps")) +
         scale_color_viridis_d(name = "Animal ID") +
         theme_minimal()
   })

   # Animal utilization distribution plot
   output$ud_plot <- renderPlot({
      req(sim_results$animal_ud)

      plot(sim_results$animal_ud, main = "Animal Utilization Distribution", col = hcl.colors(50, "YlOrRd"))
   })

   # Midge distribution plot
   output$midge_plot <- renderPlot({
      req(sim_results$midge_prediction)

      plot(sim_results$midge_prediction, main = "Predicted Midge Distribution", col = hcl.colors(50, "Blues"))
   })

   # Risk map plot
   output$risk_plot <- renderPlot({
      req(sim_results$risk_map)

      plot(sim_results$risk_map, main = "Disease Risk Map", col = hcl.colors(50, "Spectral", rev = TRUE))
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
                    subtitle = "Study area with water bodies and feeders") +
               theme_minimal()
         }, width = 10, height = 8, dpi = 300)
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
                    subtitle = paste("Tracks for", input$n_animals, "animals with", input$n_steps, "steps")) +
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
         plot(sim_results$animal_ud, main = "Animal Utilization Distribution", col = hcl.colors(50, "YlOrRd"))
         dev.off()
      }
   )

   output$download_midge <- downloadHandler(
      filename = function() {
         paste("midge-distribution-", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
         png(file, width = 3000, height = 2400, res = 300)
         plot(sim_results$midge_prediction, main = "Predicted Midge Distribution", col = hcl.colors(50, "Blues"))
         dev.off()
      }
   )

   output$download_risk <- downloadHandler(
      filename = function() {
         paste("risk-map-", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
         png(file, width = 3000, height = 2400, res = 300)
         plot(sim_results$risk_map, main = "Disease Risk Map", col = hcl.colors(50, "Spectral", rev = TRUE))
         dev.off()
      }
   )
}

# Run the application
shinyApp(ui = ui, server = server)
