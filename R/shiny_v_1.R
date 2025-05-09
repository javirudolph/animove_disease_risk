# Final Disease Risk Simulation Dashboard App
# With enhanced visualizations and risk metrics

# Load required packages
# Fixed Disease Risk Simulation Dashboard App
# With error handling for midge distribution model

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
library(cowplot)

# Source the simulation functions
# source("R/simulation_functions.R") # for local use
# when publishing:
source("simulation_functions.R")

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
         menuItem("Risk Assessment", tabName = "risk", icon = icon("exclamation-triangle"))
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
                       sliderInput("study_area_width", "Study Area Width (m):",
                                   min = 500, max = 2000, value = 1000, step = 100),
                       sliderInput("study_area_height", "Study Area Height (m):",
                                   min = 500, max = 2000, value = 1000, step = 100),
                       sliderInput("water_bodies", "Number of Water Bodies:",
                                   min = 1, max = 10, value = 3, step = 1),
                       sliderInput("n_feeders", "Number of Feeders:",
                                   min = 1, max = 20, value = 8, step = 1),
                       sliderInput("buffer_feeders", "Min Distance Between Feeders (m):",
                                   min = 10, max = 200, value = 50, step = 10)
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
                       # Add the seed input above the action button
                       numericInput("seed", "Set random seed for reproducibility:",
                                    value = 123, min = 1),
                       # Action button remains the same
                       actionButton("run_simulation", "Run Simulation",
                                    icon = icon("play"),
                                    style = "color: #fff; background-color: #28a745; border-color: #28a745; width: 100%; height: 60px; font-size: 18px")
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

         # Environment tab (all 4 rasters together)
         tabItem(tabName = "environment",
                 fluidRow(
                    box(
                       title = "Environmental Variables", width = 12, status = "primary",
                       withSpinner(plotOutput("all_env_plots", height = "800px")),
                       downloadButton("download_all_env", "Download Plot")
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
                       HTML("<p>The left plot shows the locations of individual animals in the simulation. Each color represents a different animal.</p>
                       <p>The right plot displays the animal utilization distribution, which shows the probability of finding animals in different areas of the landscape. Areas with higher values indicate higher animal activity.</p>")
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
                       <p>The right plot displays the predicted probability of midge presence across the landscape based on environmental variables.</p>")
                    )
                 )
         ),

         # Risk tab with water bodies and feeders
         tabItem(tabName = "risk",
                 fluidRow(
                    box(
                       title = "Disease Risk Map", width = 6, status = "warning",
                       withSpinner(plotOutput("risk_plot", height = "600px")),
                       downloadButton("download_risk", "Download Plot")
                    ),
                    box(
                       title = "Risk Metrics", width = 4, status = "danger",
                       withSpinner(plotOutput("risk_histogram", height = "300px")),
                       downloadButton("download_risk_metrics", "Download Risk Metrics")
                    )
                 ),
                 fluidRow(
                    box(
                       title = "Risk Feeders", width = 10, status = "danger",
                       withSpinner(tableOutput("feeder_risk_table")),
                       downloadButton("download_feeder_risk", "Download Feeder Risk Table")
                    )
                 ),
                 fluidRow(
                    box(
                       title = "Risk Interpretation", width = 12, status = "info",
                       HTML("<p>The disease risk map combines animal utilization distribution and midge distribution to identify areas of potential disease transmission.</p>
                       <p>Areas with both high animal activity and high midge presence represent higher risk for disease transmission.</p>
                       <p>The risk histogram shows the distribution of risk values across the landscape, while the table ranks feeders by their associated risk level.</p>
                       <p>This information can be used to target surveillance and control efforts in high-risk areas.</p>")
                    )
                 )
         )
      )
   )
)

# shinyApp(ui, function(...){})

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
      risk_map = NULL,
      feeder_risk = NULL
   )

   # Run simulation when button is clicked
   observeEvent(input$run_simulation, {

      # Set the seed at the beginning of the simulation
      set.seed(input$seed)

      withProgress(message = 'Running simulation...', value = 0, {

         # Step 1: Create study area
         incProgress(0.1, detail = "Creating landscape")
         sim_results$study_area <- create_study_area(xmax = input$study_area_width, ymax = input$study_area_height)

         # Step 2: Create water bodies
         incProgress(0.1, detail = "Adding water bodies")
         sim_results$water_bodies <- create_water_bodies(sim_results$study_area, n = input$water_bodies)

         # Step 3: Create feeders
         incProgress(0.1, detail = "Adding feeders")
         sim_results$feeders <- create_feeders(
            sim_results$study_area,
            n = input$n_feeders,
            water_bodies = sim_results$water_bodies,
            buffer_dist = input$buffer_feeders
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
         midge_sim <- simulate_midge_data(sim_results$env_rasters, sim_results$water_bodies)
         sim_results$midge_data <- midge_sim$midge_data

         # Step 8: Get midge samples
         incProgress(0.1, detail = "Sampling midges")
         if (input$sampling_method == "random") {
            sim_results$midge_sample <- dplyr::slice_sample(sim_results$midge_data, n = input$n_samples)
         } else { # GRTS
            grts_sample <- spsurvey::grts(sim_results$midge_data, input$n_samples)
            sim_results$midge_sample <- grts_sample$sites_base |>
               dplyr::select(presence, elevation, water_dist, veg_index, temperature, geometry)
         }

         # Step 9: Fit midge distribution model
         incProgress(0.1, detail = "Modeling midge distribution")
         midge_sdm <- fit_midge_sdm(sim_results$midge_sample, sim_results$env_rasters)
         sim_results$midge_prediction <- midge_sdm$prediction

         # Step 10: Create risk map
         incProgress(0.1, detail = "Creating risk map")
         sim_results$risk_map <- create_risk_map(sim_results$animal_ud, sim_results$midge_prediction)

         # Step 11: Calculate feeder risk
         incProgress(0.1, detail = "Calculating feeder risk")
         sim_results$feeder_risk <- calculate_feeder_risk(
            risk_surface = sim_results$risk_map,
            feeders = sim_results$feeders,
            buffer_radius = 30,
            metrics = c("point", "mean", "max", "quantile"),
            quantile_probs = c(0.5, 0.75, 0.9, 0.95)
         )

      })

   })


   # Lanscape plot
   landscape_plot <- reactive({
      req(sim_results$study_area, sim_results$water_bodies, sim_results$feeders)

      ggplot() +
         geom_sf(data = sim_results$study_area, fill = NA, color = "black") +
         geom_sf(data = sim_results$water_bodies, aes(fill = "Water Bodies", color = "Water Bodies")) +
         geom_sf(data = sim_results$feeders, aes(fill = "Feeders", color = "Feeders")) +
         # Add labels for feeder IDs
         geom_sf_label(data = sim_results$feeders,
                       aes(label = id),
                       color = "white",
                       fill = "#D95F02",
                       size = 3,
                       fontface = "bold") +
         # Set manual color and fill scales
         scale_fill_manual(name = "Features",
                           values = c("Water Bodies" ="#81C3D7",
                                      "Feeders" = "#D95F02")) +
         scale_color_manual(name = "Features",
                            values = c("Water Bodies" = "#81C3D7",
                                       "Feeders" = "#D95F02")) +
         labs(title = "Simulated Landscape",
              subtitle = "Study area with water bodies and feeders") +
         theme_minimal() +
         theme(axis.title = element_blank())
   })

   # Use the reactive for display
   output$landscape_plot <- renderPlot({
      landscape_plot()
   })

   # All environmental rasters plot (2x2 grid)
   env_plots <- reactive({
      req(sim_results$env_rasters)

      # Function to create a raster plot with specific title using tidyterra
      plot_raster <- function(raster_name, title) {
         selected_raster <- sim_results$env_rasters[[raster_name]]

         ggplot() +
            geom_spatraster(data = selected_raster) +
            # scale_fill_gradient(low = "lightgrey", high = "darkblue", na.value = "transparent") +
            scale_fill_viridis_c() +
            labs(title = title, fill = "") +
            theme_minimal() +
            theme(legend.position = "bottom",
                  plot.title = element_text(size = 12, face = "bold"),
                  axis.title = element_blank())
      }

      # Create the four plots
      p1 <- plot_raster("elevation", "Elevation (m)")
      p2 <- plot_raster("water_dist", "Distance to Water (m)")
      p3 <- plot_raster("veg_index", "Vegetation Index")
      p4 <- plot_raster("temperature", "Temperature (°C)")

      # Arrange the plots in a 2x2 grid
      cowplot::plot_grid(p1, p2, p3, p4, ncol = 2)
   })

   # Render env plots
   output$all_env_plots <- renderPlot({
      env_plots()
   })


   # Animal movement plot
   movement_plot <- reactive({
      req(sim_results$animal_tracks)

      ggplot() +
         geom_sf(data = sim_results$study_area, fill = NA, color = "black") +
         geom_sf(data = sim_results$water_bodies, color ="#81C3D7", fill ="#81C3D7", alpha = 0.6) +
         # Add labels for feeder IDs
         geom_sf_label(data = sim_results$feeders,
                       aes(label = id),
                       color = "white",
                       fill = "#D95F02",
                       size = 3,
                       fontface = "bold") +
         geom_sf(data = sim_results$animal_tracks, aes(color = factor(animal_id)), size = 1) +
         labs(title = "Animal Movement Tracks",
              subtitle = paste("Tracks for", input$n_animals, "animals with", input$n_steps, "steps")) +
         scale_color_viridis_d(name = "Animal ID", option = "H") +
         theme_minimal()
   })

   # Render the plot
   output$movement_plot <- renderPlot({
      movement_plot()
   })

   # Animal utilization distribution plot
   ud_plot <- reactive({
      req(sim_results$animal_ud)

      ggplot() +
         # Add the raster layer with geom_spatraster from tidyterra
         geom_spatraster(data = sim_results$animal_ud) +
         geom_sf(data = sim_results$study_area, fill = NA, color = "black") +
         geom_sf(data = sim_results$water_bodies, color ="#81C3D7", fill ="#81C3D7", alpha = 0.6) +
         # Add labels for feeder IDs
         geom_sf_label(data = sim_results$feeders,
                       aes(label = id),
                       color = "white",
                       fill = "#D95F02",
                       size = 3,
                       fontface = "bold") +
         scale_fill_gradientn(colors = colorRampPalette(c("white", "yellow", "orange", "red"))(100),
                              name = "Density") +
         labs(title = "Animal Utilization Distribution") +
         theme_minimal()

   })

   # Render the plot
   output$ud_plot <- renderPlot({
      ud_plot()
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

      # Get raster as a data frame for ggplot
      midge_df <- as.data.frame(sim_results$midge_prediction, xy = TRUE)
      colnames(midge_df)[3] <- "value"

      ggplot(midge_df, aes(x = x, y = y, fill = value)) +
         geom_raster() +
         geom_sf(data = sim_results$water_bodies, fill = "lightblue", color = "blue", inherit.aes = FALSE) +
         geom_sf(data = sim_results$feeders, fill = "orange", color = "red", size = 3, inherit.aes = FALSE) +
         scale_fill_gradientn(colors = colorRampPalette(c("white", "lightblue", "blue", "darkblue"))(100)) +
         labs(title = "Predicted Midge Distribution",
              subtitle = "Probability of midge presence",
              fill = "Probability") +
         theme_minimal()
   })

   # Risk map plot
   risk_plot <- reactive({
      req(sim_results$risk_map)

      ggplot() +
         geom_spatraster(data = sim_results$risk_map) +
         geom_sf(data = sim_results$study_area, fill = NA, color = "black") +
         geom_sf(data = sim_results$water_bodies, color ="#81C3D7", fill ="#81C3D7", alpha = 0.6) +
         # Add labels for feeder IDs
         geom_sf_label(data = sim_results$feeders,
                       aes(label = id),
                       color = "white",
                       fill = "#D95F02",
                       size = 3,
                       fontface = "bold") +
         scale_fill_gradientn(colors = colorRampPalette(c("white", "yellow", "orange", "red"))(100),
                              name = "Density") +
         theme_minimal() +
         labs(title = "Normalized Disease Risk Map",
              subtitle = "Animal Use × Midge Probability")
   })

   # Render the plot
   output$risk_plot <- renderPlot({
      risk_plot()
   })

   # Risk histogram
   output$risk_histogram <- renderPlot({
      req(sim_results$risk_map)

      # Get raster values
      risk_values <- terra::values(sim_results$risk_map)
      risk_values <- risk_values[!is.na(risk_values)]

      # Create data.frame for ggplot
      risk_df <- data.frame(risk = risk_values)

      # Create the histogram
      ggplot(risk_df, aes(x = risk)) +
         geom_histogram(aes(fill = after_stat(x)), bins = 30) +
         scale_fill_gradientn(colors = colorRampPalette(c("blue", "green", "yellow", "orange", "red"))(100)) +
         labs(title = "Distribution of Risk Values",
              x = "Risk Value",
              y = "Frequency") +
         theme_minimal() +
         theme(legend.position = "none")
   })

   # Feeder Risk Table
   feeder_risk_data <- reactive({
      req(sim_results$feeder_risk)
      sim_results$feeder_risk |>
         dplyr::select(feeder_id, risk_point, risk_mean, risk_max) |>
         mutate(risk_point = round(risk_point, 3),
                risk_mean = round(risk_mean, 3),
                risk_max = round(risk_max, 3)) |>
         rename(`Feeder ID` = feeder_id,
                `Risk Point` = risk_point,
                `Risk Average` = risk_mean,
                `Risk Max` = risk_max)
   })

   # Use the reactive data for the table output
   output$feeder_risk_table <- renderTable({
      feeder_risk_data()
   })

   # Download handlers for plots
   # Landscape
   output$download_landscape <- downloadHandler(
      filename = function() {
         paste("landscape-plot-", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
         ggsave(file, plot = landscape_plot(), width = 10, height = 8, dpi = 300, bg = "white")
      }
   )

   output$download_all_env <- downloadHandler(
      filename = function() {
         paste("environmental-plots-", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
         # Save the plot to a PNG file
         ggsave(file, plot = env_plots(), width = 10, height = 8, dpi = 300, bg = "white")
      }
   )

   output$download_movement <- downloadHandler(
      filename = function() {
         paste("animal-movement-", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
         ggsave(file, plot = movement_plot(), width = 10, height = 8, dpi = 300, , bg = "white")
      }
   )

   output$download_ud <- downloadHandler(
      filename = function() {
         paste("utilization-distribution-", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
         ggsave(file, plot = ud_plot(), width = 10, height = 8, dpi = 300, bg = white)
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

   output$download_midge <- downloadHandler(
      filename = function() {
         paste("midge-distribution-", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
         ggsave(file, plot = midge_plot(), width = 10, height = 8, dpi = 300, bg = white)
      }
   )

   output$download_risk <- downloadHandler(
      filename = function() {
         paste("risk-map-", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
         ggsave(file, plot = risk_plot(), width = 10, height = 8, dpi = 300, bg = white)
      }
   )

   output$download_feeder_risk <- downloadHandler(
      filename = function() {
         paste("feeder-risk-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
         write.csv(feeder_risk_data(), file, row.names = FALSE)
      }
   )
}

# Run the application
shinyApp(ui = ui, server = server)
