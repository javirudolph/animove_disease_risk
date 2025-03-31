library(shiny)
library(shinydashboard)
library(sf)
library(terra)
library(tidyterra)
library(ggplot2)
library(dplyr)
library(cowplot)

# Source the existing simulation functions
source("simulation_functions.R")

ui <- dashboardPage(
   dashboardHeader(title = "Disease Risk Simulator"),
   dashboardSidebar(
      sidebarMenu(
         menuItem("Simulation Parameters", tabName = "simulation", icon = icon("sliders")),
         menuItem("Results", tabName = "results", icon = icon("chart-area")),
         menuItem("Risk Analysis", tabName = "risk_analysis", icon = icon("virus"))
      )
   ),
   dashboardBody(
      tabItems(
         # Simulation Parameters Tab
         tabItem(tabName = "simulation",
                 fluidRow(
                    box(
                       title = "Animal Movement Parameters", status = "primary", solidHeader = TRUE,
                       numericInput("n_animals", "Number of Animals", value = 10, min = 1, max = 50),
                       numericInput("n_steps", "Simulation Steps", value = 1000, min = 100, max = 5000),
                       sliderInput("step_scale", "Step Scale", min = 10, max = 50, value = 20),
                       sliderInput("feeder_attraction", "Feeder Attraction Strength",
                                   min = 0, max = 1, value = 0.2, step = 0.1)
                    ),
                    box(
                       title = "Environmental Parameters", status = "primary", solidHeader = TRUE,
                       sliderInput("n_lakes", "Number of Lakes", min = 1, max = 10, value = 3),
                       sliderInput("n_feeders", "Number of Feeders", min = 1, max = 10, value = 5)
                    )
                 ),
                 actionButton("run_simulation", "Run Simulation", icon = icon("play"))
         ),

         # Results Tab
         tabItem(tabName = "results",
                 fluidRow(
                    box(
                       title = "Animal Movement Tracks", status = "info",
                       plotOutput("animal_tracks_plot")
                    ),
                    box(
                       title = "Utilization Distribution", status = "info",
                       plotOutput("utilization_plot")
                    )
                 ),
                 fluidRow(
                    box(
                       title = "Midge Distribution", status = "warning",
                       plotOutput("midge_distribution_plot")
                    ),
                    box(
                       title = "Disease Risk Map", status = "danger",
                       plotOutput("risk_map_plot")
                    )
                 )
         ),

         # Risk Analysis Tab
         tabItem(tabName = "risk_analysis",
                 fluidRow(
                    box(
                       title = "Risk Summary Statistics", status = "primary",
                       uiOutput("risk_summary")
                    ),
                    box(
                       title = "Spatial Risk Distribution", status = "primary",
                       plotOutput("risk_histogram")
                    )
                 ),
                 fluidRow(
                    box(
                       title = "Recommendations", status = "success",
                       htmlOutput("management_recommendations")
                    )
                 )
         )
      )
   )
)

server <- function(input, output, session) {
   # Reactive simulation function
   simulate_scenario <- eventReactive(input$run_simulation, {
      set.seed(508)  # Consistent seed for reproducibility

      # Create study area
      study_area <- create_study_area()

      # Create water bodies
      water_bodies <- create_water_bodies(study_area, n_lakes = input$n_lakes)

      # Create feeders
      feeders <- create_feeders(study_area, n = input$n_feeders)

      # Create environmental rasters
      env_rasters <- create_env_rasters(study_area)

      # Simulate animal movement
      animal_tracks <- simulate_animal_movement(
         study_area = study_area,
         water_bodies = water_bodies,
         feeders = feeders,
         n_animals = input$n_animals,
         n_steps = input$n_steps,
         step_scale = input$step_scale,
         feeder_attraction = input$feeder_attraction
      )

      # Create animal utilization distribution
      animal_ud <- create_animal_ud(animal_tracks, env_rasters)

      # Simulate midge data
      midge_sim <- simulate_midge_data(env_rasters)
      midge_data <- midge_sim$midge_data

      # Fit midge distribution model
      midge_sdm <- fit_midge_sdm(midge_data, env_rasters)
      midge_prediction <- midge_sdm$prediction

      # Create risk map
      risk_map <- create_risk_map(animal_ud, midge_prediction)

      list(
         study_area = study_area,
         water_bodies = water_bodies,
         feeders = feeders,
         animal_tracks = animal_tracks,
         animal_ud = animal_ud,
         midge_data = midge_data,
         midge_prediction = midge_prediction,
         risk_map = risk_map
      )
   })

   # Visualization Outputs
   output$animal_tracks_plot <- renderPlot({
      req(simulate_scenario())
      sim_results <- simulate_scenario()

      sim_results$animal_tracks |>
         mutate(x = st_coordinates(geometry)[,1],
                y = st_coordinates(geometry)[,2]) |>
         ggplot() +
         facet_wrap(~ factor(animal_id)) +
         geom_path(aes(x = x, y = y, group = factor(animal_id), color = factor(animal_id))) +
         theme_minimal() +
         geom_sf(data = sim_results$feeders) +
         geom_sf(data = sim_results$water_bodies)
   })

   output$utilization_plot <- renderPlot({
      req(simulate_scenario())
      sim_results <- simulate_scenario()

      ggplot() +
         geom_sf(data = sim_results$study_area) +
         tidyterra::geom_spatraster(data = sim_results$animal_ud) +
         scale_fill_viridis_c(name = "Utilization\nDistribution", option = "plasma") +
         theme_void() +
         labs(title = "Animal Utilization Distribution")
   })

   output$midge_distribution_plot <- renderPlot({
      req(simulate_scenario())
      sim_results <- simulate_scenario()

      ggplot() +
         geom_sf(data = sim_results$study_area) +
         tidyterra::geom_spatraster(data = sim_results$midge_prediction) +
         scale_fill_viridis_c(name = "Midge\nProbability", option = "plasma") +
         theme_void() +
         labs(title = "Midge Species Distribution Model")
   })

   output$risk_map_plot <- renderPlot({
      req(simulate_scenario())
      sim_results <- simulate_scenario()

      ggplot() +
         geom_sf(data = sim_results$study_area) +
         tidyterra::geom_spatraster(data = sim_results$risk_map) +
         scale_fill_viridis_c(name = "Disease Risk", option = "magma") +
         theme_void() +
         labs(title = "Normalized Disease Risk Map")
   })

   # Risk Analysis Outputs
   output$risk_summary <- renderUI({
      req(simulate_scenario())
      sim_results <- simulate_scenario()

      risk_values <- values(sim_results$risk_map)

      HTML(paste(
         "<b>Risk Summary Statistics:</b><br>",
         "Mean Risk: ", round(mean(risk_values, na.rm = TRUE), 4), "<br>",
         "Median Risk: ", round(median(risk_values, na.rm = TRUE), 4), "<br>",
         "Maximum Risk: ", round(max(risk_values, na.rm = TRUE), 4), "<br>",
         "Minimum Risk: ", round(min(risk_values, na.rm = TRUE), 4), "<br>",
         "Risk Standard Deviation: ", round(sd(risk_values, na.rm = TRUE), 4)
      ))
   })

   output$risk_histogram <- renderPlot({
      req(simulate_scenario())
      sim_results <- simulate_scenario()

      risk_values <- values(sim_results$risk_map)

      ggplot(data.frame(risk = risk_values), aes(x = risk)) +
         geom_histogram(bins = 30, fill = "steelblue", color = "black") +
         theme_minimal() +
         labs(title = "Disease Risk Distribution",
              x = "Risk Value", y = "Frequency")
   })

   output$management_recommendations <- renderUI({
      req(simulate_scenario())
      sim_results <- simulate_scenario()

      risk_values <- values(sim_results$risk_map)
      mean_risk <- mean(risk_values, na.rm = TRUE)

      recommendations <- if(mean_risk > 0.7) {
         paste(
            "<b>High Risk Recommendations:</b><br>",
            "• Implement immediate vector control measures<br>",
            "• Consider reducing animal density in high-risk areas<br>",
            "• Increase surveillance and monitoring<br>",
            "• Deploy additional protective barriers or treatments"
         )
      } else if(mean_risk > 0.4) {
         paste(
            "<b>Moderate Risk Recommendations:</b><br>",
            "• Enhanced monitoring of animal movement patterns<br>",
            "• Selective vector control in targeted areas<br>",
            "• Implement preventive vaccination strategies<br>",
            "• Regular environmental management"
         )
      } else {
         paste(
            "<b>Low Risk Recommendations:</b><br>",
            "• Continue routine surveillance<br>",
            "• Maintain current management practices<br>",
            "• Periodic risk reassessment<br>",
            "• Environmental health monitoring"
         )
      }

      HTML(recommendations)
   })
}

# Run the application
shinyApp(ui, server)
