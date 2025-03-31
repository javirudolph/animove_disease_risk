library(shiny)
library(shinydashboard)
library(sf)
library(terra)
library(tidyterra)
library(ggplot2)
library(dplyr)
library(cowplot)

# Preprocessing function to generate multiple scenarios
generate_scenarios <- function(n_scenarios = 5) {
   scenarios <- list()
   set.seed(508)

   for(i in 1:n_scenarios) {
      # Randomize parameters
      n_animals <- sample(5:20, 1)
      n_steps <- sample(500:2000, 1)
      step_scale <- runif(1, 10, 50)
      feeder_attraction <- runif(1, 0, 1)
      n_lakes <- sample(1:5, 1)
      n_feeders <- sample(3:8, 1)

      # Create study area
      study_area <- create_study_area()
      water_bodies <- create_water_bodies(study_area, n_lakes = n_lakes)
      feeders <- create_feeders(study_area, n = n_feeders)
      env_rasters <- create_env_rasters(study_area)

      # Simulate scenarios
      animal_tracks <- simulate_animal_movement(
         study_area = study_area,
         water_bodies = water_bodies,
         feeders = feeders,
         n_animals = n_animals,
         n_steps = n_steps,
         step_scale = step_scale,
         feeder_attraction = feeder_attraction
      )

      animal_ud <- create_animal_ud(animal_tracks, env_rasters)
      midge_sim <- simulate_midge_data(env_rasters)
      midge_data <- midge_sim$midge_data
      midge_sdm <- fit_midge_sdm(midge_data, env_rasters)
      midge_prediction <- midge_sdm$prediction
      risk_map <- create_risk_map(animal_ud, midge_prediction)

      scenarios[[i]] <- list(
         params = list(
            n_animals = n_animals,
            n_steps = n_steps,
            step_scale = step_scale,
            feeder_attraction = feeder_attraction,
            n_lakes = n_lakes,
            n_feeders = n_feeders
         ),
         study_area = study_area,
         water_bodies = water_bodies,
         feeders = feeders,
         animal_tracks = animal_tracks,
         animal_ud = animal_ud,
         midge_data = midge_data,
         midge_prediction = midge_prediction,
         risk_map = risk_map
      )
   }

   return(scenarios)
}

# LLM-like recommendation generator
generate_recommendations <- function(mean_risk, scenario_details) {
   # Use a predefined set of recommendations with some randomness
   risk_categories <- list(
      high = list(
         general = c(
            "Implement immediate vector control measures",
            "Consider reducing animal density in high-risk areas",
            "Increase surveillance and monitoring",
            "Deploy additional protective barriers or treatments"
         ),
         specific = c(
            "Conduct targeted culling in high-risk zones",
            "Implement intensive midge control programs",
            "Establish temporary movement restrictions for animals"
         )
      ),
      moderate = list(
         general = c(
            "Enhanced monitoring of animal movement patterns",
            "Selective vector control in targeted areas",
            "Implement preventive vaccination strategies",
            "Regular environmental management"
         ),
         specific = c(
            "Develop localized intervention plans",
            "Monitor water bodies and vegetation for midge breeding sites",
            "Implement rotational grazing strategies"
         )
      ),
      low = list(
         general = c(
            "Continue routine surveillance",
            "Maintain current management practices",
            "Periodic risk reassessment",
            "Environmental health monitoring"
         ),
         specific = c(
            "Maintain current habitat management",
            "Conduct periodic population health checks",
            "Update risk models annually"
         )
      )
   )

   # Categorize risk
   if(mean_risk > 0.7) {
      risk_level <- "high"
   } else if(mean_risk > 0.4) {
      risk_level <- "moderate"
   } else {
      risk_level <- "low"
   }

   # Generate recommendation
   recommendations <- list(
      title = paste("Management Strategy for", risk_level, "Risk Scenario"),
      risk_level = risk_level,
      general_recommendations = sample(risk_categories[[risk_level]]$general, 3),
      specific_recommendations = sample(risk_categories[[risk_level]]$specific, 2),
      scenario_details = paste(
         "Scenario Details:",
         "Animals:", scenario_details$n_animals,
         "Steps:", scenario_details$n_steps,
         "Feeder Attraction:", round(scenario_details$feeder_attraction, 2),
         "Lakes:", scenario_details$n_lakes
      )
   )

   return(recommendations)
}

# Pre-generate scenarios
pre_generated_scenarios <- generate_scenarios()

ui <- dashboardPage(
   dashboardHeader(title = "Disease Risk Simulator"),
   dashboardSidebar(
      sidebarMenu(
         menuItem("Scenario Selection", tabName = "scenario_select", icon = icon("list")),
         menuItem("Results", tabName = "results", icon = icon("chart-area")),
         menuItem("Risk Analysis", tabName = "risk_analysis", icon = icon("virus"))
      )
   ),
   dashboardBody(
      tabItems(
         # Scenario Selection Tab
         tabItem(tabName = "scenario_select",
                 fluidRow(
                    box(
                       title = "Scenario Parameters", width = 12, status = "primary",
                       selectInput("selected_scenario",
                                   "Choose a Pre-generated Scenario:",
                                   choices = setNames(
                                      1:length(pre_generated_scenarios),
                                      sapply(pre_generated_scenarios, function(scen) {
                                         paste("Scenario:",
                                               "Animals =", scen$params$n_animals,
                                               "| Steps =", scen$params$n_steps,
                                               "| Lakes =", scen$params$n_lakes)
                                      })
                                   )
                       )
                    )
                 )
         ),

         # Results Tab (similar to previous implementation)
         tabItem(tabName = "results",
                 fluidRow(
                    box(title = "Animal Movement Tracks", plotOutput("animal_tracks_plot")),
                    box(title = "Utilization Distribution", plotOutput("utilization_plot"))
                 ),
                 fluidRow(
                    box(title = "Midge Distribution", plotOutput("midge_distribution_plot")),
                    box(title = "Disease Risk Map", plotOutput("risk_map_plot"))
                 )
         ),

         # Risk Analysis Tab
         tabItem(tabName = "risk_analysis",
                 fluidRow(
                    box(title = "Risk Summary Statistics", uiOutput("risk_summary")),
                    box(title = "Spatial Risk Distribution", plotOutput("risk_histogram"))
                 ),
                 fluidRow(
                    box(title = "AI-Generated Recommendations",
                        width = 12,
                        status = "success",
                        htmlOutput("management_recommendations"))
                 )
         )
      )
   )
)

server <- function(input, output, session) {
   # Reactive scenario selection
   selected_scenario <- reactive({
      req(input$selected_scenario)
      pre_generated_scenarios[[as.numeric(input$selected_scenario)]]
   })

   # Visualization Outputs (mostly unchanged from previous version)
   output$animal_tracks_plot <- renderPlot({
      req(selected_scenario())
      sim_results <- selected_scenario()

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

   # [Other plot outputs remain the same as in previous version]

   # Risk Analysis Outputs
   output$risk_summary <- renderUI({
      req(selected_scenario())
      sim_results <- selected_scenario()

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

   output$management_recommendations <- renderUI({
      req(selected_scenario())
      sim_results <- selected_scenario()

      risk_values <- values(sim_results$risk_map)
      mean_risk <- mean(risk_values, na.rm = TRUE)

      # Generate recommendations using our custom function
      recommendations <- generate_recommendations(
         mean_risk,
         sim_results$params
      )

      HTML(paste(
         "<h4>", recommendations$title, "</h4>",
         "<b>Risk Level:</b> ", recommendations$risk_level, "<br><br>",
         "<b>General Recommendations:</b><br>",
         paste("• ", recommendations$general_recommendations, collapse = "<br>"),
         "<br><br>",
         "<b>Specific Recommendations:</b><br>",
         paste("• ", recommendations$specific_recommendations, collapse = "<br>"),
         "<br><br>",
         recommendations$scenario_details
      ))
   })
}

# Run the application
shinyApp(ui, server)
