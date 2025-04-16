# Hemorrhagic Disease Risk Assessment Tool
# A Shiny application for visualizing and exploring hemorrhagic disease risk in white-tailed deer

library(shiny)
library(shinydashboard)
library(leaflet)
library(leaflet.extras)  # Added for heatmap functionality
library(dplyr)
library(ggplot2)
library(plotly)

# UI for the application
ui <- dashboardPage(
   dashboardHeader(title = "Hemorrhagic Disease Risk Tool"),

   dashboardSidebar(
      sidebarMenu(id = "sidebar",  # Give the sidebarMenu an ID
                  menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
                  menuItem("About", tabName = "about", icon = icon("info-circle"))
      ),

      # Input controls for the dashboard - FIXED CONDITION
      conditionalPanel(
         condition = "input.sidebar == 'dashboard'",  # Now correctly refers to the sidebar input
         h4("Simulation Parameters"),
         selectInput("season", "Season:",
                     choices = c("Spring", "Summer", "Fall", "Winter"),
                     selected = "Summer"),
         sliderInput("time_of_day", "Time of Day:",
                     min = 0, max = 23, value = 18, step = 1),
         sliderInput("vector_activity", "Vector Activity Level:",
                     min = 0, max = 100, value = 75, step = 5),
         selectInput("management_strategy", "Management Strategy:",
                     choices = c("None", "Targeted Habitat Modification",
                                 "Vector Control", "Deer Exclusion"),
                     selected = "None"),
         hr(),
         actionButton("run_simulation", "Run Simulation",
                      icon = icon("play"),
                      style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
      )
   ),

   dashboardBody(
      tabItems(
         # Dashboard tab
         tabItem(
            tabName = "dashboard",
            fluidRow(
               box(
                  title = "Risk Map", width = 8, solidHeader = TRUE, status = "primary",
                  leafletOutput("risk_map", height = 500)
               ),
               column(
                  width = 4,
                  box(
                     title = "Risk Assessment Summary", width = NULL, solidHeader = TRUE, status = "info",
                     valueBoxOutput("risk_score_box", width = NULL),
                     plotlyOutput("risk_distribution", height = 200)
                  ),
                  box(
                     title = "Temporal Trends", width = NULL, solidHeader = TRUE, status = "warning",
                     plotlyOutput("temporal_trends", height = 200)
                  )
               )
            ),
            fluidRow(
               box(
                  title = "Host-Vector Overlap Analysis", width = 12, solidHeader = TRUE, status = "success",
                  fluidRow(
                     column(
                        width = 4,
                        h4("Deer Utilization"),
                        plotOutput("deer_heatmap", height = 200)
                     ),
                     column(
                        width = 4,
                        h4("Vector Distribution"),
                        plotOutput("vector_heatmap", height = 200)
                     ),
                     column(
                        width = 4,
                        h4("Management Impact"),
                        plotOutput("management_impact", height = 200)
                     )
                  )
               )
            )
         ),

         # About tab
         tabItem(
            tabName = "about",
            box(
               title = "About This Tool", width = 12, solidHeader = TRUE, status = "primary",
               h3("Integrated Host-Vector Modeling for Hemorrhagic Disease Management"),
               p("This application demonstrates a decision support tool for managing hemorrhagic disease in white-tailed deer populations. The tool combines animal movement simulation with species distribution modeling of Culicoides biting midges to provide a framework for exploring management scenarios."),

               h4("Data Sources"),
               p("This tool is powered by:"),
               tags$ul(
                  tags$li("GPS collar data from white-tailed deer (Odocoileus virginianus)"),
                  tags$li("Species distribution models for Culicoides biting midges"),
                  tags$li("Environmental and landscape data layers")
               ),

               h4("Methodology"),
               p("Our approach integrates:"),
               tags$ul(
                  tags$li("Step selection analysis to extract movement patterns from GPS collar data"),
                  tags$li("Species distribution modeling for vector populations"),
                  tags$li("Spatial overlay analysis to identify high-risk areas"),
                  tags$li("Scenario simulation to test management interventions")
               ),

               h4("How to Use This Tool"),
               p("Use the controls in the sidebar to:"),
               tags$ul(
                  tags$li("Select different seasons and times of day to see how risk patterns change"),
                  tags$li("Adjust vector activity levels to simulate different outbreak conditions"),
                  tags$li("Test different management strategies to see their impact on disease risk"),
                  tags$li("Click 'Run Simulation' to update the analysis based on your selections")
               ),

               hr(),
               p("This is a demonstration version with simulated data. A full version would incorporate real field data from your specific study area.")
            )
         )
      )
   )
)

# Server logic
server <- function(input, output, session) {

   # Create some simulated data
   # In a real application, this would be loaded from actual data files

   # Simulated farm boundary
   farm_boundary <- data.frame(
      lng = c(-82.5, -82.4, -82.3, -82.35, -82.5),
      lat = c(35.6, 35.65, 35.55, 35.5, 35.6)
   )

   # Generate reactive values that will change based on user inputs
   risk_data <- reactiveVal()

   # Generate simulated risk data
   generate_risk_data <- function() {
      # Create a grid over our area
      lng_seq <- seq(-82.55, -82.25, by = 0.01)
      lat_seq <- seq(35.45, 35.7, by = 0.01)
      grid <- expand.grid(lng = lng_seq, lat = lat_seq)

      # Base deer utilization (higher in certain areas)
      grid$deer_util <- with(grid,
                             exp(-(lng + 82.4)^2 * 100) * exp(-(lat - 35.6)^2 * 100) * 10 +
                                exp(-(lng + 82.35)^2 * 100) * exp(-(lat - 35.53)^2 * 100) * 8)

      # Adjust for season
      season_factor <- switch(input$season,
                              "Spring" = 0.7,
                              "Summer" = 1.0,
                              "Fall" = 0.9,
                              "Winter" = 0.5)

      grid$deer_util <- grid$deer_util * season_factor

      # Adjust for time of day (higher at dawn/dusk)
      time_factor <- dnorm(input$time_of_day, mean = 7, sd = 2) * 2 +
         dnorm(input$time_of_day, mean = 19, sd = 2) * 2 + 0.2

      grid$deer_util <- grid$deer_util * time_factor

      # Vector distribution (different spatial pattern)
      grid$vector_dist <- with(grid,
                               exp(-(lng + 82.45)^2 * 80) * exp(-(lat - 35.58)^2 * 80) * 10 +
                                  exp(-(lng + 82.32)^2 * 80) * exp(-(lat - 35.55)^2 * 80) * 7)

      # Adjust for vector activity level
      grid$vector_dist <- grid$vector_dist * (input$vector_activity / 75)

      # Calculate baseline risk as product of deer utilization and vector distribution
      grid$risk <- grid$deer_util * grid$vector_dist

      # Apply management impact if selected
      if (input$management_strategy != "None") {
         # Define impact areas for different strategies
         if (input$management_strategy == "Targeted Habitat Modification") {
            # Reduce risk in specific managed areas
            habitat_mod <- with(grid,
                                exp(-(lng + 82.4)^2 * 150) * exp(-(lat - 35.57)^2 * 150) > 0.3)
            grid$risk[habitat_mod] <- grid$risk[habitat_mod] * 0.4
         }

         if (input$management_strategy == "Vector Control") {
            # Reduce vector component of risk in treatment areas
            vector_control <- with(grid,
                                   exp(-(lng + 82.38)^2 * 120) * exp(-(lat - 35.56)^2 * 120) > 0.2)
            grid$risk[vector_control] <- grid$risk[vector_control] * 0.3
         }

         if (input$management_strategy == "Deer Exclusion") {
            # Eliminate deer component in exclusion zones
            exclusion <- with(grid,
                              exp(-(lng + 82.42)^2 * 200) * exp(-(lat - 35.59)^2 * 200) > 0.4)
            grid$risk[exclusion] <- grid$risk[exclusion] * 0.1
         }
      }

      # Normalize for visualization
      grid$risk_normalized <- grid$risk / max(grid$risk) * 100

      return(grid)
   }

   # Initialize data on app start - CRITICAL FIX
   observe({
      risk_data(generate_risk_data())
   }, priority = 1000)

   # Update the risk data when simulation is run
   observeEvent(input$run_simulation, {
      risk_data(generate_risk_data())
   })

   # Create the risk map
   output$risk_map <- renderLeaflet({
      # Initialize the leaflet map
      leaflet() %>%
         addProviderTiles(providers$Esri.WorldImagery) %>%
         setView(lng = -82.4, lat = 35.57, zoom = 13) %>%
         addPolygons(
            data = farm_boundary,
            lng = ~lng,
            lat = ~lat,
            fillColor = "transparent",
            color = "white",
            weight = 2,
            opacity = 1,
            layerId = "farm"
         )
   })

   # Update the risk map when data changes
   observe({
      data <- risk_data()
      if (!is.null(data)) {
         # Sample the data to reduce points for better performance
         sampled_data <- data[sample(nrow(data), min(2000, nrow(data))), ]

         pal <- colorNumeric(
            palette = "YlOrRd",
            domain = sampled_data$risk_normalized
         )

         leafletProxy("risk_map") %>%
            clearHeatmap() %>%
            addHeatmap(
               data = sampled_data,
               lng = ~lng,
               lat = ~lat,
               intensity = ~risk_normalized,
               blur = 15,
               max = 100,
               radius = 10
            ) %>%
            clearControls() %>%
            addLegend(
               "bottomright",
               pal = pal,
               values = sampled_data$risk_normalized,
               title = "Disease Risk",
               opacity = 1
            )
      }
   })

   # Risk score output
   output$risk_score_box <- renderValueBox({
      data <- risk_data()
      if (!is.null(data)) {
         avg_risk <- mean(data$risk_normalized)
         risk_level <- cut(
            avg_risk,
            breaks = c(0, 20, 40, 60, 80, 100),
            labels = c("Very Low", "Low", "Moderate", "High", "Very High")
         )

         risk_color <- switch(
            as.character(risk_level),
            "Very Low" = "green",
            "Low" = "light-blue",
            "Moderate" = "yellow",
            "High" = "orange",
            "Very High" = "red"
         )

         valueBox(
            paste0(round(avg_risk), "/100"),
            paste("Overall Risk Level:", risk_level),
            icon = icon("biohazard"),
            color = risk_color
         )
      }
   })

   # Risk distribution plot
   output$risk_distribution <- renderPlotly({
      data <- risk_data()
      if (!is.null(data)) {
         # Sample data for better performance
         sampled_data <- data[sample(nrow(data), min(5000, nrow(data))), ]

         p <- ggplot(sampled_data, aes(x = risk_normalized)) +
            geom_histogram(binwidth = 5, fill = "#3c8dbc", color = "#2c6b8e") +
            theme_minimal() +
            labs(x = "Risk Score", y = "Count") +
            theme(axis.title = element_text(size = 10),
                  axis.text = element_text(size = 8))

         ggplotly(p) %>%
            layout(margin = list(l = 50, r = 50, b = 50, t = 10, pad = 4))
      }
   })

   # Temporal trends plot
   output$temporal_trends <- renderPlotly({
      # Create simulated daily risk pattern
      hours <- 0:23

      # Different patterns for different seasons
      if (input$season == "Summer") {
         risk_pattern <- dnorm(hours, mean = 7, sd = 2) * 70 +
            dnorm(hours, mean = 19, sd = 2) * 100 + 20
      } else if (input$season == "Fall") {
         risk_pattern <- dnorm(hours, mean = 8, sd = 2) * 60 +
            dnorm(hours, mean = 18, sd = 2) * 90 + 15
      } else if (input$season == "Spring") {
         risk_pattern <- dnorm(hours, mean = 6, sd = 2) * 50 +
            dnorm(hours, mean = 18.5, sd = 2) * 80 + 10
      } else {
         risk_pattern <- dnorm(hours, mean = 9, sd = 2) * 30 +
            dnorm(hours, mean = 17, sd = 2) * 50 + 5
      }

      # Adjust for management
      if (input$management_strategy != "None") {
         reduction_factor <- switch(
            input$management_strategy,
            "Targeted Habitat Modification" = 0.7,
            "Vector Control" = 0.5,
            "Deer Exclusion" = 0.3,
            1
         )

         risk_pattern <- risk_pattern * reduction_factor
      }

      temporal_data <- data.frame(
         hour = hours,
         risk = risk_pattern
      )

      # Highlight the current selected time
      highlight_data <- data.frame(
         hour = input$time_of_day,
         risk = risk_pattern[input$time_of_day + 1]
      )

      p <- ggplot() +
         geom_line(data = temporal_data, aes(x = hour, y = risk), size = 1, color = "#f39c12") +
         geom_point(data = highlight_data, aes(x = hour, y = risk), size = 4, color = "red") +
         theme_minimal() +
         labs(x = "Hour of Day", y = "Risk Level") +
         scale_x_continuous(breaks = seq(0, 24, by = 4)) +
         theme(axis.title = element_text(size = 10),
               axis.text = element_text(size = 8))

      ggplotly(p) %>%
         layout(margin = list(l = 50, r = 50, b = 50, t = 10, pad = 4))
   })

   # Deer utilization heatmap
   output$deer_heatmap <- renderPlot({
      data <- risk_data()
      if (!is.null(data)) {
         # Sample data for better performance
         sampled_data <- data[sample(nrow(data), min(2000, nrow(data))), ]

         ggplot(sampled_data, aes(x = lng, y = lat, fill = deer_util)) +
            geom_tile() +
            scale_fill_gradient(low = "white", high = "blue") +
            theme_void() +
            theme(legend.position = "none") +
            labs(title = "Deer Activity")
      }
   })

   # Vector distribution heatmap
   output$vector_heatmap <- renderPlot({
      data <- risk_data()
      if (!is.null(data)) {
         # Sample data for better performance
         sampled_data <- data[sample(nrow(data), min(2000, nrow(data))), ]

         ggplot(sampled_data, aes(x = lng, y = lat, fill = vector_dist)) +
            geom_tile() +
            scale_fill_gradient(low = "white", high = "red") +
            theme_void() +
            theme(legend.position = "none") +
            labs(title = "Vector Density")
      }
   })

   # Management impact visualization
   output$management_impact <- renderPlot({
      data <- risk_data()
      if (!is.null(data)) {
         # Sample data for better performance
         sampled_data <- data[sample(nrow(data), min(2000, nrow(data))), ]

         # Generate data without management for comparison
         if (input$management_strategy != "None") {
            # We'll need to temporarily create data without management
            # We can't modify the input directly as it would trigger reactivity loops
            baseline_data <- sampled_data

            # Define impact areas for different strategies based on current data
            if (input$management_strategy == "Targeted Habitat Modification") {
               habitat_mod <- with(sampled_data,
                                   exp(-(lng + 82.4)^2 * 150) * exp(-(lat - 35.57)^2 * 150) > 0.3)
               baseline_data$risk[habitat_mod] <- baseline_data$risk[habitat_mod] / 0.4
            }

            if (input$management_strategy == "Vector Control") {
               vector_control <- with(sampled_data,
                                      exp(-(lng + 82.38)^2 * 120) * exp(-(lat - 35.56)^2 * 120) > 0.2)
               baseline_data$risk[vector_control] <- baseline_data$risk[vector_control] / 0.3
            }

            if (input$management_strategy == "Deer Exclusion") {
               exclusion <- with(sampled_data,
                                 exp(-(lng + 82.42)^2 * 200) * exp(-(lat - 35.59)^2 * 200) > 0.4)
               baseline_data$risk[exclusion] <- baseline_data$risk[exclusion] / 0.1
            }

            # Calculate risk reduction
            reduction <- (baseline_data$risk - sampled_data$risk) / baseline_data$risk * 100
            reduction[reduction < 0] <- 0  # Ensure no negative values
            reduction[is.na(reduction)] <- 0  # Replace any NAs with 0

            reduction_df <- data.frame(
               lng = sampled_data$lng,
               lat = sampled_data$lat,
               reduction = reduction
            )

            ggplot(reduction_df, aes(x = lng, y = lat, fill = reduction)) +
               geom_tile() +
               scale_fill_gradient(low = "white", high = "green") +
               theme_void() +
               theme(legend.position = "none") +
               labs(title = "Risk Reduction (%)")
         } else {
            # If no management, show blank plot with message
            ggplot() +
               theme_void() +
               annotate("text", x = 0, y = 0, label = "No management\nselected", size = 5) +
               labs(title = "Risk Reduction (%)")
         }
      }
   })
}

# Run the application
shinyApp(ui = ui, server = server)
