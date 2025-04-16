# Rainbow Pastel Disease Risk Simulation Shiny App
# Featuring a soft rainbow pastel color scheme

# Load required packages
library(shiny)
library(sf)
library(terra)
library(tidyterra)
library(ggplot2)
library(dplyr)
library(spsurvey)
library(shinyjs)
library(colourpicker) # for pastel colors
library(RColorBrewer)
library(viridis)

# Source the simulation functions
source("R/simulation_functions.R")

# Rainbow pastel colors
pastel_colors <- c(
   "#FFD6D6", # soft pink
   "#FFECDB", # pastel peach
   "#FFFBD6", # pastel yellow
   "#DEFFDB", # pastel mint
   "#D6F6FF", # pastel sky blue
   "#D6DBFF", # pastel lavender
   "#F6D6FF"  # pastel lilac
)

# Custom theme CSS
app_css <- "
  .sidebar {
    background-color: #FFF1F9;
    padding: 20px;
    border-radius: 15px;
    box-shadow: 0 0 10px rgba(0,0,0,0.1);
  }

  .title-panel {
    background: linear-gradient(to right, #FFD6D6, #FFECDB, #FFFBD6, #DEFFDB, #D6F6FF, #D6DBFF, #F6D6FF);
    padding: 15px;
    margin-bottom: 20px;
    border-radius: 10px;
    text-align: center;
    color: #5a5a5a;
    box-shadow: 0 0 5px rgba(0,0,0,0.1);
  }

  .btn-pastel {
    background: linear-gradient(to right, #FFD6D6, #F6D6FF);
    border: none;
    color: #5a5a5a;
    font-weight: bold;
  }

  .btn-pastel:hover {
    background: linear-gradient(to right, #F6D6FF, #FFD6D6);
    color: #3a3a3a;
  }

  .nav-tabs > li > a {
    background-color: #F6D6FF;
    color: #5a5a5a;
    border-radius: 10px 10px 0 0;
    margin-right: 5px;
  }

  .nav-tabs > li.active > a {
    background-color: #D6F6FF;
    color: #333;
    font-weight: bold;
  }

  .help-block {
    color: #8c8c8c;
    font-style: italic;
  }

  hr {
    border-color: #D6DBFF;
  }

  .slider-container .irs-bar {
    background: linear-gradient(to right, #FFD6D6, #F6D6FF);
  }

  .slider-container .irs-handle {
    background: #D6F6FF;
    border-color: #D6DBFF;
  }
"

# UI
ui <- fluidPage(
   useShinyjs(),
   tags$head(
      tags$style(app_css)
   ),

   # Title panel with rainbow gradient
   div(class = "title-panel",
       h1("Disease Risk Simulation", style = "font-family: 'Comic Sans MS', cursive, sans-serif;"),
       p("Explore how environmental factors affect disease risk", style = "font-size: 16px;")
   ),

   # Main layout
   sidebarLayout(
      sidebarPanel(
         class = "sidebar",
         # Key parameters to modify
         div(class = "slider-container",
             sliderInput("water_bodies", "Number of Water Bodies:",
                         min = 1, max = 8, value = 4, step = 1)
         ),

         div(class = "slider-container",
             sliderInput("n_feeders", "Number of Feeders:",
                         min = 5, max = 15, value = 10, step = 1)
         ),

         div(class = "slider-container",
             sliderInput("n_animals", "Number of Animals:",
                         min = 1, max = 10, value = 5, step = 1)
         ),

         div(class = "slider-container",
             sliderInput("n_steps", "Number of Movement Steps:",
                         min = 50, max = 300, value = 200, step = 50)
         ),

         actionButton("run_simulation", "Run Simulation",
                      class = "btn-pastel",
                      style = "width: 100%; margin-top: 20px; height: 50px; font-size: 18px;"),

         hr(),
         helpText("Note: Other parameters are set to default values")
      ),

      mainPanel(
         tabsetPanel(
            tabPanel("Landscape",
                     br(),
                     plotOutput("landscape_plot", height = "500px"),
                     br(),
                     div(style = "background-color: #FFFBD6; padding: 15px; border-radius: 10px;",
                         h4("Landscape Features"),
                         p("This tab shows the simulated landscape with water bodies (blue) and feeders (pink)."),
                         p("Changes to the number of water bodies and feeders affect animal movement patterns.")
                     )
            ),

            tabPanel("Animal Movement",
                     br(),
                     plotOutput("movement_plot", height = "500px"),
                     br(),
                     div(style = "background-color: #DEFFDB; padding: 15px; border-radius: 10px;",
                         h4("Animal Movement Patterns"),
                         p("This visualization shows how animals move across the landscape."),
                         p("Each color represents a different animal's movement track.")
                     )
            ),

            tabPanel("Risk Map",
                     br(),
                     plotOutput("risk_plot", height = "500px"),
                     br(),
                     div(style = "background-color: #D6F6FF; padding: 15px; border-radius: 10px;",
                         h4("Disease Risk Distribution"),
                         p("The final risk map combines animal utilization and midge distribution."),
                         p("Red areas indicate higher disease risk, while blue areas have lower risk.")
                     )
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
      midge_prediction = NULL,
      risk_map = NULL
   )

   # Run simulation when button is clicked
   observeEvent(input$run_simulation, {

      # Disable run button during simulation
      disable("run_simulation")

      # Show progress notification
      withProgress(message = 'Creating your rainbow world...', value = 0, {

         # Step 1: Create study area (fixed parameters)
         incProgress(0.1, detail = "Painting the landscape canvas")
         sim_results$study_area <- create_study_area()

         # Step 2: Create water bodies
         incProgress(0.2, detail = "Adding sparkling water bodies")
         sim_results$water_bodies <- create_water_bodies(sim_results$study_area, n = input$water_bodies)

         # Step 3: Create feeders
         incProgress(0.3, detail = "Placing colorful feeders")
         sim_results$feeders <- create_feeders(
            sim_results$study_area,
            n = input$n_feeders,
            water_bodies = sim_results$water_bodies,
            min_dist_to_water = 100,  # Default
            min_dist_between = 100    # Default
         )

         # Step 4: Create environmental rasters
         incProgress(0.4, detail = "Creating pastel environments")
         sim_results$env_rasters <- create_env_rasters(sim_results$study_area, sim_results$water_bodies)

         # Step 5: Simulate animal movement
         incProgress(0.5, detail = "Animating fluffy creatures")
         sim_results$animal_tracks <- simulate_animal_movement(
            study_area = sim_results$study_area,
            water_bodies = sim_results$water_bodies,
            feeders = sim_results$feeders,
            n_animals = input$n_animals,
            n_steps = input$n_steps
         )

         # Step 6: Create animal utilization distribution
         incProgress(0.6, detail = "Mapping animal rainbow trails")
         sim_results$animal_ud <- create_animal_ud(
            sim_results$animal_tracks,
            sim_results$study_area,
            resolution = 10,     # Default
            smoothing_factor = 9 # Default
         )

         # Step 7: Simulate midge data
         incProgress(0.7, detail = "Simulating sparkly midges")
         midge_sim <- simulate_midge_data(sim_results$env_rasters)
         sim_results$midge_data <- midge_sim$midge_data

         # Step 8: Get midge samples (using GRTS with default parameters)
         incProgress(0.8, detail = "Collecting midge data points")
         grts_sample <- spsurvey::grts(sim_results$midge_data, 100)
         sim_results$midge_sample <- grts_sample$sites_base %>%
            dplyr::select(presence, elevation, water_dist, veg_index, temperature, geometry)

         # Step 9: Fit midge distribution model
         incProgress(0.9, detail = "Creating magical midge model")
         midge_sdm <- fit_midge_sdm(sim_results$midge_sample, sim_results$env_rasters)
         sim_results$midge_prediction <- midge_sdm$prediction

         # Step 10: Create risk map
         incProgress(1.0, detail = "Finalizing rainbow risk map")
         sim_results$risk_map <- create_risk_map(sim_results$animal_ud, sim_results$midge_prediction)

      })

      # Re-enable run button
      enable("run_simulation")
   })

   # Landscape plot with pastel colors
   output$landscape_plot <- renderPlot({
      req(sim_results$study_area, sim_results$water_bodies, sim_results$feeders)

      ggplot() +
         geom_sf(data = sim_results$study_area, fill = "#FFFBF0", color = "#8C8C8C", size = 1) +
         geom_sf(data = sim_results$water_bodies, fill = "#D6F6FF", color = "#9ED2E6", size = 0.8) +
         geom_sf(data = sim_results$feeders, fill = "#FFD6D6", color = "#E6B9B9", size = 4, shape = 21) +
         labs(title = "Magical Landscape Simulation",
              subtitle = paste(input$water_bodies, "sparkling water bodies and", input$n_feeders, "colorful feeders")) +
         theme_minimal() +
         theme(
            plot.background = element_rect(fill = "#FFFDF5", color = NA),
            panel.background = element_rect(fill = "#FFFDF5", color = NA),
            plot.title = element_text(color = "#5a5a5a", size = 18, face = "bold"),
            plot.subtitle = element_text(color = "#8c8c8c", size = 12),
            panel.grid.major = element_line(color = "#F0F0F0"),
            panel.grid.minor = element_line(color = "#F8F8F8")
         )
   })

   # Animal movement plot with pastel color scheme
   output$movement_plot <- renderPlot({
      req(sim_results$animal_tracks)

      # Custom pastel palette for animal tracks
      animal_colors <- c("#FFB6C1", "#FFD700", "#98FB98", "#87CEFA", "#DDA0DD",
                         "#FFA07A", "#FFFACD", "#AFEEEE", "#D8BFD8", "#FFDAB9")

      ggplot() +
         geom_sf(data = sim_results$study_area, fill = "#FFFBF0", color = "#8C8C8C", size = 1) +
         geom_sf(data = sim_results$water_bodies, fill = "#D6F6FF", color = "#9ED2E6", size = 0.8) +
         geom_sf(data = sim_results$feeders, fill = "#FFD6D6", color = "#E6B9B9", size = 4, shape = 21) +
         geom_sf(data = sim_results$animal_tracks, aes(color = factor(animal_id)), size = 1.2, alpha = 0.7) +
         scale_color_manual(values = animal_colors) +
         labs(title = "Colorful Animal Movement Paths",
              subtitle = paste(input$n_animals, "magical creatures with", input$n_steps, "movement steps"),
              color = "Animal ID") +
         theme_minimal() +
         theme(
            plot.background = element_rect(fill = "#FFFDF5", color = NA),
            panel.background = element_rect(fill = "#FFFDF5", color = NA),
            plot.title = element_text(color = "#5a5a5a", size = 18, face = "bold"),
            plot.subtitle = element_text(color = "#8c8c8c", size = 12),
            panel.grid.major = element_line(color = "#F0F0F0"),
            panel.grid.minor = element_line(color = "#F8F8F8"),
            legend.background = element_rect(fill = "#FFFBF0"),
            legend.key = element_rect(fill = "#FFFBF0")
         )
   })

   # Risk map plot with a softer rainbow palette
   output$risk_plot <- renderPlot({
      req(sim_results$risk_map)

      # Custom rainbow pastel palette for risk map
      rainbow_pastel <- colorRampPalette(c("#FFD6D6", "#FFECDB", "#FFFBD6",
                                           "#DEFFDB", "#D6F6FF", "#D6DBFF", "#F6D6FF"))(100)

      # Plot with custom title and colors
      par(bg = "#FFFDF5", mar = c(5, 4, 4, 6) + 0.1)
      plot(sim_results$risk_map,
           main = "Rainbow Disease Risk Map",
           col = rainbow_pastel,
           axes = FALSE,
           box = FALSE)

      # Add a custom title
      title(main = "Rainbow Disease Risk Map",
            sub = "Areas with higher risk shown in warm colors",
            col.main = "#5a5a5a",
            col.sub = "#8c8c8c")
   })
}

# Run the application
shinyApp(ui = ui, server = server)
