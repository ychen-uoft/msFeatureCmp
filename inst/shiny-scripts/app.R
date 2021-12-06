# app.R
# Package: msFeatureCmp
# Author: Yijia Chen
# Date: 2021-12-04
# Version: 0.1.0

# This file contains the source for the package's shiny app.

library("shiny")

ui <- fluidPage(
  titlePanel("msFeatureCmp: Mass Spectrometry Feature Comparator"),

  # Sidebar takes user input (file locations, which functions to run, etc.)
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        # Tab for loading required files
        tabPanel("Load data",
                 br(),
                 helpText("Select files to load:"),
                 fileInput("rawDataFile", h4("Raw mzML data file:")),
                 fileInput("featureFileA", h4("First featureXML file (A):")),
                 fileInput("featureFileB", h4("Second featureXML file (B):"))),

        # Tab for running a computational function
        tabPanel("Compare",
                 br(),
                 helpText("Choose a computational operation to run:"),
                 radioButtons("cmpOptions", h4("Operations"),
                              choices = list("Compare both feature sets" = 1,
                                             "Get feature info from set A" = 2,
                                             "Get feature info from set B" = 3)),
                 numericInput("featureIdx", h5("Feature index"), value = 0),
                 actionButton("runCmp", "Run")),

        # Tab for running a visualization function
        tabPanel("Plot",
                 br(),
                 helpText("Choose a plotting operation to run:"),
                 radioButtons("pltOptions", h4("Plotting operations:"),
                              choices = list("Raw data" = 1,
                                             "Feature set A" = 2,
                                             "Feature set B" = 3,
                                             "Both feature sets" = 4,
                                             "Feature set A on raw data" = 5,
                                             "Feature set B on raw data" = 6)),
                 actionButton("runPlt", "Plot"))
      ),
      width = 5
    ),

    # Function results will be displayed in the main panel
    mainPanel(
      h4("Output", align = "center"),
      textOutput("textOutput"),
      tableOutput("tableOutput"),
      plotOutput("plotOutput"),
      width = 7
    )
  )
)

server <- function(input, output) {
  output$textOutput <- renderText({
    "Current status: OK"
  })

  # Run button clicked (in Compare)
  observeEvent(input$runCmp, {

  })

  # Plot button clicked (in Plot)
  observeEvent(input$runPlt, {
    # Raw data
    if (input$pltOptions == 1) {
      # Check if a file has been uploaded
      if (is.null(input$rawDataFile)) {
        output$textOutput <- renderText({
          "Error: no mzML file has been selected"
        })
      }
      rawDataFilePath <- input$rawDataFile$datapath
      req(rawDataFilePath)

      # Check if the file is in mzML format
      fileExt <- tools::file_ext(rawDataFilePath)
      if (fileExt != "mzML") {
        output$textOutput <- renderText({
          "Error: raw data file must be in mzML format"
        })
      }
      req(fileExt == "mzML")

      output$textOutput <- renderText({
        "Current status: plotting raw data, please wait..."
      })
      plot <- msFeatureCmp::plotRawData(rawDataFilePath)

      output$textOutput <- renderText({
        "Current status: done plotting raw data"
      })
      output$plotOutput <- renderPlot({
        plot
      })
    }
    # Feature set A
    else if (input$pltOptions == 2) {

    }
    # Feature set B
    else if (input$pltOptions == 3) {

    }
    # Both feature sets
    else if (input$pltOptions == 4) {

    }
    # Feature set A on raw data
    else if (input$pltOptions == 5) {

    }
    # Feature set B on raw data
    else if (input$pltOptions == 6) {

    }
  })

  return(invisible(NULL))
}

shinyApp(ui = ui, server = server)

# [END]
