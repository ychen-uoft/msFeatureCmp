# app.R
# Package: msFeatureCmp
# Author: Yijia Chen
# Date: 2021-12-05
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
      br(),
      verbatimTextOutput("vtextOutput"),
      br(),
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
    if (input$cmpOptions == 1) {
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

      # Do the same for feature set A
      if (is.null(input$featureFileA)) {
        output$textOutput <- renderText({
          "Error: no featureXML file for feature set A has been selected"
        })
      }
      featureFilePathA <- input$featureFileA$datapath
      req(featureFilePathA)

      fileExt <- tools::file_ext(featureFilePathA)
      if (fileExt != "featureXML") {
        output$textOutput <- renderText({
          "Error: feature set A file must in featureXML format"
        })
      }
      req(fileExt == "featureXML")

      # Do the same for feature set B
      if (is.null(input$featureFileB)) {
        output$textOutput <- renderText({
          "Error: no featureXML file for feature set B has been selected"
        })
      }
      featureFilePathB <- input$featureFileB$datapath
      req(featureFilePathB)

      fileExt <- tools::file_ext(featureFilePathB)
      if (fileExt != "featureXML") {
        output$textOutput <- renderText({
          "Error: feature set B file must in featureXML format"
        })
      }
      req(fileExt == "featureXML")

      # Interface with the main package and capture the output
      output$textOutput <- renderText({
        "Current status: comparing features between sets"
      })
      cmpResults <- capture.output(msFeatureCmp::compareFeatures(rawDataFilePath,
                                                                 featureFilePathA,
                                                                 featureFilePathB))

      output$textOutput <- renderText({
        "Current status: done comparing features between sets"
      })
      # Captured output is stored per-line in a character vector, so concatenate
      # everything into a single string for output
      displayStr <- "Feature comparison results:"
      for (line in cmpResults) {
        displayStr <- paste(displayStr, line, sep = "\n")
      }
      output$vtextOutput <- renderText({
        displayStr
      })
    }
    else if (input$cmpOptions == 2) {
      if (is.null(input$featureFileA)) {
        output$textOutput <- renderText({
          "Error: no featureXML file for feature set A has been selected"
        })
      }
      featureFilePath <- input$featureFileA$datapath
      req(featureFilePath)

      fileExt <- tools::file_ext(featureFilePath)
      if (fileExt != "featureXML") {
        output$textOutput <- renderText({
          "Error: feature set A file must in featureXML format"
        })
      }
      req(fileExt == "featureXML")

      output$textOutput <- renderText({
        "Current status: finding feature in feature set A"
      })
      featureInfo <- capture.output(msFeatureCmp::getFeatureByIdx(
        featureFilePath, input$featureIdx))

      output$textOutput <- renderText({
        "Current status: done looking for feature in feature set A.
        The three values below are retention time (sec), mass-to-charge (Th),
        and signal intensity (unitless)"
      })
      output$vtextOutput <- renderText({
        featureInfo
      })
    }
    else if (input$cmpOptions == 3) {
      if (is.null(input$featureFileB)) {
        output$textOutput <- renderText({
          "Error: no featureXML file for feature set A has been selected"
        })
      }
      featureFilePath <- input$featureFileB$datapath
      req(featureFilePath)

      fileExt <- tools::file_ext(featureFilePath)
      if (fileExt != "featureXML") {
        output$textOutput <- renderText({
          "Error: feature set A file must in featureXML format"
        })
      }
      req(fileExt == "featureXML")

      output$textOutput <- renderText({
        "Current status: finding feature in feature set A"
      })
      featureInfo <- capture.output(msFeatureCmp::getFeatureByIdx(
        featureFilePath, input$featureIdx))

      output$textOutput <- renderText({
        "Current status: done looking for feature in feature set A.
        The three values below are retention time (sec), mass-to-charge (Th),
        and signal intensity (unitless)"
      })
      output$vtextOutput <- renderText({
        featureInfo
      })
    }
  })

  # TODO: reduce duplicate code
  # Plot button clicked (in Plot)
  observeEvent(input$runPlt, {
    # Raw data
    if (input$pltOptions == 1) {
      if (is.null(input$rawDataFile)) {
        output$textOutput <- renderText({
          "Error: no mzML file has been selected"
        })
      }
      rawDataFilePath <- input$rawDataFile$datapath
      req(rawDataFilePath)

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
      # Interface with the main package and plot the results
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
      if (is.null(input$featureFileA)) {
        output$textOutput <- renderText({
          "Error: no featureXML file for feature set A has been selected"
        })
      }
      featureFilePath <- input$featureFileA$datapath
      req(featureFilePath)

      fileExt <- tools::file_ext(featureFilePath)
      if (fileExt != "featureXML") {
        output$textOutput <- renderText({
          "Error: feature set A file must in featureXML format"
        })
      }
      req(fileExt == "featureXML")

      output$textOutput <- renderText({
        "Current status: plotting feature set A, please wait..."
      })
      plot <- msFeatureCmp::plotSingleFeatureSet(featureFilePath)

      output$textOutput <- renderText({
        "Current status: done plotting feature set A"
      })
      output$plotOutput <- renderPlot({
        plot
      })
    }
    #Feature set B
    else if (input$pltOptions == 3) {
      if (is.null(input$featureFileB)) {
        output$textOutput <- renderText({
          "Error: no featureXML file for feature set B has been selected"
        })
      }
      featureFilePath <- input$featureFileB$datapath
      req(featureFilePath)

      fileExt <- tools::file_ext(featureFilePath)
      if (fileExt != "featureXML") {
        output$textOutput <- renderText({
          "Error: feature set B file must in featureXML format"
        })
      }
      req(fileExt == "featureXML")

      output$textOutput <- renderText({
        "Current status: plotting feature set B, please wait..."
      })
      plot <- msFeatureCmp::plotSingleFeatureSet(featureFilePath)

      output$textOutput <- renderText({
        "Current status: done plotting feature set B"
      })
      output$plotOutput <- renderPlot({
        plot
      })
    }
    # Both feature sets
    else if (input$pltOptions == 4) {
      # Load feature set A first
      if (is.null(input$featureFileA)) {
        output$textOutput <- renderText({
          "Error: no featureXML file for feature set A has been selected"
        })
      }
      featureFilePathA <- input$featureFileA$datapath
      req(featureFilePathA)

      fileExt <- tools::file_ext(featureFilePathA)
      if (fileExt != "featureXML") {
        output$textOutput <- renderText({
          "Error: feature set A file must in featureXML format"
        })
      }
      req(fileExt == "featureXML")

      # Next load feature set B
      if (is.null(input$featureFileB)) {
        output$textOutput <- renderText({
          "Error: no featureXML file for feature set B has been selected"
        })
      }
      featureFilePathB <- input$featureFileB$datapath
      req(featureFilePathB)

      fileExt <- tools::file_ext(featureFilePathB)
      if (fileExt != "featureXML") {
        output$textOutput <- renderText({
          "Error: feature set B file must in featureXML format"
        })
      }
      req(fileExt == "featureXML")

      # Plot the two feature sets
      output$textOutput <- renderText({
        "Current status: plotting both feature sets, please wait..."
      })
      plot <- msFeatureCmp::plotTwoFeatureSets(featureFilePathA,
                                               featureFilePathB)

      output$textOutput <- renderText({
        "Current status: done plotting both feature sets"
      })
      output$plotOutput <- renderPlot({
        plot
      })
    }
    # Feature set A on raw data
    else if (input$pltOptions == 5) {
      # Load raw data file first
      if (is.null(input$rawDataFile)) {
        output$textOutput <- renderText({
          "Error: no mzML file has been selected"
        })
      }
      rawDataFilePath <- input$rawDataFile$datapath
      req(rawDataFilePath)

      fileExt <- tools::file_ext(rawDataFilePath)
      if (fileExt != "mzML") {
        output$textOutput <- renderText({
          "Error: raw data file must be in mzML format"
        })
      }
      req(fileExt == "mzML")

      # Next load feature set A
      if (is.null(input$featureFileA)) {
        output$textOutput <- renderText({
          "Error: no featureXML file for feature set A has been selected"
        })
      }
      featureFilePath <- input$featureFileA$datapath
      req(featureFilePath)

      fileExt <- tools::file_ext(featureFilePath)
      if (fileExt != "featureXML") {
        output$textOutput <- renderText({
          "Error: feature set A file must in featureXML format"
        })
      }
      req(fileExt == "featureXML")

      # Plot feature set A on the raw data
      output$textOutput <- renderText({
        "Current status: plotting feature set A on the raw data, please wait..."
      })
      plot <- msFeatureCmp::plotFeatureSetOnRawData(rawDataFilePath,
                                                    featureFilePath)

      output$textOutput <- renderText({
        "Current status: done plotting feature set A on the raw data"
      })
      output$plotOutput <- renderPlot({
        plot
      })
    }
    # Feature set B on raw data
    else if (input$pltOptions == 6) {
      # Load raw data file first
      if (is.null(input$rawDataFile)) {
        output$textOutput <- renderText({
          "Error: no mzML file has been selected"
        })
      }
      rawDataFilePath <- input$rawDataFile$datapath
      req(rawDataFilePath)

      fileExt <- tools::file_ext(rawDataFilePath)
      if (fileExt != "mzML") {
        output$textOutput <- renderText({
          "Error: raw data file must be in mzML format"
        })
      }
      req(fileExt == "mzML")

      # Next load feature set B
      if (is.null(input$featureFileB)) {
        output$textOutput <- renderText({
          "Error: no featureXML file for feature set B has been selected"
        })
      }
      featureFilePath <- input$featureFileB$datapath
      req(featureFilePath)

      fileExt <- tools::file_ext(featureFilePath)
      if (fileExt != "featureXML") {
        output$textOutput <- renderText({
          "Error: feature set B file must in featureXML format"
        })
      }
      req(fileExt == "featureXML")

      # Plot feature set B on the raw data
      output$textOutput <- renderText({
        "Current status: plotting feature set B on the raw data, please wait..."
      })
      plot <- msFeatureCmp::plotFeatureSetOnRawData(rawDataFilePath,
                                                    featureSetFilePath)

      output$textOutput <- renderText({
        "Current status: done plotting feature set B on the raw data"
      })
      output$plotOutput <- renderPlot({
        plot
      })
    }
  })

  return(invisible(NULL))
}

shinyApp(ui = ui, server = server)

# [END]
