# app.R
# Package: msFeatureCmp
# Author: Yijia Chen
# Date: 2021-12-04
# Version: 0.1.0

# This file contains the source for the package's shiny app.

library("shiny")

ui <- fluidPage(
  shiny::titlePanel("msFeatureCmp: Mass Spectrometry Feature Comparator"),

  shiny::sidebarLayout(
    shiny::sidebarPanel("sidebar panel"),
    shiny::mainPanel(
      h4("Output", align = "center")
    )
  )
)

server <- function(input, output) {

}

shiny::shinyApp(ui = ui, server = server)

# [END]
