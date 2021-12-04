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
