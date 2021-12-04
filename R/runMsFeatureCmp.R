#' Shiny app launcher
#'
#' Launches the shiny app for the package (msFeatureCmp). The shiny code is
#' located in \code{./inst/shiny-scripts}.
#'
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' runMsFeatureCmp()
#' }
#'
#' @export
#' @import shiny
runMsFeatureCmp <- function() {
  appDir <- system.file("shiny-scripts", package = "msFeatureCmp")
  shiny::runApp(appDir, display.mode = "normal")
  return()
}
