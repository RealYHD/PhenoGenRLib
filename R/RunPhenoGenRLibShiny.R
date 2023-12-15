#' Launch Shiny App For PhenoGenRLib
#'
#' Launches Shiny app for this package. The application allows for loading
#' VCFs, linking said VCFs to their case metadata, resolving RSIDs based on
#' chromosome and position, and running associative case-control based on
#' Ensembl population. The app also allows for easier visualization of variance
#' distribution based on different nucleobases.
#'
#' @return Opens a Shiny application but does not return anything.
#'
#' @examples
#' \dontrun{
#' runPhenoGenRLib()
#' }
#'
#'
#' @export
#' @importFrom shiny runApp
runPhenoGenRLib <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "PhenoGenRLib")
  shiny::runApp(appDir, display.mode = "normal")
  return()
}

# [END]
