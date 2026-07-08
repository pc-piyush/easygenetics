#' @export

launch_app <- function() {
  #shiny::runApp("easygenetics", "easygenetics")
  #runApp('R')
   shiny::runApp(system.file('shinyApp', package='easygenetics'))
}

