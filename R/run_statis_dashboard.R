#' Launch STATIS Dual Robust Dashboard (Shiny)
#'
#' Launches an interactive Shiny dashboard that includes:
#' - Phase 1 control chart (sum of robust Mahalanobis distances)
#' - Phase 2 control chart (for new batches)
#' - GH-Biplot visualization
#'
#' @export
#'
#' @examples
#' if (interactive()) {
#'   run_statis_dashboard()
#' }

run_statis_dashboard <- function() {
  app_dir <- system.file("shiny", "full_panel_app.R", package = "robustT2")
  if (app_dir == "") {
    stop("Shiny app script not found. Did you forget to install the package with inst/shiny?")
  }

  # Crea un nuevo entorno para cargar el app
  app_env <- new.env()
  sys.source(app_dir, envir = app_env)

  # Verifica si el objeto 'app' fue definido
  if (!exists("app", envir = app_env)) {
    stop("Object 'app' not found in loaded app script.")
  }

  # Ejecuta la app desde el entorno
  shiny::runApp(app_env$app)
}

utils::globalVariables("app")

