#' Launch STATIS Dual Robust Dashboard (Shiny)
#'
#' Launches an interactive Shiny dashboard that includes:
#' - Phase 1 control chart (sum of robust Mahalanobis distances)
#' - Phase 2 control chart (for new batches)
#' - HJ-Biplot visualization
#'
#' @export
#' @return No return value, called for side effects (launches a Shiny application).
#' @examples
#' if (interactive()) {
#'   run_statis_dashboard()
#' }

run_statis_dashboard <- function() {
  app_dir <- system.file("shiny", "full_panel_app.R", package = "robustT2")
  if (app_dir == "") {
    stop("Shiny app script not found. Did you include inst/shiny/full_panel_app.R in the package?")
  }

  # Create a new environment to load the app into
  app_env <- new.env()
  sys.source(app_dir, envir = app_env)

  # Check that the 'app' object was defined
  if (!exists("app", envir = app_env)) {
    stop("Object 'app' not found in loaded app script.")
  }

  # Run the app from that environment
  shiny::runApp(app_env$app)
}

utils::globalVariables("app")
