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
#'   # Launches the full STATIS Dual Robust monitoring dashboard with:
#'   # - Phase 1 control chart
#'   # - Phase 2 control chart
#'   # - GH-Biplot
#'   #
#'   # Uses the internal dataset: datos_farma
#'   run_statis_dashboard()
#' }

run_statis_dashboard <- function() {
  app_dir <- system.file("shiny", "full_panel_app.R", package = "CCMCR")
  if (app_dir == "") {
    stop("Shiny app script not found. Did you forget to install the package with inst/?")
  }
  source(app_dir, local = TRUE)
  shiny::runApp(app)
}

utils::globalVariables("app")

