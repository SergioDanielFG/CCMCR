% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_statis_dashboard.R
\name{run_statis_dashboard}
\alias{run_statis_dashboard}
\title{Launch STATIS Dual Robust Dashboard (Shiny)}
\usage{
run_statis_dashboard()
}
\description{
Launches an interactive Shiny dashboard that includes:
\itemize{
\item Phase 1 control chart (sum of robust Mahalanobis distances)
\item Phase 2 control chart (for new batches)
\item GH-Biplot visualization
}
}
\examples{
if (interactive()) {
  run_statis_dashboard()
}
}
