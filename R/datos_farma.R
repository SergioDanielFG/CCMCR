#' Simulated Pharmaceutical Process Data
#'
#' This dataset contains simulated pharmaceutical manufacturing data,
#' including 10 in-control batches (Phase 1) and 7 batches in Phase 2
#' (4 under control, 3 out of control).
#'
#' Each batch contains 10 observations and 4 quantitative quality control variables.
#'
#' @format A data frame with 170 rows and 7 variables:
#' \describe{
#'   \item{Batch}{Batch identifier}
#'   \item{Status}{"Under Control" or "Out of Control"}
#'   \item{Fase}{Phase: "Fase 1" or "Fase 2"}
#'   \item{Concentration}{Concentration (mg/mL)}
#'   \item{Humidity}{Humidity (% w/w)}
#'   \item{Dissolution}{Dissolution time (seconds)}
#'   \item{Density}{Density (g/cm\eqn{^3})}
#' }
#'
#' @source Simulated using \code{simulate_pharma_batches()}.
#' @usage data("datos_farma")
#' @keywords datasets
#' @name datos_farma
#' @docType data
"datos_farma"
