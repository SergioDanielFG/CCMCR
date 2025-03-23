#' Simulated Pharmaceutical Process Data
#'
#' This dataset contains simulated pharmaceutical manufacturing data,
#' including 10 in-control batches and 2 out-of-control batches.
#'
#' Each batch contains 10 observations and 4 quality control variables:
#' Concentration, Humidity, Dissolution, and Density.
#'
#' @format A data frame with 120 rows and 6 variables:
#' \describe{
#'   \item{Batch}{Batch identifier}
#'   \item{Status}{"Under Control" or "Out of Control"}
#'   \item{Concentration}{Quantitative quality variable}
#'   \item{Humidity}{Quantitative quality variable}
#'   \item{Dissolution}{Quantitative quality variable}
#'   \item{Density}{Quantitative quality variable}
#' }
#' @source Simulated using \code{simulate_pharma_batches()} from the CCMCR package.
"datos_farma"
