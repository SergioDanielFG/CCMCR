#' Simulate Pharmaceutical Manufacturing Batches
#'
#' This function generates simulated pharmaceutical process data,
#' including in-control and out-of-control batches with 4 quality variables.
#'
#' @param num_lotes_control Number of in-control batches (default: 10).
#' @param num_lotes_fuera_control Number of out-of-control batches (default: 2).
#' @param obs_por_lote Number of observations per batch (default: 10).
#'
#' @return A data frame with the following variables:
#' \describe{
#'   \item{Batch}{Batch identifier}
#'   \item{Status}{"Under Control" or "Out of Control"}
#'   \item{Concentration}{Quality variable}
#'   \item{Humidity}{Quality variable}
#'   \item{Dissolution}{Quality variable}
#'   \item{Density}{Quality variable}
#' }
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom dplyr bind_rows
#'
#' @examples
#' datos <- simulate_pharma_batches()
#' head(datos)
simulate_pharma_batches <- function(num_lotes_control = 10, num_lotes_fuera_control = 2, obs_por_lote = 10) {

  # Media y covarianza para lotes bajo control
  mu_control <- c(98, 2.5, 300, 0.5)
  sigma_control <- matrix(c(2.0, 0.2, -5.0, 0.1,
                            0.2, 0.5, -2.0, 0.05,
                            -5.0, -2.0, 50.0, -0.5,
                            0.1, 0.05, -0.5, 0.02), nrow = 4, byrow = TRUE)

  # Lotes bajo control
  control_data <- dplyr::bind_rows(lapply(1:num_lotes_control, function(batch) {
    data.frame(Batch = paste0("Batch_", batch),
               Status = "Under Control",
               MASS::mvrnorm(obs_por_lote, mu_control, sigma_control))
  }))

  # Media y covarianza para lotes fuera de control (simulamos desviaciones)
  mu_out <- c(95, 4.0, 380, 0.4)
  sigma_out <- sigma_control * 1.5

  # Lotes fuera de control
  out_data <- dplyr::bind_rows(lapply((num_lotes_control + 1):(num_lotes_control + num_lotes_fuera_control), function(batch) {
    data.frame(Batch = paste0("Batch_", batch),
               Status = "Out of Control",
               MASS::mvrnorm(obs_por_lote, mu_out, sigma_out))
  }))

  # Unir los datos
  datos <- dplyr::bind_rows(control_data, out_data)
  colnames(datos)[3:6] <- c("Concentration", "Humidity", "Dissolution", "Density")

  return(datos)
}


