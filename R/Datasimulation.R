#' Simulate Pharmaceutical Manufacturing Batches
#'
#' This function generates simulated data for a pharmaceutical manufacturing process.
#' It includes both in-control and out-of-control batches with 4 continuous quality variables:
#' Concentration, Humidity, Dissolution, and Density.
#'
#' @param num_lotes_control Integer. Number of in-control batches to simulate. Default is 10.
#' @param num_lotes_fuera_control Integer. Number of out-of-control batches to simulate. Default is 2.
#' @param obs_por_lote Integer. Number of observations (samples) per batch. Default is 10.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{Batch}{Factor indicating the batch number (e.g., Batch_1, Batch_2, ...)}
#'   \item{Status}{Factor indicating whether the batch is "Under Control" or "Out of Control"}
#'   \item{Concentration}{Simulated concentration values}
#'   \item{Humidity}{Simulated humidity values}
#'   \item{Dissolution}{Simulated dissolution values}
#'   \item{Density}{Simulated density values}
#' }
#'
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom dplyr bind_rows
#'
#' @examples
#' data <- simulate_pharma_batches()
#' head(data)
simulate_pharma_batches <- function(num_lotes_control = 10, num_lotes_fuera_control = 2, obs_por_lote = 10) {

  mu_control <- c(98, 2.5, 300, 0.5)
  sigma_control <- matrix(c(2.0, 0.2, -5.0, 0.1,
                            0.2, 0.5, -2.0, 0.05,
                            -5.0, -2.0, 50.0, -0.5,
                            0.1, 0.05, -0.5, 0.02), nrow = 4, byrow = TRUE)

  control_data <- dplyr::bind_rows(lapply(1:num_lotes_control, function(batch) {
    data.frame(Batch = paste0("Batch_", batch),
               Status = "Under Control",
               MASS::mvrnorm(obs_por_lote, mu_control, sigma_control))
  }))

  mu_out <- c(95, 4.0, 380, 0.4)
  sigma_out <- sigma_control * 1.5

  out_data <- dplyr::bind_rows(lapply((num_lotes_control + 1):(num_lotes_control + num_lotes_fuera_control), function(batch) {
    data.frame(Batch = paste0("Batch_", batch),
               Status = "Out of Control",
               MASS::mvrnorm(obs_por_lote, mu_out, sigma_out))
  }))

  datos <- dplyr::bind_rows(control_data, out_data)
  colnames(datos)[3:6] <- c("Concentration", "Humidity", "Dissolution", "Density")

  # Convert 'Status' to factor
  datos$Status <- factor(datos$Status, levels = c("Under Control", "Out of Control"))

  return(datos)
}


