#' Simulate Pharmaceutical Manufacturing Batches
#'
#' This function generates simulated data for a pharmaceutical manufacturing process.
#' It includes 10 in-control batches for Phase 1 and 7 additional batches for Phase 2,
#' where 4 are under control and 3 are out of control.
#'
#' @param obs_por_lote Integer. Number of observations (samples) per batch. Default is 10.
#'
#' @return A data frame with columns: Batch, Fase, Status, Concentration, Humidity, Dissolution, Density.
#'
#' @export
#'
#' @importFrom MASS mvrnorm
#' @importFrom dplyr bind_rows

simulate_pharma_batches <- function(obs_por_lote = 10) {
  # Parámetros para lotes bajo control
  mu_control <- c(98, 2.5, 300, 0.5)
  sigma_control <- matrix(c(2.0, 0.2, -5.0, 0.1,
                            0.2, 0.5, -2.0, 0.05,
                            -5.0, -2.0, 50.0, -0.5,
                            0.1, 0.05, -0.5, 0.02),
                          nrow = 4, byrow = TRUE)

  # Simular 10 lotes bajo control (Fase 1)
  control_fase1 <- dplyr::bind_rows(lapply(1:10, function(batch) {
    data.frame(
      Batch = paste0("Batch_", batch),
      Fase = "Fase 1",
      Status = "Under Control",
      MASS::mvrnorm(obs_por_lote, mu_control, sigma_control)
    )
  }))

  # Simular 4 lotes bajo control (Fase 2)
  control_fase2 <- dplyr::bind_rows(lapply(11:14, function(batch) {
    data.frame(
      Batch = paste0("Batch_", batch),
      Fase = "Fase 2",
      Status = "Under Control",
      MASS::mvrnorm(obs_por_lote, mu_control, sigma_control)
    )
  }))

  # Parámetros para lotes fuera de control
  mu_out <- c(95, 4.0, 380, 0.4)
  sigma_out <- sigma_control * 1.5

  # Simular 3 lotes fuera de control (Fase 2)
  out_fase2 <- dplyr::bind_rows(lapply(15:17, function(batch) {
    data.frame(
      Batch = paste0("Batch_", batch),
      Fase = "Fase 2",
      Status = "Out of Control",
      MASS::mvrnorm(obs_por_lote, mu_out, sigma_out)
    )
  }))

  # Unir los datos
  datos <- dplyr::bind_rows(control_fase1, control_fase2, out_fase2)
  colnames(datos)[4:7] <- c("Concentration", "Humidity", "Dissolution", "Density")

  # Convertir columnas en factores
  datos$Batch <- factor(datos$Batch)
  datos$Status <- factor(datos$Status, levels = c("Under Control", "Out of Control"))
  datos$Fase <- factor(datos$Fase, levels = c("Fase 1", "Fase 2"))

  return(datos)
}
