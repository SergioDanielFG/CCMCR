#' Simulate Pharmaceutical Manufacturing Batches (Realistic Variability)
#'
#' Simulates pharmaceutical manufacturing batches across two phases.
#' Phase 1 includes 10 under-control batches, each with natural variability in mean and covariance.
#' Phase 2 includes 2 clean under-control batches and 3 out-of-control batches
#' with shifted mean, increased dispersion, and moderate contamination.
#'
#' The simulated data includes four quality control variables: Concentration, Humidity, Dissolution, and Density.
#'
#' @param obs_por_lote Integer. Number of observations per batch. Default is 10.
#' @param seed Optional integer. If provided, sets a random seed for reproducibility.
#'
#' @return A data frame with 150 observations and the following columns:
#' \describe{
#'   \item{Batch}{Factor. Batch identifier (Batch_1 to Batch_15).}
#'   \item{Fase}{Factor. Phase of the process: "Fase 1" or "Fase 2".}
#'   \item{Status}{Factor. Control status: "Under Control" or "Out of Control".}
#'   \item{Concentration, Humidity, Dissolution, Density}{Numeric quality control variables.}
#' }
#' @export
#'
#' @importFrom MASS mvrnorm
#' @importFrom dplyr bind_rows
#' @importFrom stats rnorm
#' @importFrom Matrix nearPD

simulate_pharma_batches <- function(obs_por_lote = 10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Reference mean and covariance under control
  mu_control <- c(98, 2.5, 300, 0.5)
  sigma_control <- matrix(c(2.0, 0.2, -5.0, 0.1,
                            0.2, 0.5, -2.0, 0.05,
                            -5.0, -2.0, 30.0, -0.5,
                            0.1, 0.05, -0.5, 0.02),
                          nrow = 4, byrow = TRUE)

  # Contamination function for Phase 2 only
  contaminate <- function(data, prop = 0.2, noise_sd = 10) {
    n <- nrow(data)
    idx_outliers <- sample(1:n, size = max(1, ceiling(prop * n)))
    if (length(idx_outliers) > 0) {
      data[idx_outliers, 1:4] <- data[idx_outliers, 1:4] + matrix(
        rnorm(length(idx_outliers) * 4, mean = 0, sd = noise_sd),
        ncol = 4
      )
    }
    return(data)
  }

  # Phase 1: 10 under-control batches with random variability
  batches_fase1 <- dplyr::bind_rows(lapply(1:10, function(batch) {
    mu_b <- mu_control + rnorm(4, 0, c(0.4, 0.05, 4, 0.015))
    noise <- matrix(rnorm(16, 0, 0.05), nrow = 4)
    sigma_b_raw <- sigma_control + (noise + t(noise)) / 2
    diag(sigma_b_raw) <- pmax(diag(sigma_b_raw), 1e-4)
    sigma_b <- as.matrix(Matrix::nearPD(sigma_b_raw, corr = FALSE)$mat)
    data_b <- MASS::mvrnorm(obs_por_lote, mu_b, sigma_b)
    data.frame(
      Batch = paste0("Batch_", batch),
      Fase = "Fase 1",
      Status = "Under Control",
      data_b
    )
  }))

  # Phase 2: 2 clean under-control batches
  control_fase2 <- dplyr::bind_rows(lapply(11:12, function(batch) {
    clean_data <- MASS::mvrnorm(obs_por_lote, mu_control, sigma_control)
    data.frame(
      Batch = paste0("Batch_", batch),
      Fase = "Fase 2",
      Status = "Under Control",
      clean_data
    )
  }))

  # Phase 2: 3 out-of-control batches with shifted mean and increased dispersion + contamination
  mu_out_fase2 <- c(90, 4.0, 315, 1.8)
  sigma_out_fase2 <- sigma_control * 0.5
  out_fase2 <- dplyr::bind_rows(lapply(13:15, function(batch) {
    base_out <- MASS::mvrnorm(obs_por_lote, mu_out_fase2, sigma_out_fase2)
    contaminated_out <- contaminate(base_out)
    data.frame(
      Batch = paste0("Batch_", batch),
      Fase = "Fase 2",
      Status = "Out of Control",
      contaminated_out
    )
  }))

  # Combine and format
  datos <- dplyr::bind_rows(batches_fase1, control_fase2, out_fase2)
  colnames(datos)[4:7] <- c("Concentration", "Humidity", "Dissolution", "Density")
  datos$Batch <- factor(datos$Batch)
  datos$Status <- factor(datos$Status, levels = c("Under Control", "Out of Control"))
  datos$Fase <- factor(datos$Fase, levels = c("Fase 1", "Fase 2"))

  # Truncate negative values for physical consistency
  datos$Concentration <- pmax(datos$Concentration, 0)
  datos$Humidity <- pmax(datos$Humidity, 0)
  datos$Dissolution <- pmax(datos$Dissolution, 0)
  datos$Density <- pmax(datos$Density, 0)

  return(datos)
}

