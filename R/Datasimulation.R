#' Simulate Pharmaceutical Manufacturing Batches (Realistic Variability)
#'
#' Simulates pharmaceutical manufacturing batches across two phases.
#' Phase 1 includes 10 under-control batches, each with natural variability in mean and covariance.
#' Phase 2 includes 2 clean under-control batches and 3 out-of-control batches
#' with shifted mean, increased dispersion, and moderate contamination.
#'
#' The simulated data includes four quality control variables: Concentration, Humidity, Dissolution, and Density.
#'
#' @param obs_per_batch Integer. Number of observations per batch. Default is 30.
#' @param seed Optional integer. If provided, sets a random seed for reproducibility.
#'
#' @return A data frame with 450 observations and the following columns:
#' \describe{
#'   \item{Batch}{Factor. Batch identifier (Batch_1 to Batch_15).}
#'   \item{Phase}{Factor. Phase of the process: "Phase 1" or "Phase 2".}
#'   \item{Status}{Factor. Control status: "Under Control" or "Out of Control".}
#'   \item{Concentration, Humidity, Dissolution, Density}{Numeric quality control variables.}
#' }
#' @export
#'
#' @importFrom MASS mvrnorm
#' @importFrom dplyr bind_rows
#' @importFrom stats rnorm
#' @importFrom Matrix nearPD
simulate_pharma_batches <- function(obs_per_batch = 30, seed = 780) {
  if (!is.null(seed)) set.seed(seed)

  # Reference mean and covariance under control
  mu_control <- c(98, 2.5, 300, 0.5)
  sigma_control <- matrix(c(2.0, 0.2, -5.0, 0.1,
                            0.2, 0.5, -2.0, 0.05,
                            -5.0, -2.0, 30.0, -0.5,
                            0.1, 0.05, -0.5, 0.02),
                          nrow = 4, byrow = TRUE)

  # Contamination function for Phase 2 only
  contaminate <- function(df, prop = 0.2, noise_sd = 10) {
    n <- nrow(df)
    idx_outliers <- sample(1:n, size = max(1, ceiling(prop * n)))
    if (length(idx_outliers) > 0) {
      df[idx_outliers, 1:4] <- df[idx_outliers, 1:4] + matrix(
        rnorm(length(idx_outliers) * 4, mean = 0, sd = noise_sd),
        ncol = 4
      )
    }
    df
  }

  # Phase 1: 10 under-control batches with random variability
  batches_phase1 <- dplyr::bind_rows(lapply(1:10, function(batch) {
    mu_b <- mu_control + rnorm(4, 0, c(0.4, 0.05, 4, 0.015))
    noise <- matrix(rnorm(16, 0, 0.05), nrow = 4)
    sigma_b_raw <- sigma_control + (noise + t(noise)) / 2
    diag(sigma_b_raw) <- pmax(diag(sigma_b_raw), 1e-4)
    sigma_b <- as.matrix(Matrix::nearPD(sigma_b_raw, corr = FALSE)$mat)
    data_b <- MASS::mvrnorm(obs_per_batch, mu_b, sigma_b)
    data.frame(
      Batch  = paste0("Batch_", batch),
      Phase  = "Phase 1",
      Status = "Under Control",
      data_b
    )
  }))

  # Phase 2: 2 clean under-control batches
  control_phase2 <- dplyr::bind_rows(lapply(11:12, function(batch) {
    clean_data <- MASS::mvrnorm(obs_per_batch, mu_control, sigma_control)
    data.frame(
      Batch  = paste0("Batch_", batch),
      Phase  = "Phase 2",
      Status = "Under Control",
      clean_data
    )
  }))

  # Phase 2: 3 out-of-control batches with shifted mean and increased dispersion + contamination
  mu_out_phase2    <- c(99, 2.8, 305, 0.8)
  sigma_out_phase2 <- sigma_control * 1.2
  out_phase2 <- dplyr::bind_rows(lapply(13:15, function(batch) {
    base_out <- MASS::mvrnorm(obs_per_batch, mu_out_phase2, sigma_out_phase2)
    contaminated_out <- contaminate(base_out)
    data.frame(
      Batch  = paste0("Batch_", batch),
      Phase  = "Phase 2",
      Status = "Out of Control",
      contaminated_out
    )
  }))

  # Combine and format
  sim_batches <- dplyr::bind_rows(batches_phase1, control_phase2, out_phase2)
  colnames(sim_batches)[4:7] <- c("Concentration", "Humidity", "Dissolution", "Density")
  sim_batches$Batch  <- factor(sim_batches$Batch)
  sim_batches$Status <- factor(sim_batches$Status, levels = c("Under Control", "Out of Control"))
  sim_batches$Phase  <- factor(sim_batches$Phase,  levels = c("Phase 1", "Phase 2"))

  # Truncate negative values for physical consistency
  sim_batches$Concentration <- pmax(sim_batches$Concentration, 0)
  sim_batches$Humidity      <- pmax(sim_batches$Humidity, 0)
  sim_batches$Dissolution   <- pmax(sim_batches$Dissolution, 0)
  sim_batches$Density       <- pmax(sim_batches$Density, 0)

  return(sim_batches)
}
