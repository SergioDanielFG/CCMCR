#' Classical Hotelling T2 Chart - Phase 1 (Contaminated Batches)
#'
#' Applies the classical Hotelling T2 methodology to Phase 1 data (contaminated under-control batches),
#' using standard mean and covariance estimators without robust corrections.
#'
#' @param data A data frame containing Phase 1 data (contaminated under-control batches).
#' @param variables A character vector with the names of the quantitative variables to be used.
#'
#' @return A list containing:
#' \describe{
#'   \item{center}{Classical center (mean vector) estimated from Phase 1 batches.}
#'   \item{covariance}{Classical covariance matrix estimated from Phase 1 batches.}
#'   \item{batch_statistics}{A data frame with Batch and T2_Stat (classical Hotelling-type T2 statistic).}
#'   \item{threshold}{Chi-squared control limit at the 0.9973 quantile, degrees of freedom equal to the number of variables.}
#' }
#'
#' @importFrom stats cov qchisq
#' @export
#'
#' @examples
#' # Simulate pharmaceutical manufacturing batches
#' datos <- simulate_pharma_batches()
#'
#' # Phase 1 analysis: use Phase 1 data
#' phase1_data <- subset(datos, Fase == "Fase 1")
#'
#' # Apply classical Hotelling T2 methodology
#' t2_result <- hotelling_t2_phase1(
#'   data = phase1_data,
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density")
#' )
#'
#' # View main outputs
#' t2_result$batch_statistics
#' t2_result$threshold

hotelling_t2_phase1 <- function(data, variables) {
  batches <- unique(data$Batch)
  p <- length(variables)

  # Compute classical center (mean vector) and covariance matrix
  classical_center <- colMeans(data[, variables])
  classical_covariance <- cov(data[, variables])

  # Initialize data frame for batch statistics
  batch_statistics <- data.frame(Batch = character(), T2_Stat = numeric(), stringsAsFactors = FALSE)

  # Compute Hotelling T2 statistic per batch
  for (batch in batches) {
    subset_batch <- data[data$Batch == batch, variables]
    n_b <- nrow(subset_batch)
    batch_mean <- colMeans(subset_batch)
    diff <- batch_mean - classical_center
    T2_b <- n_b * t(diff) %*% solve(classical_covariance) %*% diff
    batch_statistics <- rbind(batch_statistics, data.frame(
      Batch = as.character(batch),
      T2_Stat = as.numeric(T2_b),
      stringsAsFactors = FALSE
    ))
  }

  # Compute Chi-squared control limit
  threshold <- qchisq(0.9973, df = p)

  return(list(
    center = classical_center,
    covariance = classical_covariance,
    batch_statistics = batch_statistics,
    threshold = threshold
  ))
}
