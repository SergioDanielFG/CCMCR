#' Classical Hotelling T2 Chart - Phase 1
#'
#' Applies the classical Hotelling T2 methodology to Phase 1 data,
#' using the sample mean and covariance matrix.
#'
#' @param data A data frame containing Phase 1 data (under control batches).
#' @param variables A character vector with the names of the quantitative variables to be used.
#'
#' @return A list with:
#'   - center: classical mean vector
#'   - covariance: classical covariance matrix
#'   - batch_statistics: data frame with T2_Stat per batch
#'   - threshold: Chi-squared control limit (0.9973 quantile)
#' @export
#' @keywords internal
hotelling_t2_phase1 <- function(data, variables) {
  batches <- unique(data$Batch)
  p <- length(variables)

  classical_center <- colMeans(data[, variables])
  classical_covariance <- cov(data[, variables])

  batch_statistics <- data.frame(Batch = character(), T2_Stat = numeric(), stringsAsFactors = FALSE)

  for (batch in batches) {
    subset_batch <- data[data$Batch == batch, variables]
    n_b <- nrow(subset_batch)
    batch_mean <- colMeans(subset_batch)
    diff <- batch_mean - classical_center
    T2_b <- n_b * t(diff) %*% solve(classical_covariance) %*% diff
    batch_statistics <- rbind(batch_statistics,
                              data.frame(Batch = as.character(batch), T2_Stat = as.numeric(T2_b)))
  }

  threshold <- qchisq(0.9973, df = p)

  return(list(
    center = classical_center,
    covariance = classical_covariance,
    batch_statistics = batch_statistics,
    threshold = threshold
  ))
}


#' Classical Hotelling T2 Chart - Phase 2
#'
#' Evaluates new batches (Phase 2) using T2 statistics based on Phase 1 estimators.
#'
#' @param new_data A data frame with new batches to evaluate (Phase 2).
#' @param variables Character vector of quantitative variables.
#' @param center Mean vector from Phase 1.
#' @param covariance Covariance matrix from Phase 1.
#'
#' @return A list with:
#'   - batch_statistics: data frame with T2_Stat per new batch
#'   - threshold: Chi-squared control limit (0.9973 quantile)
#' @export
#' @keywords internal
hotelling_t2_phase2 <- function(new_data, variables, center, covariance) {
  batches <- unique(new_data$Batch)
  p <- length(variables)

  batch_statistics <- data.frame(Batch = character(), T2_Stat = numeric(), stringsAsFactors = FALSE)

  for (batch in batches) {
    subset_batch <- new_data[new_data$Batch == batch, variables]
    n_b <- nrow(subset_batch)
    batch_mean <- colMeans(subset_batch)
    diff <- batch_mean - center
    T2_b <- n_b * t(diff) %*% solve(covariance) %*% diff
    batch_statistics <- rbind(batch_statistics,
                              data.frame(Batch = as.character(batch), T2_Stat = as.numeric(T2_b)))
  }

  threshold <- qchisq(0.9973, df = p)

  return(list(
    batch_statistics = batch_statistics,
    threshold = threshold
  ))
}

#' @examples
#' sim_batches <- simulate_pharma_batches()
#' phase1_data <- subset(sim_batches, Phase == "Phase 1")
#' phase2_data <- subset(sim_batches, Phase == "Phase 2")
#'
#' # Phase 1 classical Hotelling T2
#' t2_phase1 <- hotelling_t2_phase1(
#'   data = phase1_data,
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density")
#' )
#'
#' # Phase 2 evaluation using Phase 1 estimators
#' t2_phase2 <- hotelling_t2_phase2(
#'   new_data = phase2_data,
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density"),
#'   center = t2_phase1$center,
#'   covariance = t2_phase1$covariance
#' )
