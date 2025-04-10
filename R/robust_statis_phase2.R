#' Robust STATIS Dual - Phase 2 (New Batches Evaluation)
#'
#' Applies the robust STATIS control chart methodology to new batches using the compromise matrix
#' and robust global center obtained in Phase 1.
#'
#' @param new_data A data frame containing new batches to evaluate (usually out-of-control).
#' @param variables Character vector of variable names to use.
#' @param medians Named numeric vector of medians from Phase 1 (one per variable).
#' @param mads Named numeric vector of MADs from Phase 1 (one per variable).
#' @param compromise_matrix Robust compromise covariance matrix from Phase 1.
#' @param global_center Robust global center (vector) from Phase 1.
#'
#' @return A list with:
#' \describe{
#'   \item{standardized_data}{Data frame with standardized values and distances.}
#'   \item{chi_stats_by_batch}{Chi-square statistics aggregated by batch.}
#'   \item{threshold}{Chi-square threshold used.}
#' }
#' @export
#' @importFrom stats mahalanobis median mad aggregate
#'
#' @examples
#' data <- simulate_pharma_batches()
#' phase1 <- robust_statis_phase1(
#' data,
#' variables = c("Concentration", "Humidity", "Dissolution", "Density"))
#' new_data <- subset(data, Status == "Out of Control")
#' result_phase2 <- robust_statis_phase2(
#'   new_data = new_data,
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density"),
#'   medians = phase1$global_medians,
#'   mads = phase1$global_mads,
#'   compromise_matrix = phase1$compromise_matrix,
#'   global_center = phase1$global_center
#' )
robust_statis_phase2 <- function(new_data,
                                 variables,
                                 medians,
                                 mads,
                                 compromise_matrix,
                                 global_center) {
  stopifnot(all(variables %in% colnames(new_data)))
  stopifnot(all(variables %in% names(medians)))
  stopifnot(all(variables %in% names(mads)))

  standardized_data <- new_data

  for (var in variables) {
    standardized_data[[var]] <- (new_data[[var]] - medians[[var]]) / mads[[var]]
  }

  X_phase2 <- as.matrix(standardized_data[, variables])
  distances <- mahalanobis(
    x = X_phase2,
    center = global_center,
    cov = compromise_matrix
  )

  standardized_data$Robust_STATIS_Distance <- distances

  chi_stats_by_batch <- aggregate(Robust_STATIS_Distance ~ Batch, data = standardized_data, FUN = sum)
  threshold <- qchisq(0.9973, df = length(variables) * 10)

  return(list(
    standardized_data = standardized_data,
    chi_stats_by_batch = chi_stats_by_batch,
    threshold = threshold
  ))
}
