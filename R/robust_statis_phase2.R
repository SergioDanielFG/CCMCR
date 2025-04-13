#' Robust STATIS Dual - Phase 2 (New Batches Evaluation)
#'
#' Applies the robust STATIS control chart methodology to new batches using the compromise matrix
#' and robust global center obtained in Phase 1. This function computes robust Mahalanobis distances
#' for each observation in new batches (Phase 2), standardizes the data using robust estimates
#' (median and MAD) from Phase 1, and then aggregates the distances per batch as the sum of squared
#' distances. This sum serves as the Phase 2 control statistic.
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
#'   \item{standardized_data}{Data frame with standardized values and distances per observation.}
#'   \item{chi_stats_by_batch}{Chi-squared statistics per batch, calculated as the sum of squared distances.}
#'   \item{threshold}{Chi-squared control limit based on degrees of freedom = number of variables * number of observations per batch.}
#' }
#'
#' @details
#' This approach is aligned with Phase 1 methodology, where each batch is evaluated using the
#' total contribution of its observations. The Mahalanobis distances are computed with respect
#' to the global robust center and the compromise matrix obtained in Phase 1.
#'
#' @export
#' @importFrom stats mahalanobis qchisq aggregate
#'
#' @examples
#' data("datos_farma")
#' phase1 <- robust_statis_phase1(
#'   subset(datos_farma, Status == "Under Control"),
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density")
#' )
#' new_data <- subset(datos_farma, Status == "Out of Control")
#' result_phase2 <- robust_statis_phase2(
#'   new_data = new_data,
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density"),
#'   medians = phase1$global_medians,
#'   mads = phase1$global_mads,
#'   compromise_matrix = phase1$compromise_matrix,
#'   global_center = phase1$global_center
#' )
#' result_phase2$standardized_data
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
