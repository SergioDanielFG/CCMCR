#' Robust STATIS Dual - Phase 2 (New Batches Evaluation)
#'
#' Applies the robust STATIS control chart methodology to new batches using the compromise matrix
#' and robust global center obtained in Phase 1. Each batch is summarized by a robust Hotelling T²-type
#' statistic using the standardized batch mean and the Phase 1 compromise matrix.
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
#'   \item{standardized_data}{Data frame with standardized values.}
#'   \item{chi_stats_by_batch}{Chi-squared statistics per batch using Hotelling-type T² formulation.}
#'   \item{threshold}{Chi-squared control limit based on degrees of freedom = number of variables.}
#' }
#'
#' @export
#' @importFrom stats qchisq
#'
#' @examples
#' data("datos_farma")
#' phase1 <- robust_statis_phase1(
#'   subset(datos_farma, Fase == "Fase 1" & Status == "Under Control"),
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density")
#' )
#' new_data <- subset(datos_farma, Fase == "Fase 2")
#' result_phase2 <- robust_statis_phase2(
#'   new_data = new_data,
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density"),
#'   medians = phase1$global_medians,
#'   mads = phase1$global_mads,
#'   compromise_matrix = phase1$compromise_matrix,
#'   global_center = phase1$global_center
#' )
#' result_phase2$chi_stats_by_batch
#' result_phase2$threshold

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

  # Calcular estadístico Hotelling T² por lote
  batches <- unique(standardized_data$Batch)
  chi_stats_by_batch <- data.frame(Batch = character(), Chi2_Stat = numeric())

  for (batch in batches) {
    subset_batch <- standardized_data[standardized_data$Batch == batch, variables]
    n_b <- nrow(subset_batch)
    x_b <- colMeans(subset_batch)
    diff <- x_b - global_center
    T2_b <- n_b * t(diff) %*% solve(compromise_matrix) %*% diff
    chi_stats_by_batch <- rbind(chi_stats_by_batch, data.frame(
      Batch = batch,
      Chi2_Stat = as.numeric(T2_b)
    ))
  }

  # Umbral basado en distribución Chi²
  threshold <- qchisq(0.9973, df = length(variables))

  return(list(
    standardized_data = standardized_data,
    chi_stats_by_batch = chi_stats_by_batch,
    threshold = threshold
  ))
}
