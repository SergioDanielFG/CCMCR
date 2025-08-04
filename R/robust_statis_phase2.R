#' Robust STATIS Dual - Phase 2 (New Batches Evaluation)
#'
#' Applies the robust STATIS Dual control chart methodology to evaluate new batches,
#' using the compromise matrix and the global robust center obtained in Phase 1.
#' Each batch is summarized using a robust Hotelling-type \( T^2 \) statistic.
#'
#' @param new_data A data frame containing the new batches to evaluate.
#' @param variables Character vector with the names of the variables to be used.
#' @param medians Named numeric vector containing the global medians obtained in Phase 1.
#' @param mads Named numeric vector containing the scaled MADs obtained in Phase 1.
#' @param compromise_matrix Robust compromise matrix computed in Phase 1.
#' @param global_center Robust global center obtained in Phase 1.
#'
#' @return A list containing:
#' \describe{
#'   \item{standardized_data}{Data frame with the new batches standardized using the global medians and scaled MADs.}
#'   \item{t2_stats_by_batch}{Data frame with the Hotelling-type \( T^2 \) statistics per batch.}
#'   \item{threshold}{Control limit based on the Chi-squared distribution (0.9973 quantile, degrees of freedom equal to the number of variables).}
#' }
#'
#' @export
#' @keywords internal
#' @importFrom stats qchisq
#'
#' @examples
#' # Simulate new pharmaceutical manufacturing batches
#' datos <- simulate_pharma_batches()
#'
#' # Phase 1 analysis: use only Phase 1 and under control batches
#' phase1_data <- subset(datos, Fase == "Fase 1" & Status == "Under Control")
#' phase1 <- robust_statis_phase1(
#'   data = phase1_data,
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density")
#' )
#'
#' # Phase 2 analysis: evaluate new batches (Phase 2)
#' new_data <- subset(datos, Fase == "Fase 2")
#' result_phase2 <- robust_statis_phase2(
#'   new_data = new_data,
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density"),
#'   medians = phase1$global_medians,
#'   mads = phase1$global_mads,
#'   compromise_matrix = phase1$compromise_matrix,
#'   global_center = phase1$global_center
#' )
#'
#' # View main outputs
#' result_phase2$t2_stats_by_batch
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

  # Global standardization
  for (var in variables) {
    standardized_data[[var]] <- (new_data[[var]] - medians[[var]]) / mads[[var]]
  }

  batches <- unique(standardized_data$Batch)
  t2_stats_by_batch <- data.frame(Batch = character(), T2_Stat = numeric())

  for (batch in batches) {
    subset_batch <- standardized_data[standardized_data$Batch == batch, variables]
    n_b <- nrow(subset_batch)
    x_b <- rrcov::CovMcd(subset_batch, alpha = 0.75)@center
    diff <- x_b - global_center
    T2_b <- n_b * t(diff) %*% solve(compromise_matrix) %*% diff
    t2_stats_by_batch <- rbind(t2_stats_by_batch, data.frame(
      Batch = batch,
      T2_Stat = as.numeric(T2_b)
    ))
  }

  threshold <- qchisq(0.9973, df = length(variables))

  return(list(
    standardized_data = standardized_data,
    t2_stats_by_batch = t2_stats_by_batch,
    threshold = threshold
  ))
}
