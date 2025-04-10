#' Plot Control Chart for Phase 1 - Robust STATIS Dual
#'
#' Generates a control chart using the Chi-squared statistics for each batch
#' from the result of `robust_statis_phase1()`.
#'
#' @param phase1_result A list returned by `robust_statis_phase1()`, including `batch_statistics`.
#'
#' @return A ggplot2 object with the control chart.
#' @import ggplot2
#' @importFrom stats qchisq
#' @export
#'
#' @examples
#' data <- simulate_pharma_batches()
#' data_phase1 <- subset(data, Status == "Under Control")
#' result <- robust_statis_phase1(data_phase1,
#'            variables = c("Concentration", "Humidity", "Dissolution", "Density"))
#' plot_statis_phase1_chart(result)
plot_statis_phase1_chart <- function(phase1_result) {
  df <- phase1_result$batch_statistics
  num_vars <- ncol(phase1_result$compromise_matrix)
  chi_threshold <- qchisq(0.9973, df = num_vars * 10)

  ggplot(df, aes(x = Batch, y = Chi2_Stat, group = 1)) +
    geom_point(size = 2) +
    geom_line() +
    geom_hline(yintercept = chi_threshold, linetype = "dashed", color = "red") +
    labs(
      title = "Robust STATIS Dual Control Chart - Phase 1",
      x = "Batch",
      y = "Chi-Squared Statistic"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
}
