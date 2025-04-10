#' Plot STATIS Dual Robust Control Chart - Phase 2 (All Batches)
#'
#' Plots the Chi-squared statistics for both Phase 1 and Phase 2 batches.
#'
#' @param phase1_result List from robust_statis_phase1().
#' @param phase2_result List from robust_statis_phase2().
#'
#' @return A ggplot2 object.
#' @export
#' @import ggplot2
#'
#' @examples
#' data <- simulate_pharma_batches()
#' phase1 <- robust_statis_phase1(
#'   data[data$Status == "Under Control", ],
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density"))
#' new_data <- subset(data, Status == "Out of Control")
#' phase2 <- robust_statis_phase2(
#'   new_data = new_data,
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density"),
#'   medians = phase1$global_medians,
#'   mads = phase1$global_mads,
#'   compromise_matrix = phase1$compromise_matrix,
#'   global_center = phase1$global_center)
#' plot_statis_phase2_chart(phase1, phase2)
plot_statis_phase2_chart <- function(phase1_result, phase2_result) {
  df1 <- phase1_result$batch_statistics
  df2 <- phase2_result$chi_stats_by_batch

  df1$Status <- "Under Control"
  df2$Status <- "Out of Control"

  # Igualar columnas
  df1 <- df1[, c("Batch", "Chi2_Stat", "Status")]
  df2 <- df2[, c("Batch", "Robust_STATIS_Distance", "Status")]
  colnames(df2)[colnames(df2) == "Robust_STATIS_Distance"] <- "Chi2_Stat"

  combined <- rbind(df1, df2)

  ggplot(combined, aes(x = Batch, y = Chi2_Stat, color = Status, group = Status)) +
    geom_point(size = 3) +
    geom_line() +
    geom_hline(yintercept = phase2_result$threshold, linetype = "dashed", color = "red") +
    labs(title = "Robust STATIS Dual Control Chart - Phase 2",
         x = "Batch", y = "Chi-Squared Statistic") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
}
