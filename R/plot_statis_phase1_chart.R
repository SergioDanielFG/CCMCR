#' Plot Control Chart - Sum of Robust Mahalanobis Distances (Phase 1)
#'
#' Generates a control chart using the **sum** of robust Mahalanobis distances for each batch,
#' based on the output from `robust_statis_phase1()`.
#'
#' @param chi_sum_df A data frame with columns `Batch` and `Robust_STATIS_Distance`,
#'                   usually obtained by aggregating `phase1_result$standardized_data`.
#' @param num_vars Integer. Number of variables used in the multivariate analysis (to compute the Chi² threshold).
#' @param title Optional string. Plot title.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @import ggplot2
#' @importFrom stats qchisq
#'
#' @examples
#' data("datos_farma")
#' phase1_result <- robust_statis_phase1(
#'   subset(datos_farma, Status == "Under Control"),
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density")
#' )
#' chi_sum_df <- aggregate(Robust_STATIS_Distance ~ Batch,
#'                         data = phase1_result$standardized_data, sum)
#' plot_statis_phase1_chart_sum(chi_sum_df, num_vars = 4)

plot_statis_phase1_chart_sum <- function(chi_sum_df, num_vars, title = "Robust STATIS Dual Control Chart - Phase 1 ") {
  chi_threshold <- qchisq(0.9973, df = num_vars * 10)

  ggplot(chi_sum_df, aes(x = Batch, y = Robust_STATIS_Distance, group = 1)) +
    geom_point(size = 3, color = "#0072B2") +
    geom_line(color = "#0072B2", linewidth = 0.8) +
    geom_hline(yintercept = chi_threshold, linetype = "dashed", color = "red", linewidth = 0.8) +
    labs(
      title = title,
      x = "Batch",
      y = "Chi-Squared Statistic "
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11)
    )
}

utils::globalVariables("Robust_STATIS_Distance")
