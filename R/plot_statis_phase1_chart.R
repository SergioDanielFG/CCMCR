#' Plot Control Chart - Robust STATIS Dual (Phase 1)
#'
#' Plots the Chi-squared statistic per batch using the robust Hotelling T² statistic
#' calculated in `robust_statis_phase1()`. The control limit is based on a Chi-squared
#' distribution with degrees of freedom equal to the number of variables.
#'
#' @param batch_statistics A data frame with columns `Batch` and `Chi2_Stat`,
#'                         typically from `phase1_result$batch_statistics`.
#' @param num_vars Integer. Number of variables used in the multivariate analysis (to compute the Chi² threshold).
#' @param title Optional string. Plot title.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @import ggplot2
#' @importFrom stats qchisq
#' @importFrom forcats fct_inorder
#'
#' @examples
#' data("datos_farma")
#' phase1_result <- robust_statis_phase1(
#'   subset(datos_farma, Fase == "Fase 1" & Status == "Under Control"),
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density")
#' )
#' plot_statis_phase1_chart(
#'   batch_statistics = phase1_result$batch_statistics,
#'   num_vars = 4
#' )

plot_statis_phase1_chart <- function(batch_statistics, num_vars,
                                     title = "Robust STATIS Dual Control Chart - Phase 1") {
  chi_threshold <- qchisq(0.9973, df = num_vars)

  batch_statistics$Batch <- forcats::fct_inorder(batch_statistics$Batch)

  ggplot(batch_statistics, aes(x = Batch, y = Chi2_Stat, group = 1)) +
    geom_point(size = 3, color="#00C8D7" ) +
    geom_line(linewidth = 0.8, color="#00C8D7") +
    geom_text(
      aes(label = round(Chi2_Stat, 1)),
      color = "black", vjust = -0.8, size = 3.2, show.legend = FALSE
    ) +
    geom_hline(yintercept = chi_threshold, linetype = "dashed", color = "red", linewidth = 0.8) +
    annotate("text",
             x = Inf, y = chi_threshold,
             label = paste0("UCL = ", round(chi_threshold, 1)),
             hjust = 1.2, vjust = -0.5, color = "red", size = 4) +
    labs(
      title = title,
      x = "Batch",
      y = expression(chi^2 ~ "(Hotelling~T^2~per~Batch)")
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11)
    )
}

utils::globalVariables(c("Chi2_Stat", "Batch"))
