#' Plot Control Chart - Robust STATIS Dual (Phase 1)
#'
#' Plots the Chi-squared statistic per batch based on the **sum of robust Mahalanobis distances**
#' calculated in `robust_statis_phase1()`.
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
#'
#' @examples
#' data("datos_farma")
#'
#' # Filtrar lotes bajo control de la Fase 1
#' phase1_result <- robust_statis_phase1(
#'   subset(datos_farma, Fase == "Fase 1" & Status == "Under Control"),
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density")
#' )
#'
#' # Graficar los estadísticos Chi² por lote
#' plot_statis_phase1_chart(
#'   batch_statistics = phase1_result$batch_statistics,
#'   num_vars = 4
#' )


plot_statis_phase1_chart <- function(batch_statistics, num_vars, title = "Robust STATIS Dual Control Chart - Phase 1") {
  # Supone 10 observaciones por lote
  chi_threshold <- qchisq(0.9973, df = num_vars * 10)

  ggplot(batch_statistics, aes(x = Batch, y = Chi2_Stat, group = 1)) +
    geom_point(size = 3, color = "#0072B2") +
    geom_line(color = "#0072B2", linewidth = 0.8) +

    # etiquetas sobre cada punto
    geom_text(
      aes(label = round(Chi2_Stat, 1)),
      color = "black",
      vjust = -0.8,
      size = 3.2,
      show.legend = FALSE
    ) +

    # Línea de umbral + etiqueta del valor del umbral
    geom_hline(yintercept = chi_threshold, linetype = "dashed", color = "red", linewidth = 0.8) +
    annotate("text",
             x = Inf, y = chi_threshold,
             label = paste0("UCL = ", round(chi_threshold, 1)),
             hjust = 1.2, vjust = -0.5, color = "red", size = 4) +

    labs(
      title = title,
      x = "Batch",
      y = expression(chi^2 ~ "(Sum of Squared Robust Mahalanobis Distances per Batch)")
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11)
    )
}

utils::globalVariables(c("Chi2_Stat", "Batch"))
