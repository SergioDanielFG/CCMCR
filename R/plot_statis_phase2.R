#' Plot STATIS Dual Robust Control Chart - Phase 2 (All Batches)
#'
#' Plots the Chi-squared statistics for both Phase 1 and Phase 2 batches,
#' using the robust Hotelling T² statistic per batch. Batches from Phase 1 are assumed
#' to be under control and those from Phase 2 are evaluated against the control threshold
#' obtained from a chi-squared distribution with df equal to the number of variables.
#'
#' @param phase1_result A list returned by `robust_statis_phase1()`, which includes
#' `batch_statistics` with Chi-squared values for Phase 1 batches.
#' @param phase2_result A list returned by `robust_statis_phase2()`, including
#' `chi_stats_by_batch` with Hotelling T² values and a control `threshold`.
#' @param title Optional string. Plot title.
#'
#' @return A ggplot2 object representing the control chart with all batches (Fase 1 y Fase 2).
#' @export
#'
#' @import ggplot2
#' @importFrom stats qchisq
#' @importFrom forcats fct_inorder
#'
#' @examples
#' data("datos_farma")
#' phase1 <- robust_statis_phase1(
#'   subset(datos_farma, Fase == "Fase 1" & Status == "Under Control"),
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density")
#' )
#' phase2 <- robust_statis_phase2(
#'   new_data = subset(datos_farma, Fase == "Fase 2"),
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density"),
#'   medians = phase1$global_medians,
#'   mads = phase1$global_mads,
#'   compromise_matrix = phase1$compromise_matrix,
#'   global_center = phase1$global_center
#' )
#' plot_statis_phase2_chart(phase1, phase2)

plot_statis_phase2_chart <- function(phase1_result, phase2_result,
                                     title = "Robust STATIS Dual Control Chart - All Batches") {
  df1 <- phase1_result$batch_statistics
  df2 <- phase2_result$chi_stats_by_batch

  df1$Status <- "Under Control"
  df2$Status <- "Out of Control"

  df1 <- df1[, c("Batch", "Chi2_Stat", "Status")]
  df2 <- df2[, c("Batch", "Chi2_Stat", "Status")]

  combined <- rbind(df1, df2)
  combined$Batch <- forcats::fct_inorder(combined$Batch)

  ggplot(combined, aes(x = Batch, y = Chi2_Stat, color = Status, group = 1)) +
    geom_point(size = 3 ) +
    geom_line(linewidth = 0.8 ) +

    # Etiquetas para lotes bajo control
    geom_text(
      data = subset(combined, Status == "Under Control"),
      aes(label = round(Chi2_Stat, 1)),
      color = "black", vjust = 2, size = 3.2, show.legend = FALSE
    ) +

    # Etiquetas para lotes fuera de control
    geom_text(
      data = subset(combined, Status == "Out of Control"),
      aes(label = round(Chi2_Stat, 1)),
      color = "#B22222", fontface = "bold", vjust = 2, size = 3.5, show.legend = FALSE
    ) +

    geom_hline(yintercept = phase2_result$threshold, linetype = "dashed", color = "red", linewidth = 0.8) +
    annotate("text",
             x = Inf, y = phase2_result$threshold,
             label = paste0("UCL = ", round(phase2_result$threshold, 1)),
             hjust = 1.1, vjust = -1.2, color = "red", size = 4) +

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

utils::globalVariables(c("Status", "Batch", "Chi2_Stat"))
