#' Plot STATIS Dual Robust Control Chart - Phase 2 Only
#'
#' Plots the robust Hotelling T² statistics for Phase 2 batches only,
#' using the results from the robust STATIS Dual method.
#'
#' @param phase2_result A list returned by `robust_statis_phase2()`, including
#' `t2_stats_by_batch` with Hotelling T² values and a control `threshold`.
#' @param title Optional string. Plot title.
#'
#' @return A ggplot2 object representing the control chart for Phase 2 batches.
#' @export
#'
#' @import ggplot2
#' @importFrom forcats fct_inorder
#'
#' @examples
#' datos <- simulate_pharma_batches()
#' phase1 <- robust_statis_phase1(
#'   data = subset(datos, Fase == "Fase 1" & Status == "Under Control"),
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density")
#' )
#' phase2 <- robust_statis_phase2(
#'   new_data = subset(datos, Fase == "Fase 2"),
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density"),
#'   medians = phase1$global_medians,
#'   mads = phase1$global_mads,
#'   compromise_matrix = phase1$compromise_matrix,
#'   global_center = phase1$global_center
#' )
#' plot_statis_phase2_chart(phase2_result = phase2)

plot_statis_phase2_chart <- function(phase2_result,
                                     title = "Robust STATIS Dual Control Chart - Phase 2") {
  df2 <- phase2_result$t2_stats_by_batch

  # Verificar o construir la columna Status usando el umbral
  if (!"Status" %in% colnames(df2)) {
    df2$Status <- ifelse(df2$T2_Stat > phase2_result$threshold, "Out of Control", "Under Control")
  }


  df2 <- df2[, c("Batch", "T2_Stat", "Status")]
  df2$Batch <- forcats::fct_inorder(df2$Batch)
  df2$Status <- factor(df2$Status, levels = c("Under Control", "Out of Control"))

  ggplot(df2, aes(x = Batch, y = T2_Stat, color = Status, group = 1)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.8) +
    geom_text(
      aes(label = round(T2_Stat, 1)),
      color = "black", vjust = 2, size = 3.2, show.legend = FALSE
    ) +
    geom_hline(yintercept = phase2_result$threshold, linetype = "dashed", color = "red", linewidth = 0.8) +
    annotate("text",
             x = Inf, y = phase2_result$threshold,
             label = paste0("UCL = ", round(phase2_result$threshold, 1)),
             hjust = 1.1, vjust = -1.2, color = "red", size = 4) +
    scale_color_manual(values = c("Under Control" = "#00C8D7", "Out of Control" = "#B22222")) +
    labs(
      title = title,
      x = "Batch",
      y = expression(T^2 ~ "(Hotelling~per~Batch)"),
      color = "Status"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11)
    )
}

utils::globalVariables(c("Status", "Batch", "T2_Stat"))
