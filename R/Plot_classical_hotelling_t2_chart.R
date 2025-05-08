#' Plot Classical Hotelling T2 Control Chart
#'
#' Plots the classical Hotelling T2 statistics per batch with a uniform color line.
#' Batches are evaluated against a control threshold obtained from
#' the chi-squared distribution with degrees of freedom equal to the number of variables.
#'
#' @param t2_statistics A data frame with columns `Batch` and `T2_Stat`.
#' @param num_vars Integer. Number of variables used in the multivariate analysis (to compute the Chi² threshold).
#' @param title Optional string. Plot title.
#'
#' @return A ggplot2 object representing the control chart.
#' @export
#'
#' @import ggplot2
#' @importFrom stats qchisq
#' @importFrom forcats fct_inorder
#'
#' @examples
#' # Simulate pharmaceutical manufacturing batches
#' datos <- simulate_pharma_batches()
#'
#' # Phase 1 analysis: use Phase 1 data
#' phase1_data <- subset(datos, Fase == "Fase 1")
#'
#' # Apply classical Hotelling T2 methodology
#' t2_result <- hotelling_t2_phase1(
#'   data = phase1_data,
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density")
#' )
#'
#' # Plot classical Hotelling T2 control chart
#' plot_classical_hotelling_t2_chart(
#'   t2_statistics = t2_result$batch_statistics,
#'   num_vars = 4
#' )

plot_classical_hotelling_t2_chart <- function(t2_statistics, num_vars,
                                              title = "Classical Hotelling T2 Control Chart") {
  chi_threshold <- qchisq(0.9973, df = num_vars)

  t2_statistics$Batch <- forcats::fct_inorder(t2_statistics$Batch)

  g <- ggplot(t2_statistics, aes(x = Batch, y = T2_Stat, group = 1)) +
    geom_point(size = 3, color = "#00C8D7") +
    geom_line(linewidth = 0.8, color = "#00C8D7") +

    # Labels for batches
    geom_text(
      aes(label = round(T2_Stat, 1)),
      vjust = -0.8, size = 3.2, show.legend = FALSE
    ) +

    geom_hline(yintercept = chi_threshold, linetype = "dashed", color = "red", linewidth = 0.8) +
    annotate("text",
             x = Inf, y = chi_threshold,
             label = paste0("UCL = ", round(chi_threshold, 1)),
             hjust = 1.1, vjust = -0.5, color = "red", size = 4) +

    labs(
      title = title,
      x = "Batch",
      y = expression(T^2 ~ "(per~Batch)")
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11)
    )

  return(g)
}

utils::globalVariables(c("T2_Stat", "Batch"))
