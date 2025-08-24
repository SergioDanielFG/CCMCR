#' Plot Classical Hotelling T2 Control Chart - Phase 2
#'
#' Plots the classical Hotelling T² statistics per batch for Phase 2 data,
#' using the reference mean and covariance matrix estimated from Phase 1.
#' Batches are color-coded by control status ("Under Control" = blue, "Out of Control" = red).
#'
#' @param t2_statistics A data frame with columns `Batch`, `T2_Stat`, and `Status`.
#' @param num_vars Integer. Number of variables used in the multivariate analysis (degrees of freedom for Chi²).
#' @param title Optional string. Plot title.
#'
#' @return A ggplot2 object with the Phase 2 control chart.
#' @export
#'
#' @import ggplot2
#' @importFrom stats qchisq cov
#' @importFrom forcats fct_inorder
#'
#' @examples
#' # Simulate pharmaceutical manufacturing batches
#' sim_batches <- simulate_pharma_batches()
#'
#' # Split by phase
#' phase1_data <- subset(sim_batches, Phase == "Phase 1")
#' phase2_data <- subset(sim_batches, Phase == "Phase 2")
#'
#' # Fit Phase 1 classical estimators
#' t2_phase1 <- hotelling_t2_phase1(
#'   data = phase1_data,
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density")
#' )
#'
#' # Evaluate Phase 2 batches
#' t2_phase2 <- hotelling_t2_phase2(
#'   new_data = phase2_data,
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density"),
#'   center = t2_phase1$center,
#'   covariance = t2_phase1$covariance
#' )
#'
#' # Combine with status for plotting
#' status_info <- phase2_data[!duplicated(phase2_data$Batch), "Status"]
#' t2_phase2_plot <- cbind(t2_phase2$batch_statistics, Status = status_info)
#'
#' # Plot Phase 2 control chart
#' plot_classical_hotelling_t2_phase2_chart(
#'   t2_statistics = t2_phase2_plot,
#'   num_vars = 4
#' )
plot_classical_hotelling_t2_phase2_chart <- function(t2_statistics, num_vars,
                                                     title = "Classical Hotelling T2 Control Chart (Phase 2)") {
  chi_threshold <- qchisq(0.9973, df = num_vars)
  t2_statistics$Batch <- forcats::fct_inorder(t2_statistics$Batch)

  g <- ggplot(t2_statistics, aes(x = Batch, y = T2_Stat, group = 1, color = Status)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.8) +

    # Labels for under-control batches
    geom_text(
      data = subset(t2_statistics, Status == "Under Control"),
      aes(label = round(T2_Stat, 1)),
      color = "black", vjust = 2, size = 3.2, show.legend = FALSE
    ) +

    # Labels for out-of-control batches
    geom_text(
      data = subset(t2_statistics, Status == "Out of Control"),
      aes(label = round(T2_Stat, 1)),
      color = "black", vjust = 2, size = 3.5, show.legend = FALSE
    ) +

    geom_hline(yintercept = chi_threshold, linetype = "dashed", color = "red", linewidth = 0.8) +
    annotate("text",
             x = Inf, y = chi_threshold,
             label = paste0("UCL = ", round(chi_threshold, 1)),
             hjust = 1.1, vjust = -1.2, color = "red", size = 4) +

    scale_color_manual(values = c("Under Control" = "#00C8D7", "Out of Control" = "#B22222")) +

    labs(
      title = title,
      x = "Batch",
      y = expression(T^2 ~ "(Classical)"),
      color = "Status"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11)
    )

  return(g)
}

utils::globalVariables(c("T2_Stat", "Batch", "Status"))
