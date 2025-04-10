#' Biplot of Robust STATIS Dual Compromise
#'
#' Generates a biplot using the compromise matrix from STATIS Dual analysis. Projects the robust batch centers and variable axes onto the principal dimensions of the compromise.
#'
#' @param phase1_result Result from `robust_statis_phase1()`, must include `compromise_matrix` and `robust_means`.
#' @param dims Integer vector of length 2 indicating the dimensions to plot (default: c(1, 2)).
#'
#' @return A ggplot2 object representing the biplot.
#' @export
#'
#' @import ggplot2
#' @importFrom grid unit
#'
#' @examples
#' data <- simulate_pharma_batches()
#' phase1 <- robust_statis_phase1(
#'   data[data$Status == "Under Control", ],
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density"))
#' plot_statis_biplot(phase1)
plot_statis_biplot <- function(phase1_result, dims = c(1, 2)) {
  stopifnot(length(dims) == 2)

  compromise <- phase1_result$compromise_matrix
  centers <- do.call(rbind, phase1_result$robust_means)
  batches <- names(phase1_result$robust_means)

  # Eigen-decomposition
  eig <- eigen(compromise)
  eig_vectors <- eig$vectors
  eig_values <- eig$values

  # Projections
  variable_coords <- eig_vectors[, dims]
  rownames(variable_coords) <- colnames(compromise)

  projected_centers <- as.matrix(centers) %*% eig_vectors[, dims]
  rownames(projected_centers) <- batches

  df_vars <- as.data.frame(variable_coords)
  df_vars$Variable <- rownames(df_vars)

  df_batches <- as.data.frame(projected_centers)
  df_batches$Batch <- rownames(df_batches)

  # Percent variance explained
  var_explained <- round(100 * eig_values[dims] / sum(eig_values), 1)

  ggplot() +
    # Variable vectors
    geom_segment(data = df_vars,
                 aes(x = 0, y = 0, xend = V1, yend = V2),
                 arrow = arrow(length = unit(0.2, "cm")),
                 color = "darkred", linewidth = 0.8) +
    geom_text(data = df_vars, aes(x = V1, y = V2, label = Variable),
              color = "darkred", size = 4, hjust = 1.1) +

    # Batch centers
    geom_point(data = df_batches, aes(x = V1, y = V2), color = "#0072B2", size = 2.5) +
    geom_text(data = df_batches, aes(x = V1, y = V2, label = Batch),
              hjust = -0.2, size = 3, color = "#0072B2") +

    labs(title = "Biplot - Robust STATIS Dual Compromise",
         x = paste0("Dim ", dims[1], " (", var_explained[1], "%)"),
         y = paste0("Dim ", dims[2], " (", var_explained[2], "%)")) +
    theme_minimal(base_size = 13)
}
utils::globalVariables(c("V1", "V2", "Variable"))
