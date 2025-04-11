#' GH-Biplot of Robust STATIS Dual Compromise (Galindo-Hernández Biplot)
#'
#' Generates a GH-Biplot (Galindo-Hernández) using the compromise matrix obtained
#' from robust STATIS Dual. It represents variables as vectors and batch centers as points.
#'
#' This is the most common form of biplot used in multivariate analysis (PCA, STATIS).
#'
#' @param phase1_result Result from `robust_statis_phase1()`, must include `compromise_matrix`,
#'   `robust_means`, and `batch_statistics`.
#' @param dims Dimensions to plot (default: c(1, 2)).
#' @param color_by Optional: "none" (default), "weight" (from STATIS), or "distance".
#' @param highlight_batches Optional: Vector of batch names to highlight.
#'
#' @return A ggplot2 object representing the GH-Biplot.
#' @export
#'
#' @import ggplot2
#' @importFrom grid unit
#'
#' @examples
#' data("datos_farma")
#' phase1 <- robust_statis_phase1(
#'   subset(datos_farma, Status == "Under Control"),
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density"))
#'
#' # Basic GH-Biplot
#' plot_statis_biplot(phase1)
#'
#' # Colored by STATIS weights
#' plot_statis_biplot(phase1, color_by = "weight")
#'
#' # Colored by Chi² distance
#' plot_statis_biplot(phase1, color_by = "distance")
#'
#' # Highlight specific batches
#' plot_statis_biplot(phase1, highlight_batches = c("Batch_7", "Batch_10"))
plot_statis_biplot <- function(phase1_result,
                               dims = c(1, 2),
                               color_by = c("none", "weight", "distance"),
                               highlight_batches = NULL) {
  stopifnot(length(dims) == 2)
  color_by <- match.arg(color_by)

  compromise <- phase1_result$compromise_matrix
  centers <- do.call(rbind, phase1_result$robust_means)
  batches <- names(phase1_result$robust_means)

  eig <- eigen(compromise)
  eig_vectors <- eig$vectors
  eig_values <- eig$values

  variable_coords <- eig_vectors[, dims]
  rownames(variable_coords) <- colnames(compromise)

  projected_centers <- as.matrix(centers) %*% eig_vectors[, dims]
  rownames(projected_centers) <- batches

  df_vars <- as.data.frame(variable_coords)
  df_vars$Variable <- rownames(df_vars)

  df_batches <- as.data.frame(projected_centers)
  df_batches$Batch <- rownames(df_batches)

  var_explained <- round(100 * eig_values[dims] / sum(eig_values), 1)

  if (color_by == "weight") {
    weights <- phase1_result$batch_statistics$Weight
    names(weights) <- phase1_result$batch_statistics$Batch
    df_batches$Color <- weights[df_batches$Batch]
  } else if (color_by == "distance") {
    dists <- phase1_result$batch_statistics$Chi2_Stat
    names(dists) <- phase1_result$batch_statistics$Batch
    df_batches$Color <- dists[df_batches$Batch]
  } else {
    df_batches$Color <- "#0072B2"
  }

  df_batches$Size <- ifelse(df_batches$Batch %in% highlight_batches, 4, 2.5)

  g <- ggplot() +
    geom_segment(data = df_vars,
                 aes(x = 0, y = 0, xend = V1, yend = V2),
                 arrow = arrow(length = unit(0.2, "cm")),
                 color = "darkred", linewidth = 0.8) +
    geom_text(data = df_vars,
              aes(x = V1, y = V2, label = Variable),
              color = "darkred", size = 4, hjust = 1.1) +
    geom_point(data = df_batches,
               aes(x = V1, y = V2, color = Color, size = Size)) +
    geom_text(data = df_batches,
              aes(x = V1, y = V2, label = Batch),
              hjust = -0.2, size = 3, color = "black") +
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +
    geom_vline(xintercept = 0, linetype = "solid", color = "black") +
    labs(
      title = "GH-Biplot - Robust STATIS Dual Compromise",
      x = paste0("Dim ", dims[1], " (", var_explained[1], "%)"),
      y = paste0("Dim ", dims[2], " (", var_explained[2], "%)")
    ) +
    coord_equal() +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      legend.position = if (color_by == "none") "none" else "right"
    ) +
    scale_size_identity()



  if (color_by != "none") {
    g <- g + scale_color_viridis_c(name = color_by)
  }

  return(g)
}

utils::globalVariables(c("V1", "V2", "Variable", "Color", "Batch", "Size"))
