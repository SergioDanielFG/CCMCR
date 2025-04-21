#' GH-Biplot of Robust STATIS Dual Compromise (Galindo-Hernández Biplot)
#'
#' Generates a GH-Biplot (Galindo-Hernández) using the compromise matrix obtained
#' from robust STATIS Dual (without trace normalization). It represents variables
#' as vectors (axes) and batch centers as projected points in a reduced 2D space.
#'
#' This biplot enables visual interpretation of how batch centers relate to
#' the compromise structure derived from robust covariance matrices. Variables
#' are represented as vectors from the origin, while batch centers are projected
#' using the eigenvectors of the compromise matrix.
#'
#' @param phase1_result Result from `robust_statis_phase1()`, must include `compromise_matrix`,
#'   `robust_means`, `batch_statistics`, and `global_center`.
#' @param dims Numeric vector of length 2 indicating the dimensions to plot (default: c(1, 2)).
#' @param color_by Optional string to color batches: "none" (default), "weight" (from STATIS weights), or "distance" (Chi² stat).
#' @param highlight_batches Optional: Vector of batch names to enlarge and emphasize.
#'
#' @return A ggplot2 object representing the GH-Biplot:
#'   - Arrows: variables (based on eigenvectors of the compromise matrix)
#'   - Dots: robust centers of batches (projected)
#'   - Red point: global robust center (compromise)
#' @export
#'
#' @import ggplot2
#' @importFrom grid unit
#' @importFrom ggrepel geom_text_repel
#' @importFrom forcats fct_inorder
#'
#' @examples
#' data("datos_farma")
#' phase1_result <- robust_statis_phase1(
#'   subset(datos_farma, Fase == "Fase 1" & Status == "Under Control"),
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density")
#' )
#'
#' # Biplot básico
#' plot_statis_biplot(phase1_result)
#'
#' # Biplot coloreado por pesos del STATIS
#' plot_statis_biplot(phase1_result, color_by = "weight")
#'
#' # Biplot coloreado por estadístico Chi² robusto por lote
#' plot_statis_biplot(phase1_result, color_by = "distance")
#'
#' # Biplot destacando algunos lotes
#' plot_statis_biplot(phase1_result, highlight_batches = c("Batch_1", "Batch_10"))
plot_statis_biplot <- function(phase1_result,
                               dims = c(1, 2),
                               color_by = c("none", "weight", "distance"),
                               highlight_batches = NULL) {
  stopifnot(length(dims) == 2)
  color_by <- match.arg(color_by)

  compromise <- phase1_result$compromise_matrix
  centers <- do.call(rbind, phase1_result$robust_means)
  batches <- names(phase1_result$robust_means)
  compromise_center <- phase1_result$global_center

  if (!all.equal(compromise, t(compromise), tolerance = 1e-8)) {
    warning("Compromise matrix is not symmetric. Eigen decomposition may not be reliable.")
  }

  eig <- eigen(compromise)
  eig_vectors <- eig$vectors
  eig_values <- eig$values

  variable_coords <- eig_vectors[, dims]
  rownames(variable_coords) <- colnames(compromise)

  projected_centers <- as.matrix(centers) %*% eig_vectors[, dims]
  rownames(projected_centers) <- batches

  projected_compromise <- as.numeric(compromise_center %*% eig_vectors[, dims])
  df_compromise <- data.frame(V1 = projected_compromise[1],
                              V2 = projected_compromise[2],
                              Label = "Compromise")

  df_vars <- as.data.frame(variable_coords)
  df_vars$Variable <- rownames(df_vars)

  df_batches <- as.data.frame(projected_centers)
  colnames(df_batches)[1:2] <- c("V1", "V2")
  df_batches$Batch <- rownames(df_batches)
  df_batches$Batch <- forcats::fct_inorder(df_batches$Batch)

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
    geom_hline(yintercept = 0, linetype = "solid", color = "grey30") +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey30") +

    # Flechas de variables
    geom_segment(data = df_vars,
                 aes(x = 0, y = 0, xend = V1 * 1.2, yend = V2 * 1.2),
                 arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
                 color = "brown", linewidth = 1) +

    # Etiquetas de variables
    ggrepel::geom_text_repel(data = df_vars,
                             aes(x = V1 * 1.3, y = V2 * 1.3, label = Variable),
                             color = "brown", size = 4.5, fontface = "bold") +

    # Puntos de lotes
    geom_point(data = df_batches,
               aes(x = V1, y = V2, color = Color, size = Size)) +

    # Etiquetas de lotes
    ggrepel::geom_text_repel(data = df_batches,
                             aes(x = V1, y = V2, label = Batch),
                             size = 3, color = "black") +

    # Solo punto rojo del compromiso
    geom_point(data = df_compromise,
               aes(x = V1, y = V2),
               color = "red", size = 4) +

    # Etiqueta debajo del punto rojo
    geom_text(data = df_compromise,
              aes(x = V1, y = V2, label = Label),
              color = "red", fontface = "bold", size = 4, vjust = 1.8) +

    # Escalas y estilos
    { if (color_by != "none") scale_color_viridis_c(
      name = ifelse(color_by == "weight", "STATIS Weight", "Chi^2 Distance"),
      option = "C", end = 0.9
    ) else scale_color_identity() } +
    scale_size_identity() +

    labs(
      title = "GH-Biplot - Robust STATIS Dual Compromise",
      x = paste0("Dim ", dims[1], " (", var_explained[1], "%)"),
      y = paste0("Dim ", dims[2], " (", var_explained[2], "%)")
    ) +
    coord_cartesian(expand = TRUE) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      legend.position = if (color_by == "none") "none" else "right",
      aspect.ratio = 0.5
    )

  return(g)
}

utils::globalVariables(c("V1", "V2", "Variable", "Color", "Batch", "Size", "Label"))

