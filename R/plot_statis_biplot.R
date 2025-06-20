#' HJ-Biplot of Robust STATIS Dual Compromise (Galindo-Villardón)
#'
#' Generates an HJ-Biplot using the compromise matrix obtained
#' from robust STATIS Dual. Individuals (batch centers) are projected as G = U D,
#' and variables as H = V D, where D is the diagonal matrix of square roots of eigenvalues.
#'
#' @param phase1_result Result from `robust_statis_phase1()`.
#' @param dims Dimensions to plot (default: c(1, 2)).
#' @param color_by One of "none", "weight", or "distance" for coloring batches.
#' @param highlight_batches Optional vector of batch names to emphasize.
#'
#' @return ggplot2 object with HJ-Biplot.
#' @export
#' @importFrom stats setNames

#' @import ggplot2
#' @importFrom grid unit
#' @importFrom ggrepel geom_text_repel
#' @importFrom forcats fct_inorder
#'
#' @examples
#' datos <- simulate_pharma_batches()
#' phase1 <- robust_statis_phase1(
#'   data = subset(datos, Fase == "Fase 1" & Status == "Under Control"),
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density")
#' )
#' plot_statis_hj_biplot(phase1)

plot_statis_hj_biplot <- function(phase1_result,
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

  # HJ-Biplot scaling: G = U D and H = V D
  D_k <- diag(sqrt(eig_values[dims]))
  V_k <- eig_vectors[, dims]

  G <- as.matrix(centers) %*% V_k %*% D_k  # batch centers
  H <- V_k %*% D_k                         # variable vectors

  rownames(G) <- batches
  rownames(H) <- colnames(compromise)

  df_vars <- as.data.frame(H)
  df_vars$Variable <- rownames(df_vars)

  df_batches <- as.data.frame(G)
  colnames(df_batches) <- c("V1", "V2")
  df_batches$Batch <- forcats::fct_inorder(rownames(df_batches))

  var_explained <- round(100 * eig_values[dims] / sum(eig_values), 1)

  if (color_by == "weight") {
    weights <- setNames(phase1_result$batch_statistics$Weight,
                        phase1_result$batch_statistics$Batch)
    df_batches$Color <- weights[as.character(df_batches$Batch)]
  } else if (color_by == "distance") {
    dists <- setNames(phase1_result$batch_statistics$Chi2_Stat,
                      phase1_result$batch_statistics$Batch)
    df_batches$Color <- dists[as.character(df_batches$Batch)]
  } else {
    df_batches$Color <- "#0072B2"
  }

  df_batches$Size <- ifelse(df_batches$Batch %in% highlight_batches, 4, 2.5)

  g <- ggplot() +
    geom_hline(yintercept = 0, color = "grey30") +
    geom_vline(xintercept = 0, color = "grey30") +
    geom_segment(data = df_vars,
                 aes(x = 0, y = 0, xend = V1 * 1.2, yend = V2 * 1.2),
                 arrow = arrow(length = unit(0.25, "cm")), color = "brown", linewidth = 1) +
    geom_text(data = df_vars,
              aes(x = V1 * 1.3, y = V2 * 1.3, label = Variable),
              color = "brown", size = 4.5, fontface = "bold") +
    geom_point(data = df_batches,
               aes(x = V1, y = V2, color = Color, size = Size)) +
    ggrepel::geom_text_repel(data = df_batches,
                             aes(x = V1, y = V2, label = Batch),
                             size = 3, color = "black") +
    { if (color_by != "none") scale_color_viridis_c(
      name = ifelse(color_by == "weight", "STATIS Weight", "Chi^2 Distance"),
      option = "C", end = 0.9
    ) else scale_color_identity() } +
    scale_size_identity() +
    labs(
      title = "HJ-Biplot - Robust STATIS Dual Compromise",
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

utils::globalVariables(c("V1", "V2", "Variable", "Color", "Batch", "Size"))
