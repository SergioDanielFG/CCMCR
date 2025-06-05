#' STATIS Dual Biplot Projection for Phase 2 Batches
#'
#' Projects Phase 2 batch centers onto the reduced space defined by the compromise matrix
#' obtained in Phase 1 (from `robust_statis_phase1()`).
#' Only Phase 2 batches are displayed, enabling visual detection of shifts or anomalies
#' relative to the robust reference structure.
#'
#' This is not a GH-Biplot in sentido estricto, ya que no se recalcula una nueva base.
#' En su lugar, utiliza la base de la matriz de compromiso obtenida en la Fase 1 y proyecta
#' los nuevos lotes sobre ella, manteniendo la estructura estadística robusta original.
#'
#' @param phase1_result A list returned by `robust_statis_phase1()`, which includes
#' `compromise_matrix` and `global_center`.
#' @param phase2_result A list returned by `robust_statis_phase2()`, which includes
#' `standardized_data` (Phase 2 batches already standardized).
#' @param dims Numeric vector of length 2 indicating the dimensions to plot (default: c(1, 2)).
#'
#' @return A ggplot2 object showing Phase 2 batch centers projected onto the STATIS space.
#' @export
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom forcats fct_inorder
#' @importFrom grid unit
#'
#' @examples
#' # Simulate data
#' datos <- simulate_pharma_batches()
#'
#' # Apply Phase 1 (STATIS Dual robust)
#' phase1 <- robust_statis_phase1(
#'   data = subset(datos, Fase == "Fase 1" & Status == "Under Control"),
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density")
#' )
#'
#' # Apply Phase 2
#' phase2 <- robust_statis_phase2(
#'   new_data = subset(datos, Fase == "Fase 2"),
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density"),
#'   medians = phase1$global_medians,
#'   mads = phase1$global_mads,
#'   compromise_matrix = phase1$compromise_matrix,
#'   global_center = phase1$global_center
#' )
#'
#' # Plot projection biplot of Phase 2 batches
#' plot_statis_biplot_projection(phase1_result = phase1, phase2_result = phase2)

plot_statis_biplot_projection <- function(phase1_result,
                                          phase2_result,
                                          dims = c(1, 2)) {
  stopifnot(length(dims) == 2)

  compromise <- phase1_result$compromise_matrix
  compromise_center <- phase1_result$global_center
  eig <- eigen(compromise)
  eig_vectors <- eig$vectors
  eig_values <- eig$values
  var_explained <- round(100 * eig_values[dims] / sum(eig_values), 1)

  # Project Phase 2 batch centers
  phase2_data <- phase2_result$standardized_data
  variables <- colnames(compromise)

  centers_phase2 <- do.call(rbind, lapply(split(phase2_data, phase2_data$Batch), function(df) {
    colMeans(df[, variables, drop = FALSE])
  }))

  projected_phase2 <- centers_phase2 %*% eig_vectors[, dims]
  df_phase2 <- as.data.frame(projected_phase2)
  colnames(df_phase2) <- c("V1", "V2")
  df_phase2$Batch <- rownames(df_phase2)
  df_phase2$Batch <- forcats::fct_inorder(df_phase2$Batch)

  # Plot
  g <- ggplot(df_phase2, aes(x = V1, y = V2)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey30") +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey30") +
    geom_point(color = "#B22222", size = 3.5) +
    ggrepel::geom_text_repel(aes(label = Batch),
                             size = 3.5, color = "black") +
    labs(
      title = "Projection of Phase 2 Batch Centers onto STATIS Compromise Space",
      x = paste0("Dim ", dims[1], " (", var_explained[1], "%)"),
      y = paste0("Dim ", dims[2], " (", var_explained[2], "%)")
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11)
    )

  return(g)
}

utils::globalVariables(c("V1", "V2", "Batch"))
