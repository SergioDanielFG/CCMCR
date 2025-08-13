#' HJ-Biplot Projection - Robust STATIS Dual (Fase 2)
#'
#' Projects new batches from Phase 2 into the HJ-Biplot space defined by the robust compromise matrix
#' and eigen decomposition from Phase 1.
#'
#' @details
#' This implementation follows the HJ-Biplot formulation of Galindo-Villardón (1986).
#' The compromise matrix \eqn{C}, being symmetric and positive semidefinite, is
#' decomposed via an eigen decomposition, not a generic Singular Value Decomposition (SVD).
#' The square roots of the eigenvalues (\eqn{\sigma_i = \sqrt{\lambda_i}}) are used
#' by convention to build the biplot scaling. This ensures consistency with the
#' robust STATIS Dual methodology and avoids confusion with standard SVD of rectangular matrices.
#'
#' @param phase1_result Result from `robust_statis_phase1()`.
#' @param phase2_result Result from `robust_statis_phase2()` (must include standardized_data and chi_stats_by_batch).
#' @param dims Dimensions to plot (default: c(1, 2)).
#'
#' @return A ggplot2 object with the projected HJ-Biplot for Phase 2 batches.
#' @export
#' @importFrom stats aggregate
#' @examples
#' datos <- simulate_pharma_batches()
#' phase1_data <- subset(datos, Fase == "Fase 1" & Status == "Under Control")
#' phase2_data <- subset(datos, Fase == "Fase 2")
#'
#' phase1 <- robust_statis_phase1(
#'   data = phase1_data,
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density")
#' )
#'
#' phase2 <- robust_statis_phase2(
#'   new_data = phase2_data,
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density"),
#'   medians = phase1$global_medians,
#'   mads = phase1$global_mads,
#'   compromise_matrix = phase1$compromise_matrix,
#'   global_center = phase1$global_center
#' )
#'
#' plot_statis_biplot_projection(phase1, phase2)

plot_statis_biplot_projection <- function(phase1_result, phase2_result, dims = c(1, 2)) {
  stopifnot(length(dims) == 2)

  # --- 1. Eigen descomposición de la matriz de compromiso ---
  compromise <- phase1_result$compromise_matrix
  eig <- eigen(compromise)
  eig_vectors <- eig$vectors
  eig_values  <- eig$values

  D_k <- diag(sqrt(eig_values[dims]))
  V_k <- eig_vectors[, dims]

  # --- 2. Calcular centros multivariados por lote (fase 2) ---
  selected_vars <- colnames(compromise)
  phase2_means <- aggregate(
    phase2_result$standardized_data[, selected_vars],
    by = list(Batch = phase2_result$standardized_data$Batch),
    FUN = mean
  )
  rownames(phase2_means) <- phase2_means$Batch
  G2 <- as.matrix(phase2_means[, selected_vars]) %*% V_k %*% D_k

  df_batches <- data.frame(G2, Batch = rownames(phase2_means))
  colnames(df_batches)[1:2] <- c("V1", "V2")

  # --- 3. Clasificación según estadístico T² robusto ---
  stats_phase2 <- phase2_result$t2_stats_by_batch
  ucl <- phase2_result$threshold
  df_batches$Status <- ifelse(stats_phase2$T2_Stat > ucl,
                              "Out of Control", "Under Control")

  # --- 4. Variables como vectores desde el origen ---
  H <- V_k %*% D_k
  df_vars <- data.frame(H, Variable = selected_vars)
  colnames(df_vars)[1:2] <- c("V1", "V2")

  # --- 5. Varianza explicada ---
  var_explained <- round(100 * eig_values[dims] / sum(eig_values), 1)

  # --- 6. Márgenes automáticos ---
  x_range <- range(c(df_batches$V1, df_vars$V1 * 1.4), na.rm = TRUE)
  y_range <- range(c(df_batches$V2, df_vars$V2 * 1.4), na.rm = TRUE)
  x_margin <- diff(x_range) * 0.15
  y_margin <- diff(y_range) * 0.15
  x_limits <- c(x_range[1] - x_margin, x_range[2] + x_margin)
  y_limits <- c(y_range[1] - y_margin, y_range[2] + y_margin)

  # --- 7. Construcción del gráfico ---
  g <- ggplot() +
    geom_hline(yintercept = 0, color = "grey40") +
    geom_vline(xintercept = 0, color = "grey40") +

    # Flechas de variables
    geom_segment(data = df_vars,
                 aes(x = 0, y = 0, xend = V1 * 1.4, yend = V2 * 1.4),
                 arrow = arrow(length = unit(0.25, "cm")),
                 color = "brown", linewidth = 1.1) +
    ggrepel::geom_text_repel(
      data = df_vars,
      aes(x = V1 * 2, y = V2 * 2, label = Variable),
      color = "brown", size = 5, fontface = "bold",
      segment.color = NA, max.overlaps = Inf
    ) +

    # Puntos de lotes
    geom_point(data = df_batches,
               aes(x = V1, y = V2, color = Status, size = Status)) +
    ggrepel::geom_text_repel(data = df_batches,
                             aes(x = V1, y = V2, label = Batch),
                             size = 4, color = "black") +

    scale_color_manual(values = c("Under Control" = "#0072B2",
                                  "Out of Control" = "brown")) +
    scale_size_manual(values = c("Under Control" = 3,
                                 "Out of Control" = 4)) +

    labs(
      title = "HJ-Biplot Projection - Phase 2 Batches",
      x = paste0("Dim ", dims[1], " (", var_explained[1], "%)"),
      y = paste0("Dim ", dims[2], " (", var_explained[2], "%)"),
      color = "Batch Status",
      size  = "Batch Status"
    ) +
    coord_cartesian(xlim = x_limits, ylim = y_limits) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 11)
    )

  return(g)
}

utils::globalVariables(c("V1", "V2", "Variable", "Batch", "Status"))
