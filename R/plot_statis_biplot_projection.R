#' HJ-Biplot Projection - Robust STATIS Dual (Fase 2)
#'
#' Projects new batches from Phase 2 into the HJ-Biplot space defined by the robust compromise matrix
#' and eigen decomposition from Phase 1.
#'
#' @param phase1_result Result from `robust_statis_phase1()`.
#' @param phase2_result Result from `robust_statis_phase2()` (must include standardized_data and chi_stats_by_batch).
#' @param dims Dimensions to plot (default: c(1, 2)).
#'
#' @return A ggplot2 object with the projected HJ-Biplot for Phase 2 batches.
#' @export
#'
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

  # 1. Matriz de compromiso y descomposición espectral
  compromise <- phase1_result$compromise_matrix
  eig <- eigen(compromise)
  eig_vectors <- eig$vectors
  eig_values <- eig$values

  D_k <- diag(sqrt(eig_values[dims]))
  V_k <- eig_vectors[, dims]

  # 2. Calcular centros multivariados por lote usando las variables del compromiso
  selected_vars <- colnames(compromise)
  phase2_means <- aggregate(
    phase2_result$standardized_data[, selected_vars],
    by = list(Batch = phase2_result$standardized_data$Batch),
    FUN = mean
  )
  rownames(phase2_means) <- phase2_means$Batch
  G2 <- as.matrix(phase2_means[, selected_vars]) %*% V_k %*% D_k

  # 3. Variables como vectores desde el origen
  H <- V_k %*% D_k
  rownames(H) <- selected_vars

  df_vars <- as.data.frame(H)
  df_vars$Variable <- rownames(df_vars)

  df_batches <- as.data.frame(G2)
  colnames(df_batches) <- c("V1", "V2")
  df_batches$Batch <- forcats::fct_inorder(rownames(df_batches))

  # 4. Clasificación según distancia Chi² robusta
  stats_phase2 <- phase2_result$chi_stats_by_batch
  ucl <- phase2_result$threshold

  df_batches$Status <- ifelse(
    stats_phase2$Chi2_Stat > ucl,
    "Out of Control", "Under Control"
  )

  df_batches$Color <- ifelse(df_batches$Status == "Out of Control", "brown", "#0072B2")
  df_batches$Size <- ifelse(df_batches$Status == "Out of Control", 4, 3)

  var_explained <- round(100 * eig_values[dims] / sum(eig_values), 1)

  # 5. Márgenes automáticos
  x_range <- range(df_batches$V1, na.rm = TRUE)
  y_range <- range(df_batches$V2, na.rm = TRUE)
  x_margin <- diff(x_range) * 0.15
  y_margin <- diff(y_range) * 0.15
  x_limits <- c(x_range[1] - x_margin, x_range[2] + x_margin)
  y_limits <- c(y_range[1] - y_margin, y_range[2] + y_margin)

  # 6. Construcción del gráfico
  g <- ggplot() +
    geom_hline(yintercept = 0, color = "grey40") +
    geom_vline(xintercept = 0, color = "grey40") +
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
    geom_point(data = df_batches,
               aes(x = V1, y = V2, color = Color, size = Size)) +
    ggrepel::geom_text_repel(data = df_batches,
                             aes(x = V1, y = V2, label = Batch),
                             size = 4, color = "black", fontface = "plain") +
    scale_color_identity() +
    scale_size_identity() +
    labs(
      title = "HJ-Biplot Projection - Phase 2 Batches",
      x = paste0("Dim ", dims[1], " (", var_explained[1], "%)"),
      y = paste0("Dim ", dims[2], " (", var_explained[2], "%)")
    ) +
    coord_cartesian(xlim = x_limits, ylim = y_limits) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      legend.position = "none"
    )

  return(g)
}

utils::globalVariables(c("V1", "V2", "Variable", "Color", "Batch", "Size", "Status"))
