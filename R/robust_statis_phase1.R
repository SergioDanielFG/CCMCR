#' Robust STATIS Dual - Phase 1 (Under Control Batches)
#'
#' Applies robust STATIS Dual methodology to Phase 1 data (under control batches).
#'
#' @param data A data frame containing the process data with batch information.
#' @param variables Character vector of column names to be used in the analysis.
#'
#' @return A list containing:
#' \describe{
#'   \item{compromise_matrix}{The compromise matrix (robust)}
#'   \item{global_center}{Global robust center}
#'   \item{global_medians}{Global medians per variable}
#'   \item{global_mads}{Global MADs per variable}
#'   \item{batch_statistics}{Data frame with Batch, Chi2_Stat and Weight}
#' }
#' @importFrom stats median mad mahalanobis
#' @importFrom rrcov CovMcd
#' @export
#'
#' @examples
#' data <- simulate_pharma_batches()
#' result <- robust_statis_phase1(
#'   data,
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density"))
robust_statis_phase1 <- function(data, variables) {
  batches <- unique(data$Batch)
  p <- length(variables)

  # Estandarización robusta global
  global_medians <- apply(data[, variables], 2, median)
  global_mads <- apply(data[, variables], 2, mad)

  standardized_data <- data
  for (var in variables) {
    standardized_data[[var]] <- (data[[var]] - global_medians[[var]]) / global_mads[[var]]
  }

  # Calcular medias robustas y covarianzas robustas por lote
  cov_matrices <- list()
  robust_means <- list()

  for (batch in batches) {
    subset_batch <- standardized_data[standardized_data$Batch == batch, variables]
    if (nrow(subset_batch) <= p) next
    cov_est <- rrcov::CovMcd(subset_batch)
    cov_matrices[[batch]] <- cov_est@cov
    robust_means[[batch]] <- cov_est@center
  }

  # Eliminar lotes inválidos (NULL)
  valid_batches <- names(robust_means)
  cov_matrices <- cov_matrices[valid_batches]
  robust_means <- robust_means[valid_batches]
  k_valid <- length(valid_batches)

  # Mostrar lotes usados
  message("Valid batches used in compromise construction: ", paste(valid_batches, collapse = ", "))

  # Calcular matriz de similitud (producto de Hilbert-Schmidt)
  S <- matrix(0, nrow = k_valid, ncol = k_valid)
  for (i in 1:k_valid) {
    for (j in 1:k_valid) {
      S[i, j] <- sum(cov_matrices[[i]] * cov_matrices[[j]])
    }
  }

  # Obtener pesos desde primer autovector
  eig <- eigen(S)
  weights <- eig$vectors[, 1]
  weights <- weights / sum(weights)
  names(weights) <- valid_batches

  # Compromise matrix
  compromise_matrix <- Reduce("+", Map(function(w, mat) w * mat, weights, cov_matrices))

  # Centro global robusto
  global_center <- Reduce("+", Map(function(w, mu) w * mu, weights, robust_means))

  # Estadísticos Chi-cuadrado por lote
  chi2_stats <- sapply(valid_batches, function(batch) {
    mu <- robust_means[[batch]]
    mahalanobis(mu, global_center, compromise_matrix)
  })

  # Tabla de estadísticas por lote
  batch_statistics <- data.frame(
    Batch = factor(valid_batches, levels = sort(valid_batches)),
    Chi2_Stat = chi2_stats,
    Weight = weights[valid_batches]
  )

  return(list(
    compromise_matrix = compromise_matrix,
    global_center = global_center,
    global_medians = global_medians,
    global_mads = global_mads,
    batch_statistics = batch_statistics
  ))
}
