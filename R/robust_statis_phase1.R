#' Robust STATIS Dual - Phase 1 (Under Control Batches)
#'
#' Applies robust STATIS Dual methodology to Phase 1 data (under control batches),
#' using robust standardization (median and MAD) by batch. Computes Hotelling-type T² statistics
#' per batch using robust batch centers and the compromise matrix.
#'
#' @param data A data frame containing the process data with batch information.
#' @param variables Character vector of column names to be used in the analysis.
#'
#' @return A list containing:
#' \describe{
#'   \item{compromise_matrix}{The compromise matrix (robust)}
#'   \item{global_center}{Global robust center}
#'   \item{batch_statistics}{Data frame with Batch, Chi2_Stat and Weight}
#'   \item{batch_medians}{Named list of per-batch medians per variable}
#'   \item{batch_mads}{Named list of per-batch MADs per variable}
#'   \item{global_medians}{Global medians per variable (for phase 2)}
#'   \item{global_mads}{Global MADs per variable (for phase 2)}
#' }
#' @importFrom stats median mad
#' @importFrom rrcov CovMcd
#' @export
#'
#' @examples
#' data("datos_farma")
#' phase1_data <- subset(datos_farma, Fase == "Fase 1" & Status == "Under Control")
#' result <- robust_statis_phase1(
#'   data = phase1_data,
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density")
#' )
#' result$compromise_matrix
#' result$batch_statistics
#' result$standardized_data

robust_statis_phase1 <- function(data, variables) {
  batches <- unique(data$Batch)
  p <- length(variables)

  standardized_data <- data
  batch_medians <- list()
  batch_mads <- list()

  # Estandarización robusta por lote
  for (batch in batches) {
    rows <- data$Batch == batch
    batch_data <- data[rows, variables]

    medians <- apply(batch_data, 2, median)
    mads <- apply(batch_data, 2, mad, constant = 1)

    batch_medians[[batch]] <- medians
    batch_mads[[batch]] <- mads

    for (var in variables) {
      standardized_data[rows, var] <- (data[rows, var] - medians[var]) / mads[var]
    }
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

  # Eliminar lotes inválidos
  valid_batches <- names(robust_means)
  cov_matrices <- cov_matrices[valid_batches]
  robust_means <- robust_means[valid_batches]
  k_valid <- length(valid_batches)

  message("Valid batches used in compromise construction: ", paste(valid_batches, collapse = ", "))

  # Matriz de similitud (Hilbert-Schmidt)
  S <- matrix(0, nrow = k_valid, ncol = k_valid)
  for (i in 1:k_valid) {
    for (j in 1:k_valid) {
      S[i, j] <- sum(cov_matrices[[i]] * cov_matrices[[j]])
    }
  }

  # Pesos desde el primer autovector
  eig <- eigen(S)
  weights <- eig$vectors[, 1]
  weights <- weights / sum(weights)
  names(weights) <- valid_batches

  # Matriz de compromiso
  compromise_matrix <- Reduce("+", Map(function(w, mat) w * mat, weights, cov_matrices))

  # Centro robusto global
  global_center <- Reduce("+", Map(function(w, mu) w * mu, weights, robust_means))

  # Calcular estadístico T² robusto por lote (Hotelling-type)
  chi2_stats <- data.frame(Batch = character(), Chi2_Stat = numeric(), Weight = numeric())

  for (batch in valid_batches) {
    n_b <- sum(data$Batch == batch)
    x_b <- robust_means[[batch]]
    diff <- x_b - global_center
    T2_b <- n_b * t(diff) %*% solve(compromise_matrix) %*% diff
    chi2_stats <- rbind(chi2_stats, data.frame(
      Batch = batch,
      Chi2_Stat = as.numeric(T2_b),
      Weight = weights[batch]
    ))
  }

  batch_statistics <- chi2_stats

  # Mediana y MAD global para fase 2
  global_medians <- apply(data[, variables], 2, median)
  global_mads <- apply(data[, variables], 2, mad, constant = 1)

  return(list(
    compromise_matrix = compromise_matrix,
    global_center = global_center,
    batch_statistics = batch_statistics,
    batch_medians = batch_medians,
    batch_mads = batch_mads,
    global_medians = global_medians,
    global_mads = global_mads,
    robust_means = robust_means,
    standardized_data = standardized_data,
    robust_covariances = cov_matrices
  ))
}
