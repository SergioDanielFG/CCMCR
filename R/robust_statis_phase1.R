#' Robust STATIS Dual - Phase 1 (Under Control Batches)
#'
#' Applies the Robust STATIS Dual methodology to Phase 1 data (under control batches),
#' using robust batch-wise standardization (median and MAD scaled with constant 1.4826).
#' Covariance matrices are robustly estimated using the MCD method
#' and used directly (without trace normalization) to construct the compromise matrix.
#'
#' @param data A data frame containing the process data with batch information.
#' @param variables Character vector with the names of the variables to be used in the analysis.
#'
#' @return A list containing:
#' \describe{
#'   \item{compromise_matrix}{Robust compromise matrix (without trace normalization)}
#'   \item{global_center}{Global robust center of the batches}
#'   \item{batch_statistics}{Data frame with Batch, Chi2_Stat (Hotelling-type statistic), and Weight}
#'   \item{batch_medians}{List of medians per batch and variable}
#'   \item{batch_mads}{List of scaled MADs per batch and variable (constant 1.4826)}
#'   \item{global_medians}{Global medians per variable (for use in Phase 2)}
#'   \item{global_mads}{Global scaled MADs per variable (constant 1.4826)}
#'   \item{robust_means}{List of robust centers of each batch (estimated by MCD)}
#'   \item{standardized_data}{Data set standardized batch by batch}
#'   \item{robust_covariances}{List of robust covariance matrices per batch}
#'   \item{similarity_matrix}{Hilbert-Schmidt similarity matrix between batches}
#'   \item{statis_weights}{Weights obtained from the first eigenvector of the similarity matrix}
#'   \item{first_eigenvector}{First eigenvector of the similarity matrix (unnormalized)}
#' }
#'
#' @importFrom stats median mad
#' @importFrom rrcov CovMcd
#' @export
#'
#' @examples
#' # Simulate new pharmaceutical manufacturing batches
#' datos <- simulate_pharma_batches()
#'
#' # Select only Phase 1 under control batches
#' phase1_data <- subset(datos, Fase == "Fase 1" & Status == "Under Control")
#'
#' # Apply robust STATIS Dual methodology
#' result <- robust_statis_phase1(
#'   data = phase1_data,
#'   variables = c("Concentration", "Humidity", "Dissolution", "Density")
#' )
#'
#' # View main outputs
#' result$compromise_matrix
#' result$batch_statistics
#' result$robust_covariances
#' result$similarity_matrix
#' result$statis_weights
#' result$robust_means

robust_statis_phase1 <- function(data, variables) {
  batches <- unique(data$Batch)
  p <- length(variables)

  standardized_data <- data
  batch_medians <- list()
  batch_mads <- list()

  # Robust standardization by batch: median and MAD (scaled with constant 1.4826)
  for (batch in batches) {
    rows <- data$Batch == batch
    batch_data <- data[rows, variables]

    medians <- apply(batch_data, 2, median)
    mads <- apply(batch_data, 2, mad)  # Default scaling: constant = 1.4826

    batch_medians[[batch]] <- medians
    batch_mads[[batch]] <- mads

    for (var in variables) {
      standardized_data[rows, var] <- (data[rows, var] - medians[var]) / mads[var]
    }
  }

  # Computation of robust centers and covariance matrices per batch
  cov_matrices <- list()
  robust_means <- list()

  for (batch in batches) {
    subset_batch <- standardized_data[standardized_data$Batch == batch, variables]
    if (nrow(subset_batch) <= p) next
    cov_est <- rrcov::CovMcd(subset_batch)

    cov_matrices[[batch]] <- cov_est@cov
    robust_means[[batch]] <- cov_est@center
  }

  # Only keep valid batches
  valid_batches <- names(robust_means)
  cov_matrices <- cov_matrices[valid_batches]
  robust_means <- robust_means[valid_batches]
  k_valid <- length(valid_batches)

  message("Valid batches used for compromise construction: ", paste(valid_batches, collapse = ", "))

  # Construction of the Hilbert-Schmidt similarity matrix
  S <- matrix(0, nrow = k_valid, ncol = k_valid)
  for (i in seq_len(k_valid)) {
    for (j in seq_len(k_valid)) {
      S[i, j] <- sum(cov_matrices[[i]] * cov_matrices[[j]])
    }
  }

  # Eigen decomposition
  eig <- eigen(S)
  first_eigenvector <- eig$vectors[, 1]
  weights <- first_eigenvector / sum(first_eigenvector)
  names(first_eigenvector) <- valid_batches
  names(weights) <- valid_batches

  # Construction of the compromise matrix and the global center
  compromise_matrix <- Reduce("+", Map(function(w, mat) w * mat, weights, cov_matrices))
  global_center <- Reduce("+", Map(function(w, mu) w * mu, weights, robust_means))

  # Computation of Hotelling-type Chi² statistics per batch
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

  # Computation of global medians and MADs for Phase 2
  global_medians <- apply(data[, variables], 2, median)
  global_mads <- apply(data[, variables], 2, mad)  # Scaling: constant = 1.4826

  return(list(
    compromise_matrix = compromise_matrix,
    global_center = global_center,
    batch_statistics = chi2_stats,
    batch_medians = batch_medians,
    batch_mads = batch_mads,
    global_medians = global_medians,
    global_mads = global_mads,
    robust_means = robust_means,
    standardized_data = standardized_data,
    similarity_matrix = S,
    robust_covariances = cov_matrices,
    statis_weights = weights,
    first_eigenvector = first_eigenvector
  ))
}
