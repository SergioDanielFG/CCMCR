% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plot_classical_hotelling_t2_chart.R
\name{plot_classical_hotelling_t2_chart}
\alias{plot_classical_hotelling_t2_chart}
\title{Plot Classical Hotelling T2 Control Chart}
\usage{
plot_classical_hotelling_t2_chart(
  t2_statistics,
  num_vars,
  title = "Classical Hotelling T2 Control Chart"
)
}
\arguments{
\item{t2_statistics}{A data frame with columns \code{Batch} and \code{T2_Stat}.}

\item{num_vars}{Integer. Number of variables used in the multivariate analysis (to compute the Chi² threshold).}

\item{title}{Optional string. Plot title.}
}
\value{
A ggplot2 object representing the control chart.
}
\description{
Plots the classical Hotelling T2 statistics per batch with a uniform color line.
Batches are evaluated against a control threshold obtained from
the chi-squared distribution with degrees of freedom equal to the number of variables.
}
\examples{
# Simulate pharmaceutical manufacturing batches
datos <- simulate_pharma_batches()

# Phase 1 analysis: use Phase 1 data
phase1_data <- subset(datos, Fase == "Fase 1")

# Apply classical Hotelling T2 methodology
t2_result <- hotelling_t2_phase1(
  data = phase1_data,
  variables = c("Concentration", "Humidity", "Dissolution", "Density")
)

# Plot classical Hotelling T2 control chart
plot_classical_hotelling_t2_chart(
  t2_statistics = t2_result$batch_statistics,
  num_vars = 4
)
}
