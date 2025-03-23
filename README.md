# CCMCR: Multivariate Control Chart with Chi-Square and Robust Estimators

⚠️ **Preliminary version**: This package currently includes simulation tools for pharmaceutical quality control processes. Robust multivariate control charts will be added in upcoming versions.

---

## Overview

`CCMCR` is an R package developed to support research on statistical process control in pharmaceutical manufacturing. It provides tools to simulate batches with multivariate quality variables under controlled and out-of-control conditions.

---

## Installation

You can install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("SergioDanielFG/CCMCR")
---

## Included functions

- `simulate_pharma_batches()`: Simulates pharmaceutical batches with 4 variables: Concentration, Humidity, Dissolution, and Density.
- Dataset `datos_farma`: Contains a simulated dataset with 10 batches under control and 2 batches out of control.

---

## Example

```r
library(CCMCR)

# Simulate data
datos <- simulate_pharma_batches()

# View first few rows
head(datos)

Future development

The package will include robust estimators (e.g., MCD, MVE) and the implementation of robust multivariate control charts for real pharmaceutical datasets.

Author

Sergio Daniel Frutos Galarza
Doctoral Candidate - University of Salamanca
Collaborators: Omar Ruiz Barzola, Purificación Galindo Villardón

License

MIT © 2025 Sergio Daniel Frutos Galarza
