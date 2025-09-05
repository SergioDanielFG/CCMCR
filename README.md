# robustT2: Robust Hotelling-Type T2 Control Chart Based on STATIS Dual

## Overview

`robustT2` is an R package designed for robust multivariate statistical process control.  
It implements methods based on the **STATIS Dual** approach to monitor batch-based industrial processes involving multiple correlated quality variables.

The package provides:
- Construction of a robust compromise covariance matrix using **Minimum Covariance Determinant (MCD)** estimators.  
- Robust Hotelling-type T² statistics for anomaly detection.  
- Phase II monitoring using standardized Mahalanobis distances projected onto the compromise structure.  
- Visualization tools through robust biplots (GH-Biplot and HJ-Biplot) and an interactive **Shiny dashboard**.  

An internal dataset (`pharma_data`) is included for reproducibility and demonstration.

---

## Installation

You can install the development version directly from GitHub:

``r

# install.packages("devtools")
devtools::install_github("SergioDanielFG/robustT2")

## Authors

Sergio Daniel Frutos Galarza
PhD Candidate in Multivariate Statistics
University of Salamanca (USAL)

Contributors

Omar Ruiz Barzola

Dr. Purificación Galindo Villardón

## License  

MIT © 2025 Sergio Daniel Frutos Galarza
