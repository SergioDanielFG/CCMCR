# CCMCR: Chi-Square Control Chart with Multivariate and Robust Estimators

⚠️ **Versión preliminar**: Este paquete incluye herramientas de simulación para procesos de control de calidad en la industria farmacéutica. Las funciones para gráficos de control multivariantes robustos están en desarrollo y serán añadidas en versiones futuras.

---

## Descripción general

`CCMCR` es un paquete en R diseñado para apoyar investigaciones en control estadístico de procesos multivariados. Proporciona herramientas para simular lotes de producción con múltiples variables de calidad, tanto en condiciones bajo control como fuera de control, con posibles contaminaciones o alteraciones.

El paquete servirá como base para comparar gráficos de control clásicos (Hotelling T²) y gráficos robustos basados en estimadores como **MCD** y técnicas de tres vías como **STATIS Dual**.

---

## Instalación

Puedes instalar la versión de desarrollo directamente desde GitHub:

```r
# install.packages("devtools")
devtools::install_github("SergioDanielFG/CCMCR")
```

---

## Funciones incluidas

- `simulate_pharma_batches()`: Simula lotes farmacéuticos con 4 variables cuantitativas:
  - `Concentration`
  - `Humidity`
  - `Dissolution`
  - `Density`

- Dataset `datos_farma`: Contiene un conjunto de datos simulados con 10 lotes bajo control (Fase 1) y 5 lotes (Fase 2) que incluyen tanto lotes bajo control como fuera de control.

---

## Ejemplo básico

```r
library(CCMCR)

# Simular datos
datos <- simulate_pharma_batches()

# Ver primeras filas
head(datos)
```

---

## Desarrollo futuro

El paquete incorporará:

- Estimadores robustos: **MCD**, **MVE**
- Gráficos de control multivariantes robustos tipo Hotelling T²
- Visualizaciones interactivas (Shiny Dashboard)
- Métodos basados en STATIS Dual para construir matrices de compromiso robustas

---

## Autor

**Sergio Daniel Frutos Galarza**  
Doctorando en Estadística Multivariante  
Universidad de Salamanca (USAL)

**Colaboradores**  
- Omar Ruiz Barzola  
- Dra. Purificación Galindo Villardón

---

## Licencia

MIT © 2025 Sergio Daniel Frutos Galarza
