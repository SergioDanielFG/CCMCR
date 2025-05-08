# app.R - STATIS Dual Robust Control Dashboard (FUNCIONAL)

library(shiny)
library(ggplot2)
library(shinythemes)
library(CCMCR)

ui <- fluidPage(
  theme = shinytheme("cosmo"),
  titlePanel("STATIS Dual Robust Control Dashboard"),

  sidebarLayout(
    sidebarPanel(
      selectInput("color_by", "\u25B6 Color by:",
                  choices = c("none", "weight", "distance"),
                  selected = "none")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Robust Phase 1", plotOutput("phase1_plot", height = "500px")),
        tabPanel("Robust Phase 2", plotOutput("phase2_plot", height = "500px")),
        tabPanel("GH-Biplot", plotOutput("biplot", height = "500px")),
        tabPanel("Classical Hotelling T²", plotOutput("classic_plot", height = "500px"))
      )
    )
  )
)

server <- function(input, output, session) {
  vars <- c("Concentration", "Humidity", "Dissolution", "Density")

  # Usar simulación directamente
  datos_simulados <- simulate_pharma_batches()

  # Robust Phase 1 (solo lotes bajo control de Fase 1)
  phase1_data <- subset(datos_simulados, Fase == "Fase 1" & Status == "Under Control")
  phase1 <- robust_statis_phase1(data = phase1_data, variables = vars)

  # Robust Phase 2 (lotes nuevos de Fase 2)
  phase2_data <- subset(datos_simulados, Fase == "Fase 2")
  phase2 <- robust_statis_phase2(
    new_data = phase2_data,
    variables = vars,
    medians = phase1$global_medians,
    mads = phase1$global_mads,
    compromise_matrix = phase1$compromise_matrix,
    global_center = phase1$global_center
  )

  # Classical Hotelling T² (con todos los lotes de Fase 1, contaminados)
  classical <- hotelling_t2_phase1(data = subset(datos_simulados, Fase == "Fase 1"), variables = vars)

  # Gráfico robusto Fase 1
  output$phase1_plot <- renderPlot({
    plot_statis_phase1_chart(
      batch_statistics = phase1$batch_statistics,
      num_vars = length(vars)
    )
  })

  # Gráfico robusto Fase 2 (Fase 1 + Fase 2)
  output$phase2_plot <- renderPlot({
    plot_statis_phase2_chart(
      phase1_result = phase1,
      phase2_result = phase2
    )
  })

  # GH-Biplot
  output$biplot <- renderPlot({
    plot_statis_biplot(
      phase1_result = phase1,
      color_by = input$color_by
    )
  })

  # Gráfico clásico T²
  output$classic_plot <- renderPlot({
    plot_classical_hotelling_t2_chart(
      t2_statistics = classical$batch_statistics,
      num_vars = length(vars)
    )
  })
}

# Lanzar aplicación
app <- shinyApp(ui = ui, server = server)

