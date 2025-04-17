library(shiny)
library(ggplot2)
library(CCMCR)
library(shinythemes)

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
        tabPanel("Phase 1", plotOutput("phase1_plot", height = "500px")),
        tabPanel("Phase 2", plotOutput("phase2_plot", height = "500px")),
        tabPanel("Biplot", plotOutput("biplot", height = "500px"))
      )
    )
  )
)

server <- function(input, output, session) {
  # Cargar datos
  data("datos_farma", package = "CCMCR")
  vars <- c("Concentration", "Humidity", "Dissolution", "Density")

  # Fase 1: entrenamiento con los 10 lotes bajo control
  phase1_data <- subset(datos_farma, Fase == "Fase 1" & Status == "Under Control")
  phase1 <- CCMCR::robust_statis_phase1(phase1_data, variables = vars)

  # Fase 2: evaluación de los 7 lotes nuevos (Fase 2)
  phase2_data <- subset(datos_farma, Fase == "Fase 2")
  phase2 <- CCMCR::robust_statis_phase2(
    new_data = phase2_data,
    variables = vars,
    medians = phase1$global_medians,
    mads = phase1$global_mads,
    compromise_matrix = phase1$compromise_matrix,
    global_center = phase1$global_center
  )

  # Salida: Gráfico de control Fase 1
  output$phase1_plot <- renderPlot({
    CCMCR::plot_statis_phase1_chart(
      batch_statistics = phase1$batch_statistics,
      num_vars = length(vars)
    )
  })

  # Salida: Gráfico de control Fase 2
  output$phase2_plot <- renderPlot({
    CCMCR::plot_statis_phase2_chart(
      phase1_result = phase1,
      phase2_result = phase2
    )
  })

  # Salida: GH-Biplot interactivo
  output$biplot <- renderPlot({
    CCMCR::plot_statis_biplot(
      phase1_result = phase1,
      color_by = input$color_by
    )
  })
}

# Exportar la app
app <- shinyApp(ui = ui, server = server)
