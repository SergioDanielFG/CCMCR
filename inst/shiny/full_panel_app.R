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
  data("datos_farma")
  vars <- c("Concentration", "Humidity", "Dissolution", "Density")

  # Fase 1
  phase1 <- CCMCR::robust_statis_phase1(
    subset(datos_farma, Status == "Under Control"),
    variables = vars
  )

  # Fase 2
  phase2 <- CCMCR::robust_statis_phase2(
    new_data = subset(datos_farma, Status == "Out of Control"),
    variables = vars,
    medians = phase1$global_medians,
    mads = phase1$global_mads,
    compromise_matrix = phase1$compromise_matrix,
    global_center = phase1$global_center
  )

  # Gráfico de control - Fase 1 (usando Chi2_Stat directamente)
  output$phase1_plot <- renderPlot({
    CCMCR::plot_statis_phase1_chart(
      batch_statistics = phase1$batch_statistics,
      num_vars = length(vars)
    )
  })

  # Gráfico de control - Fase 2
  output$phase2_plot <- renderPlot({
    CCMCR::plot_statis_phase2_chart(
      phase1_result = phase1,
      phase2_result = phase2
    )
  })

  # Biplot con opción de color dinámico
  output$biplot <- renderPlot({
    CCMCR::plot_statis_biplot(
      phase1_result = phase1,
      color_by = input$color_by
    )
  })
}

# Define la app como objeto exportable
app <- shinyApp(ui = ui, server = server)
