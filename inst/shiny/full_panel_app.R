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
        tabPanel("HJ-Biplot (Fase 1)", plotOutput("biplot", height = "500px")),
        tabPanel("Projection Biplot (Fase 2)", plotOutput("projection_biplot", height = "500px")),
        tabPanel("Hotelling T² Clásico - Fase 1", plotOutput("classic_phase1", height = "500px")),
        tabPanel("Hotelling T² Clásico - Fase 2", plotOutput("classic_phase2", height = "500px"))
      )
    )
  )
)

server <- function(input, output, session) {
  vars <- c("Concentration", "Humidity", "Dissolution", "Density")

  # Datos simulados
  datos_simulados <- simulate_pharma_batches()

  # Fase 1 robusta
  phase1_data <- subset(datos_simulados, Fase == "Fase 1" & Status == "Under Control")
  phase1 <- reactive({
    robust_statis_phase1(data = phase1_data, variables = vars)
  })

  # Fase 2 robusta
  phase2_data <- subset(datos_simulados, Fase == "Fase 2")
  phase2 <- reactive({
    robust_statis_phase2(
      new_data = phase2_data,
      variables = vars,
      medians = phase1()$global_medians,
      mads = phase1()$global_mads,
      compromise_matrix = phase1()$compromise_matrix,
      global_center = phase1()$global_center
    )
  })

  # Hotelling clásico - Fase 1
  classical_f1 <- hotelling_t2_phase1(
    data = subset(datos_simulados, Fase == "Fase 1"),
    variables = vars
  )

  # Hotelling clásico - Fase 2
  classical_f2 <- hotelling_t2_phase2(
    new_data = phase2_data,
    variables = vars,
    center = classical_f1$center,
    covariance = classical_f1$covariance
  )

  # Estado de Fase 2
  status_f2 <- unique(phase2_data[, c("Batch", "Status")])
  classical_f2$batch_statistics <- merge(classical_f2$batch_statistics, status_f2, by = "Batch")

  # Gráfico robusto Fase 1
  output$phase1_plot <- renderPlot({
    plot_statis_phase1_chart(batch_statistics = phase1()$batch_statistics, num_vars = length(vars))
  })

  # Gráfico robusto Fase 2
  output$phase2_plot <- renderPlot({
    plot_statis_phase2_chart(phase2_result = phase2())
  })

  # HJ-Biplot Fase 1
  output$biplot <- renderPlot({
    plot_statis_hj_biplot(phase1_result = phase1(), color_by = input$color_by)
  })

  # Biplot de proyección Fase 2
  output$projection_biplot <- renderPlot({
    plot_statis_biplot_projection(
      phase1_result = phase1(),
      phase2_result = phase2()
    )
  })

  # Gráfico clásico T² - Fase 1
  output$classic_phase1 <- renderPlot({
    plot_classical_hotelling_t2_chart(t2_statistics = classical_f1$batch_statistics, num_vars = length(vars))
  })

  # Gráfico clásico T² - Fase 2
  output$classic_phase2 <- renderPlot({
    plot_classical_hotelling_t2_phase2_chart(t2_statistics = classical_f2$batch_statistics, num_vars = length(vars))
  })
}

app <- shinyApp(ui = ui, server = server)
