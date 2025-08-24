library(shiny)
library(ggplot2)
library(shinythemes)
library(robustT2)

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
        tabPanel("HJ-Biplot (Phase 1)", plotOutput("biplot", height = "500px")),
        tabPanel("Projection Biplot (Phase 2)", plotOutput("projection_biplot", height = "500px")),
        tabPanel("Hotelling T² Classic - Phase 1", plotOutput("classic_phase1", height = "500px")),
        tabPanel("Hotelling T² Classic - Phase 2", plotOutput("classic_phase2", height = "500px"))
      )
    )
  )
)

server <- function(input, output, session) {
  vars <- c("Concentration", "Humidity", "Dissolution", "Density")

  # Simulated data
  sim_batches <- simulate_pharma_batches()

  # Robust Phase 1
  phase1_data <- subset(sim_batches, Phase == "Phase 1" & Status == "Under Control")
  phase1 <- reactive({
    robust_statis_phase1(data = phase1_data, variables = vars)
  })

  # Robust Phase 2
  phase2_data <- subset(sim_batches, Phase == "Phase 2")
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

  # Classical Hotelling - Phase 1
  classical_f1 <- hotelling_t2_phase1(
    data = subset(sim_batches, Phase == "Phase 1"),
    variables = vars
  )

  # Classical Hotelling - Phase 2
  classical_f2 <- hotelling_t2_phase2(
    new_data = phase2_data,
    variables = vars,
    center = classical_f1$center,
    covariance = classical_f1$covariance
  )

  # Phase 2 status (to color classical Phase 2 chart)
  status_f2 <- unique(phase2_data[, c("Batch", "Status")])
  classical_f2$batch_statistics <- merge(classical_f2$batch_statistics, status_f2, by = "Batch")

  # Robust Phase 1 chart
  output$phase1_plot <- renderPlot({
    plot_statis_phase1_chart(batch_statistics = phase1()$batch_statistics, num_vars = length(vars))
  })

  # Robust Phase 2 chart
  output$phase2_plot <- renderPlot({
    plot_statis_phase2_chart(phase2_result = phase2())
  })

  # HJ-Biplot Phase 1
  output$biplot <- renderPlot({
    plot_statis_hj_biplot(phase1_result = phase1(), color_by = input$color_by)
  })

  # Projection Biplot Phase 2
  output$projection_biplot <- renderPlot({
    plot_statis_biplot_projection(
      phase1_result = phase1(),
      phase2_result = phase2()
    )
  })

  # Classical T² - Phase 1
  output$classic_phase1 <- renderPlot({
    plot_classical_hotelling_t2_chart(t2_statistics = classical_f1$batch_statistics, num_vars = length(vars))
  })

  # Classical T² - Phase 2
  output$classic_phase2 <- renderPlot({
    plot_classical_hotelling_t2_phase2_chart(t2_statistics = classical_f2$batch_statistics, num_vars = length(vars))
  })
}

app <- shinyApp(ui = ui, server = server)
