library(shiny)
library(bslib)

ui <- fluidPage(
  withMathJax(),
  # App title
  titlePanel("E-optimal Design"),
  
  # Sidebar + main layout
  sidebarLayout(
    sidebarPanel(
      h3("Model Settings"),
      numericInput("poly", "poly", value = 1, min = 1, max = 6),
      input_switch("log", "Include Log Term", value = 0),
      numericInput("u", "u", value = 1, min = 0, max = 4),
      numericInput("v", "v", value = 1, min = 0, max = 4),
      h3("Design Requirements"),
      conditionalPanel("input.log == 0",
                       numericInput("upp", "Upper Bound", value = 1, min = -10, max = 10),
                       numericInput("low", "Lower Bound", value = -1, min = -10, max = 10),
      ),
      conditionalPanel("input.log == 1",
                       numericInput("upp", "Upper Bound", value = 80, min = 1, max = 200),
                       numericInput("low", "Lower Bound", value = 20, min = 1, max = 200),
      ),
      numericInput("num", "Number of Support Points", value = 2, min = 2, max = 20),
      selectInput("ctype", "Design Criterion", choices = c("D", "E", "A"), selected = "E"),
      h3("PSO Settings"),
      numericInput("nswarm",  "Swarm Size", value = 32, min = 4, max = 256),
      numericInput("maxiter", "Number of Iterations", value = 100, min = 20, max = 1000),
      actionButton("go", "Run")
    ),
    
    mainPanel(
      h3("Model"),
      uiOutput("modelform"),
      h3("Design"),
      tableOutput("design"),
      h3("Equivalence Plot"),
      imageOutput("eqvplot")
    )
  )
)