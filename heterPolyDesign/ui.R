library(shiny)

ui <- fluidPage(
  
  # App title
  titlePanel("My Shiny App"),
  
  # Sidebar + main layout
  sidebarLayout(
    sidebarPanel(
      numericInput("poly", "poly", value = 1, min = 1, max = 6),
      numericInput("u", "u", value = 1, min = 0, max = 4),
      numericInput("v", "v", value = 1, min = 0, max = 4),
      numericInput("upp", "Upper Bound", value = -1, min = -10, max = 10),
      numericInput("low", "Lower Bound", value = 1, min = 10, max = 10),
      numericInput("num", "Number of Support Points", value = 2, min = 2, max = 20),
      selectInput("ctype", "Design Criterion", choices = c("D", "E", "A"), selected = "E"),
      actionButton("go", "Run")
    ),
    
    mainPanel(
      h3("Model"),
      textOutput("greeting"),
      h3("Results"),
      textOutput("greeting"),
      verbatimTextOutput("summary")
    )
  )
)