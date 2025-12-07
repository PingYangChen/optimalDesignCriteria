library(shiny)

library(Rcpp)
library(RcppArmadillo)
library(globpso)

sourceCpp("../src/crit.cpp")
source("../tool.R")

server <- function(input, output, session) {
  
  # Reactive expression
  data_reactive <- reactive({
    list(
      name = input$name,
      squared = input$num^2
    )
  })
  
  # Output for greeting
  output$greeting <- renderText({
    paste("Hello,", data_reactive()$name, "!")
  })
  
  # Output for number summary
  output$summary <- renderPrint({
    list(
      original_number = input$num,
      squared_number  = data_reactive()$squared
    )
  })
}