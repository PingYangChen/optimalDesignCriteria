library(shiny)
library(globpso)
source("../tool.R")
library(Rcpp)
library(RcppArmadillo)
sourceCpp("../src/crit.cpp")
source("../Approx_HeterPoly.R")

#shiny::runApp("heterPolyDesign")

server <- function(input, output, session) {
  
  appData <- reactiveValues(pOrder = NULL, vParam = NULL, dSpace = NULL, n = NULL, ctype = NULL, nswarm = NULL, maxiter = NULL, out = NULL)
  
  observeEvent(input$go, { 
    appData$pOrder <- 0:input$poly
    if (input$log) {
      appData$pOrder <- c(appData$pOrder, -Inf)
    }
    appData$dSpace <- c(input$low, input$upp)
    appData$vParam <- c(input$u, input$v)
    appData$n <- input$num
    appData$ctype <- input$ctype
    appData$nswarm <- input$nswarm
    appData$maxiter <- input$maxiter
    
    appData$out <- genHeterPolyE(nSupp = appData$n, designSpace = appData$dSpace, polyOrder = appData$pOrder, 
                         varParam = appData$vParam, designType = input$ctype, 
                         eqvData = TRUE, nswarm = appData$nswarm, maxiter = appData$maxiter, seed = NULL, verbose = FALSE)
    
  })

  output$eqvplot <- renderImage({
    
    fw <- 480
    fh <- 400
    outfile <- tempfile(fileext='.png')
    png(outfile, width = fw, height = fh)
    
    if (!is.null(appData$out)) {
      plot(appData$out$Equiv$grid, appData$out$Equiv$DD, type = "l", col = "blue", xlab = "x", ylab = "d.d.")
      abline(h = 0, lty = 3)
      #title("Equivalence Plot", line = 1.6)
    }
    
    dev.off()
    list(src = outfile, contentType = 'image/png', width = fw, height = fh, alt = "")
  }, deleteFile = TRUE)
  
  output$design <- renderTable({
    
    if (!is.null(appData$out)) {
      cbind("X" = rownames(appData$out$optDesign), round(appData$out$optDesign, 4))
    }
    
  })
  
  output$modelform <- renderUI({
    polyOrder <- 0:input$poly
    mterms <- sprintf('\\theta_{%d} x^{%d}', polyOrder, polyOrder)
    mterms[which(mterms == '\\theta_{0} x^{0}')] <- '\\theta_{0}'
    mterms[which(mterms == '\\theta_{1} x^{1}')] <- '\\theta_{1} x'
    if (input$log) {
      mterms <- c(mterms, '\\log{(x)}')
    }
    mstr <- paste0(mterms, collapse = " + ")
    withMathJax(sprintf("$$f(x) = %s \\\\ \\text{ and } \\lambda(x) = [x - (%.1f)]^{%.1f}[(%.1f) - x]^{%.1f}$$", mstr, input$low, input$u, input$upp, input$v))
  })
  

}