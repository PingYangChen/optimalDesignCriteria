library(globpso)
source("tool.R")
library(Rcpp)
library(RcppArmadillo)
sourceCpp("src/crit.cpp")

# Wrapper function for optimizing optimal Design of Heteroscedastic Polynomial Model
objWrapper <- function(x, n, d, ctype, poly, varpar, a, b) {
  tmp <- deisgnCritCpp(x, n, d, ctype, poly, varpar[1], varpar[2], a, b)
  tmp$crit
}

#' Particle Swarm Optimization Algorithms for Finding E-optimal Design of Heteroscedastic Polynomial Model
#' @param nSupp A integer number of the number of support points. The default is \code{2}.
#' @param designSpace The vector of the design space, \code{c(lower, upper)} and \code{lower < upper}.  The default is \code{c(-1, 1)}.
#' @param polyOrder The vector of the order of polynomial model terms.  The default is \code{c(0, 1)}.
#' \code{0} indicates the intercept term. \code{-Inf} indicates the $\log(x)$ term. Otherwise, set \code{a = numerical value} for $x^a$ term.
#' @param varParam The vector of the variance parameters, \code{c(u, v)} and \code{u >= 0} and \code{v >= 0}.
#' @param designType \code{"D"}, \code{"A"} or \code{"E"}.  Default is \code{"E"}.
#' @param eqvData The logical value controls if the function returns the data for drawing the equivalence plot. The default is \code{TRUE}.
#' @param nswarm A integer number of swarm size in PSO algorithm. The default is \code{32}.
#' @param maxiter A integer number of maximal PSO iterations. The default is \code{100}.
#' @param seed The random seed that controls initial swarm of PSO. The default is \code{NULL}.
#' @param verbose The logical value controls if PSO would reports the updating progress. The default is \code{TRUE}.
#' @references 
#' \enumerate{
#'   \item Dette, H., 1993. A note on E-optimal designs for weighted polynomial regression. The Annals of Statistics 21, 767–771.
#' }
#' @name genHeterPolyE
#' @rdname genHeterPolyE
#' @export
genHeterPolyE <- function(nSupp = 2, designSpace = c(-1, 1), polyOrder = c(0, 1), 
                          varParam = c(1, 1), designType = "E", 
                          eqvData = TRUE, nswarm = 32, maxiter = 100, seed = NULL, verbose = TRUE) {
  # test
  # nSupp=2;designSpace=c(-1,1);polyOrder=c(0,1);varParam=c(1,1);designType="E";eqvData=TRUE;seed=NULL;verbose=TRUE
  
  ### Check input qualites
  stopifnot(
    nSupp >= 2,
    designSpace[1] < designSpace[2],
    sum(is.infinite(polyOrder)) < 2,
    all(varParam >= 0)
  )
  
  ctype <- switch(designType, "D" = 0, "A" = 1, "E" = 2, -1)
  if (ctype == -1) {
    error('designType should be one of c("D", "A", "E")')
  }
  
  infoL <- list(
    "n" = nSupp,
    "d" = 1,
    "dsp" = designSpace,
    "poly" = polyOrder,
    "param" = varParam, # u, v
    "ctype" = ctype
  )

  # Define the search domain 
  low_bound <- c(rep(infoL$dsp[1], infoL$n), rep(pi/2, infoL$n - 2),  0.0)
  upp_bound <- c(rep(infoL$dsp[2], infoL$n),   rep(pi, infoL$n - 2), pi/2)
  
  alg_setting <- getPSOInfo(nSwarm = nswarm, maxIter = maxiter, psoType = "basic")
  
  # Run PSO for finding E-optimal design for heteroscedastic polynomial model
  set.seed(seed)
  res <- globpso(objFunc = objWrapper, lower = low_bound, upper = upp_bound, PSO_INFO = alg_setting, verbose = verbose,
                 n = infoL$n, d = infoL$d, ctype = infoL$ctype, poly = infoL$poly, 
                 varpar = infoL$param, a = infoL$dsp[1], b = infoL$dsp[2])
  
  # Draw the equivalence plot
  if (eqvData) {
    grid <- seq(infoL$dsp[1], infoL$dsp[2], length = 100)
    eqv  <- equivCpp(grid, res$par, infoL$n, infoL$d, ctype = infoL$ctype, 
                     poly = infoL$poly, u = infoL$param[1], v = infoL$param[2], 
                     a = infoL$dsp[1], b = infoL$dsp[2]) 
  }
  
  # Output the results
  return(list(
    "requirements" = infoL,
    "pso_results" = res,
    "optDesign" = designM2V(res$par, infoL$n, infoL$d),
    "CritVal" = res$val,
    "Equiv" = eqv
  ))
}





