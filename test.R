library(Rcpp)
library(RcppArmadillo)
library(globpso)

sourceCpp("src/crit.cpp")
source("tool.R")

critFunc <- function(x, n, d, poly, varpar, a, b) {
  tmp <- deisgnCritCpp(x, n, d, poly, varpar[1], varpar[2], a, b)
  tmp$crit
}

nSupp <- 2
vp <- c(0.5, 1.5)
# The search domain is [-5, 5]^3
upp_bound <- c(rep(1, nSupp), rep(pi, nSupp - 2), pi/2)
low_bound <- c(rep(-1, nSupp), rep(pi/2, nSupp - 2),  0.0)

# Run PSO for this optimization problem
# Also input the enviorment variable, the location shift 'loc_shift'
res <- globpso(objFunc = critFunc, lower = low_bound, upper = upp_bound, 
               n = nSupp, d = 1, poly = c(0, 1), varpar = vp, a = -1, b = 1)
res$par
res$val

designM2V(res$par, 2, 1)


grid <- seq(-1, 1, length = 100)
eqv  <- equivCpp(grid, res$par, nSupp, 1, 2, poly = c(0, 1), u = vp[1], v = vp[2], a = -1, b = 1) 

plot(eqv$grid, eqv$DD, type = "l", col = "blue")
abline(h = 0, lty = 3)
  



