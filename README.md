# Using Particle Swarm Optimization Algorithms for Finding E-optimal Design of Heteroscedastic Polynomial Model

Please find the example in the `exmaple.R` file.

```{r}
# Import the main function
source("Approx_HeterPoly.R")

### Case 1: 
### The Polynomial Model is f(x) = 1 + x

# Set the design requirements
n <- 2                # the number of support points
dSpace <- c(-1, 1)    # the design space, c(lower, upper)
pOrder <- c(0, 1)     # the order of polynomial model terms
                      # 0: intercept
                      # -Inf: log(x)
                      # a: x^a, a is numerical value
vParam <- c(1.5, 0.5) # the variance parameters, c(u, v)


# Run PSO for Finding E-optimal Design of Heteroscedastic Polynomial Model
out <- genHeterPolyE(nSupp = n, designSpace = dSpace, polyOrder = pOrder, 
              varParam = vParam, designType = "E", 
              eqvData = TRUE, nswarm = 32, maxiter = 100, seed = NULL, verbose = TRUE)

# The PSO-generated Design
out$optDesign

# The corresponding E-criterion value 
# (NEGATIVE of the minimal eignevalue of FIM)
out$CritVal

# Draw the equivalence plot
plot(out$Equiv$grid, out$Equiv$DD, type = "l", col = "blue")
abline(h = 0, lty = 3)
```


```{r}
### Case 2: 
### The Polynomial Model is f(x) = 1 + x^(-2) + x^(-1) + x^(-0.5) + log(x)

# Set the design requirements
n <- 5                # the number of support points
dSpace <- c(20, 80)    # the design space, c(lower, upper)
pOrder <- c(0, -2, -1, -0.5, -Inf)     # the order of polynomial model terms
# 0: intercept
# Inf: log(x)
# a: x^a, a is numerical value
vParam <- c(1.5, 0.5) # the variance parameters, c(u, v)


# Run PSO for Finding E-optimal Design of Heteroscedastic Polynomial Model
out <- genHeterPolyE(nSupp = n, designSpace = dSpace, polyOrder = pOrder, 
                     varParam = vParam, designType = "E", 
                     eqvData = TRUE, nswarm = 32, maxiter = 100, seed = NULL, verbose = TRUE)

# The PSO-generated Design
out$optDesign

# The corresponding E-criterion value 
# (NEGATIVE of the minimal eignevalue of FIM)
out$CritVal

# Draw the equivalence plot
plot(out$Equiv$grid, out$Equiv$DD, type = "l", col = "blue")
abline(h = 0, lty = 3)
```