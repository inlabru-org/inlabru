\donttest{
if (require("INLA", quietly = TRUE)) {

# Simulate some covariates x and observations y
input.df <- data.frame(x=cos(1:10))
input.df <- within(input.df, y <- 5 + 2*x + rnorm(10, mean=0, sd=0.1))

# Fit a Gaussian likelihood model
fit <- bru(y ~ x + Intercept, "gaussian", input.df)

# Obtain summary
fit$summary.fixed
}

  
if (require("INLA", quietly = TRUE)) {
  
# Alternatively, we can use the like() function to construct the likelihood:

lik = like(family = "gaussian", data = input.df)
fit <- bru(y ~ x + Intercept, lik)
fit$summary.fixed

}
  
# An important addition to the INLA methodology is bru's ability to use
# non-linear predictors. Such a predictor can be formulated via like()'s 
# \code{formula} parameter. For instance

if (require("INLA", quietly = TRUE)) {
    
z = 2
input.df <- within(input.df, y <- 5 + exp(z)*x + rnorm(10, mean=0, sd=0.1))
lik = like(family = "gaussian", data = input.df, formula = y ~ exp(z)*x + Intercept, E = 10000)
fit <- bru( ~ z + Intercept, lik)

# Check the result (z posterior should be around 2)
fit$summary.fixed
}
  
}
