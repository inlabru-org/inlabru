\donttest{
if (bru_safe_inla()) {
  # Simulate some covariates x and observations y
  input.df <- data.frame(x = cos(1:10))
  input.df <- within(input.df, {
    y <- 5 + 2 * x + rnorm(10, mean = 0, sd = 0.1)
  })

  # Fit a Gaussian likelihood model
  fit <- bru(y ~ x + Intercept(1), family = "gaussian", data = input.df)

  # Obtain summary
  fit$summary.fixed
}


if (bru_safe_inla()) {
  # Alternatively, we can use the bru_obs() function to construct the likelihood:

  lik <- bru_obs(family = "gaussian",
              formula = y ~ x + Intercept,
              data = input.df)
  fit <- bru(~ x + Intercept(1), lik)
  fit$summary.fixed
}

# An important addition to the INLA methodology is bru's ability to use
# non-linear predictors. Such a predictor can be formulated via bru_obs()'s
# \code{formula} parameter. The z(1) notation is needed to ensure that
# the z component should be interpreted as single latent variable and not
# a covariate:

if (bru_safe_inla()) {
  z <- 2
  input.df <- within(input.df, {
    y <- 5 + exp(z) * x + rnorm(10, mean = 0, sd = 0.1)
  })
  lik <- bru_obs(
    family = "gaussian", data = input.df,
    formula = y ~ exp(z) * x + Intercept
  )
  fit <- bru(~ z(1) + Intercept(1), lik)

  # Check the result (z posterior should be around 2)
  fit$summary.fixed
}
}
