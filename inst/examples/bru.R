# Simulate some covariates x and observations y
input.df <- data.frame(x=cos(1:10))
input.df <- within(input.df, y <- 5 + 2*cos(1:10) + rnorm(10, mean=0, sd=0.1))

# Fit a Gaussian likelihood model
fit <- bru(y ~ x, "gaussian", input.df)

# Obtain summary
summary(fit)
