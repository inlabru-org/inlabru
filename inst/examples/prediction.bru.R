\donttest{
# Generate some data
input.df <- data.frame(x=cos(1:10))
input.df <- within(input.df, y <- 5 + 2*cos(1:10) + rnorm(10, mean=0, sd=0.1))

# Fit a model with fixed effect 'x' and intercept 'Intercept'

fit <- bru(y ~ x, "gaussian", input.df)

# Predict posterior statistics of 'x'

xpost = predict(fit, formula = ~ x)

# The result is a data.frame inheriting from class 'prediction'

class(xpost)

# The statistics include mean, standard deviation, the 2.5% quantile, the median,
# the 97.5% quantile, minimum and maximum sample drawn from the posterior as well as
# the coefficient of variation and the variance.

xpost

# For a single variable like 'x' the default plotting method invoked by gg() will
# show these statisics in a fashion similar to a box plot:

ggplot() + gg(xpost)


# The predict function can also be used to simulatenneously estimate posteriors
# of multiple variables:

xipost = predict(fit, formula = ~ data.frame(post = c(x, Intercept)))
xipost
                 
# If we still want a plot in the previous style we have to set the 'bar' parameter to TRUE:

rownames(xipost) = c("x","Intercept")
ggplot() + gg(xipost, bar = TRUE)
                 
}