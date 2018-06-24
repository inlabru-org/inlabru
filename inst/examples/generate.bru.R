\donttest{
if (require("INLA", quietly = TRUE)) {
  
# Generate data for a simple linear model

input.df <- data.frame(x=cos(1:10))
input.df <- within(input.df, y <- 5 + 2*cos(1:10) + rnorm(10, mean=0, sd=0.1))

# Fit the model

fit <- bru(y ~ xeff(map = x, model = "linear"), "gaussian", input.df)
summary(fit)

# Generate samples for some predefined x

df = data.frame(x = seq(-4, 4, by = 0.1))
smp = generate(fit, df, ~ xeff + Intercept, n.samples = 10)

# Plot the resulting realizations

plot(df$x, smp[[1]], type = "l")
for (k in 2:length(smp)) points(df$x, smp[[k]], type = "l")

# We can also draw samples form the joint posterior

df = data.frame(x = 1)
smp = generate(fit, df, ~ data.frame(xeff, Intercept), n.samples = 10)
smp[[1]]

# ... and plot them

plot(do.call(rbind, smp))

}
}
