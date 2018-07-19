\donttest{
if (require("INLA", quietly = TRUE)) {

# Load 1D Poisson process data

data(Poisson2_1D, package = "inlabru")

# Take a look at the point (and frequency) data

ggplot(pts2) + 
  geom_histogram(aes(x = x), binwidth = 55/20, boundary = 0, fill = NA, color = "black") +
  geom_point(aes(x), y = 0, pch = "|", cex = 4) + 
  coord_fixed(ratio = 1)

# Fit an LGCP model with  and SPDE component

x <- seq(0, 55, length = 20)
mesh1D <- inla.mesh.1d(x, boundary = "free")
mdl <- x ~ spde1D(map = x, model = inla.spde2.matern(mesh1D)) + Intercept
fit <- lgcp(mdl, pts2, domain = list(x = c(0,55)))

# Calculate and plot the posterior range 

range = spde.posterior(fit, "spde1D", "range")
plot(range)

# Calculate and plot the posterior log range 

lrange = spde.posterior(fit, "spde1D", "log.range")
plot(lrange)

# Calculate and plot the posterior variance

variance = spde.posterior(fit, "spde1D", "variance")
plot(variance)

# Calculate and plot the posterior log variance

lvariance = spde.posterior(fit, "spde1D", "log.variance")
plot(lvariance)

# Calculate and plot the posterior Matern correlation

matcor = spde.posterior(fit, "spde1D", "matern.correlation")
plot(matcor)

# Calculate and plot the posterior Matern covariance

matcov = spde.posterior(fit, "spde1D", "matern.covariance")
plot(matcov)

}
}
