\donttest{
if (bru_safe_inla() &&
    require("sn", quietly = TRUE) &&
    require("ggplot2", quietly = TRUE) &&
    require("terra", quietly = TRUE) &&
    require("sf", quietly = TRUE)) {

  # Load the Gorilla data

  gorillas <- gorillas_sf

  # Plot the Gorilla nests, the mesh and the survey boundary

  ggplot() +
    gg(gorillas$mesh) +
    gg(gorillas$nests) +
    gg(gorillas$boundary, alpha = 0.1)

  # Define SPDE prior

  matern <- INLA::inla.spde2.pcmatern(
    gorillas$mesh,
    prior.sigma = c(0.1, 0.01),
    prior.range = c(0.01, 0.01)
  )

  # Define domain of the LGCP as well as the model components (spatial SPDE effect and Intercept)

  cmp <- geometry ~ field(geometry, model = matern) + Intercept(1)

  # Fit the model, with "eb" instead of full Bayes
  fit <- lgcp(
    cmp,
    data = gorillas$nests,
    samplers = gorillas$boundary,
    domain = list(geometry = gorillas$mesh),
    options = list(control.inla = list(int.strategy = "eb"))
  )

  # Once we obtain a fitted model the predict function can serve various purposes.
  # The most basic one is to determine posterior statistics of a univariate
  # random variable in the model, e.g. the intercept

  icpt <- predict(fit, NULL, ~ c(Intercept = Intercept_latent))
  plot(icpt)

  # The formula argument can take any expression that is valid within the model, for
  # instance a non-linear transformation of a random variable

  exp.icpt <- predict(fit, NULL, ~ c(
    "Intercept" = Intercept_latent,
    "exp(Intercept)" = exp(Intercept_latent)
  ))
  plot(exp.icpt, bar = TRUE)

  # The intercept is special in the sense that it does not depend on other variables
  # or covariates. However, this is not true for the smooth spatial effects 'field'.
  # In order to predict 'field' we have to define where (in space) to predict. For
  # this purpose, the second argument of the predict function can take \code{data.frame}
  # objects as well as Spatial objects. For instance, we might want to predict
  # 'field' at the locations of the mesh vertices. Using

  vrt <- fm_vertices(gorillas$mesh, format = "sf")

  # we obtain these vertices as a SpatialPointsDataFrame

  ggplot() +
    gg(gorillas$mesh) +
    gg(vrt, color = "red")

  # Predicting 'mySmooth' at these locations works as follows

  field <- predict(fit, vrt, ~field)

  # Note that just like the input also the output will be a sf object with
  # points and that the predicted statistics are simply added as columns

  class(field)
  head(vrt)
  head(field)

  # Plotting the mean, for instance, at the mesh node is straight forward

  ggplot() +
    gg(gorillas$mesh) +
    gg(field, aes(color = mean), size = 2)

  # However, we are often interested in a spatial field and thus a linear interpolation,
  # which can be achieved by using the gg mechanism for meshes

  ggplot() +
    gg(gorillas$mesh, color = field$mean)

  # Alternatively, we can predict the spatial field at a grid of locations, e.g. a
  # sf object with a grid of points covering the relevant part of mesh

  pxl <- fm_pixels(gorillas$mesh, format = "sf", mask = gorillas$boundary)
  field2 <- predict(fit, pxl, ~field)

  # This will give us a sf with the columns we are looking for

  head(field2)
  ggplot() +
    gg(gorillas$boundary) +
    gg(data = field2, geom = "tile")
}
}
