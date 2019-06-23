basic_intercept_testdata <- function() {
  set.seed(123)
  data.frame(
    Intercept = 1,
    y = rnorm(100)
  )
}

basic_fixed_effect_testdata <- function() {
  cbind(
    basic_intercept_testdata(),
    data.frame(x1 = rnorm(100))
  )
}



gorillas_lgcp_2d_testdata <- function() {
  data(gorillas, package = "inlabru")
  matern <- inla.spde2.pcmatern(gorillas$mesh, prior.sigma = c(0.1, 0.01), prior.range = c(5, 0.01))
  cmp <- coordinates ~ mySmooth(map = coordinates, model = matern) + Intercept

  fit <- lgcp(cmp, gorillas$nests,
    samplers = gorillas$boundary,
    options = list(control.inla = list(int.strategy = "eb"))
  )

  list(
    gorillas = gorillas,
    matern = matern,
    cmp = cmp,
    fit = fit
  )
}
