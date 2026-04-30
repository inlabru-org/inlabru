test_that("Joint gcpo", {
  skip_on_cran()
  local_bru_safe_inla()
  skip_if_not_installed("sf")

  ## partition
  cvpart <- cv_hex(
    gorillas_sf$boundary,
    cellsize = 0.5,
    n_group = 3,
    resolution = c(95, 80)
  )
  cvpart$block_ID <- seq_len(nrow(cvpart))
  cvpart$group <- NULL
  nblock <- nrow(cvpart)

  ## split
  gorillas_nests_major <-
    gorillas_sf$nests[gorillas_sf$nests$group == "major", ]
  gorillas_nests_minor <-
    gorillas_sf$nests[gorillas_sf$nests$group == "minor", ]

  # blocks for each nest
  a <- sf::st_intersects(gorillas_nests_major, cvpart)
  if (!all(vapply(a, function(x) (length(x) == 1), logical(1)))) {
    stop("Point in none or multiple polygons")
  }
  gorillas_nests_major$.block <- unlist(a)

  b <- sf::st_intersects(gorillas_nests_minor, cvpart)
  if (!all(vapply(b, function(x) (length(x) == 1), logical(1)))) {
    stop("Point in none or multiple polygons")
  }
  gorillas_nests_minor$.block <- unlist(b)

  # model
  matern <- INLA::inla.spde2.pcmatern(
    gorillas_sf$mesh,
    prior.range = c(0.1, 0.01),
    prior.sigma = c(1, 0.01)
  )

  cmp <- ~
    Common(geometry, model = matern) +
      Difference(geometry, model = matern) +
      Intercept(1)

  fml.major <- geometry ~ Intercept + Common + Difference / 2
  fml.minor <- geometry ~ Intercept + Common - Difference / 2

  lik_major <- bru_obs("cp",
    formula = fml.major,
    samplers = cvpart,
    data = gorillas_nests_major,
    domain = list(geometry = gorillas_sf$mesh),
    control.gcpo = list(
      enable = TRUE,
      type.cv = "joint",
      num.level.sets = -1
    ),
    tag = "major"
  )
  lik_minor <- bru_obs("cp",
    formula = fml.minor,
    samplers = cvpart,
    data = gorillas_nests_minor,
    domain = list(geometry = gorillas_sf$mesh),
    control.gcpo = list(
      enable = TRUE,
      type.cv = "joint",
      num.level.sets = -1
    ),
    tag = "minor"
  )

  c_g <- bru_obs_control_gcpo(
    c(lik_major, lik_minor),
    control.gcpo = bru_options("control.gcpo")
  )
  expect_length(
    c_g$friends,
    length(bru_index(lik_major)) + length(bru_index(lik_minor))
  )
})
