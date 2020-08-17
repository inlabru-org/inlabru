context("Latent models - 2D SPDE - Group parameter (test_latent_spde2D.R)")

latent_spde2D_group_testdata <- function(num.threads = 1,
                                         tolerance = NULL,
                                         h = 0.005) {
  disable_PROJ6_warnings()
  set.seed(123)

  # Load and reduce data set
  data(mrsea)
  mrsea <- mrsea_rebuild_CRS(mrsea)
  mrsea$points <- mrsea$points[mrsea$points$season == 1 |
    mrsea$points$season == 2, ]
  mrsea$samplers <- mrsea$samplers[mrsea$samplers$season == 1 |
    mrsea$samplers$season == 2, ]

  # The estimation is numerically unreliable when the spatial
  # domain is represented in metres, and has been seen to produce
  # different results on different systems (e.g. Travis CI).
  # Transform m to km:
    crs_km <- inla.crs_set_lengthunit(mrsea$mesh$crs, "km")
  mrsea$mesh <- inla.spTransform(mrsea$mesh, crs_km)
  mrsea$samplers <- sp::spTransform(mrsea$samplers, crs_km)

  # Integration points
  ips <- ipoints(mrsea$samplers, group = "season")

  # Run the model
  matern <- inla.spde2.pcmatern(mrsea$mesh,
    prior.sigma = c(0.1, 0.01),
    prior.range = c(10, 0.01)
  )

  cmp <- coordinates + season ~
  mySmooth(
    main = coordinates, model = matern,
    group = season, ngroup = 2
  ) + Intercept(1)
  fit <- lgcp(cmp, mrsea$points,
    ips = ips,
    options = list(control.inla = list(int.strategy = "eb",
                                       tolerance = tolerance,
                                       h = h),
                   num.threads = num.threads)
  )

  data <-
  list(
    mrsea = mrsea,
    matern = matern,
    cmp = cmp,
    fit = fit
  )
  data
}

test_that("Latent models: SPDE with group parameter (spatiotemporal)", {
  skip_on_cran()
  disable_PROJ6_warnings()
  library(INLA)
  data <- latent_spde2D_group_testdata(num.threads = 1, h = 0.005)

  # Check Intercept
  expect_equal(data$fit$summary.fixed["Intercept", "mean"], -9.231591, midtol)

  # Check SPDE
  expect_equal(
    data$fit$summary.random$mySmooth$mean[c(1, 250, 550)],
    c(-1.27410, -1.96819, 0.7479571), midtol
  )
  expect_equal(
    data$fit$summary.random$mySmooth$sd[c(1, 250, 550)],
    c(1.439, 1.635, 0.935), midtol
  )
})
