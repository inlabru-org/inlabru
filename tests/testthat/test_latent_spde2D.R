local_bru_testthat_setup()

latent_spde2D_group_testdata <- function(num.threads = NULL,
                                         tolerance = NULL,
                                         h = 0.005) {
  set.seed(123)

  # Load and reduce data set
  data(mrsea, package = "inlabru")
  mrsea <- local_mrsea_rebuild_CRS(mrsea, use_km = TRUE)
  mrsea$points <- mrsea$points[mrsea$points$season == 1 |
    mrsea$points$season == 2, ]
  mrsea$samplers <- mrsea$samplers[mrsea$samplers$season == 1 |
    mrsea$samplers$season == 2, ]

  # Integration points
  ips <- ipoints(mrsea$samplers, domain = mrsea$mesh, group = "season")

  # Run the model
  matern <- INLA::inla.spde2.pcmatern(mrsea$mesh,
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
    options = list(
      control.inla = list(
        int.strategy = "eb",
        tolerance = tolerance,
        h = h
      ),
      num.threads = num.threads
    )
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
  local_bru_safe_inla()
  expect_warning(
    {
      data <- latent_spde2D_group_testdata(h = 0.005)
    },
    "export to PROJ failed: generic error of unknown origin"
  )

  # Check Intercept
  expect_equal(
    data$fit$summary.fixed["Intercept", "mean"],
    -8.8628,
    tolerance = midtol
  )

  # Check SPDE
  expect_equal(
    data$fit$summary.random$mySmooth$mean[c(1, 250, 550)],
    c(-0.8247674, -2.3758650, 0.9492320),
    tolerance = midtol
  )
  expect_equal(
    data$fit$summary.random$mySmooth$sd[c(1, 250, 550)],
    c(0.9502163, 1.1058397, 0.6451496),
    tolerance = midtol
  )
  expect_error(spde.posterior(data$fit, "mySmooth", what = "range"), NA)
})
