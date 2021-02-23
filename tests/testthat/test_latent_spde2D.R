local_bru_testthat_setup()


test_that("Georeferenced data with sp", {
  skip_on_cran()
  local_bru_safe_inla()

  set.seed(123)
  mydata <- expand.grid(
    Easting = seq(5, 45, by = 20),
    Northing = seq(10, 30, by = 10),
    KEEP.OUT.ATTRS = FALSE
  )
  mydata[["obs"]] <- (mydata$Easting - 20) / 10 + rnorm(NROW(mydata))
  coordinates(mydata) <- c("Easting", "Northing")

  mesh <- INLA::inla.mesh.2d(
    loc = mydata,
    offset = 5,
    max.edge = 4,
    n = 16
  )

  matern <- INLA::inla.spde2.pcmatern(
    mesh,
    prior.sigma = c(10, 0.01),
    prior.range = c(4, 0.01)
  )

  # Check that mistaken empty or unnamed arguments are detected
  cmp <- obs ~ Intercept(1) + field(coordinates, model = matern, )
  expect_error(
    component_list(cmp),
    "Unnamed arguments detected in component .* position\\(s\\) 3"
  )

  cmp <- obs ~ Intercept(1) + field(coordinates, model = matern)

  fit <- bru(
    cmp,
    data = mydata,
    options = list(
      control.inla = list(
        int.strategy = "eb"
      )
    )
  )

  # Check Intercept
  expect_equal(
    fit$summary.fixed["Intercept", "mean"],
    0.5398535,
    tolerance = midtol
  )

  # Check SPDE
  expect_equal(
    fit$summary.random$field$mean[mesh$idx$loc[1:3]],
    c(-2.6003077, -0.2699909, 3.5188725),
    tolerance = midtol
  )


  # Check that explicit access to the data object works
  cmp <- obs ~ Intercept(1) + field(coordinates(.data.), model = matern)

  fit <- bru(
    cmp,
    data = mydata,
    options = list(
      control.inla = list(
        int.strategy = "eb"
      )
    )
  )

  # Check Intercept
  expect_equal(
    fit$summary.fixed["Intercept", "mean"],
    0.5398535,
    tolerance = midtol
  )

  # Check SPDE
  expect_equal(
    fit$summary.random$field$mean[mesh$idx$loc[1:3]],
    c(-2.6003077, -0.2699909, 3.5188725),
    tolerance = midtol
  )


  # Check that explicit access to the data object works
  cmp <- obs ~ Intercept(1) + field(cbind(
    as.data.frame(.data.)$Easting,
    as.data.frame(.data.)$Northing
  ),
  model = matern
  )

  fit <- bru(
    cmp,
    data = mydata,
    options = list(
      control.inla = list(
        int.strategy = "eb"
      )
    )
  )

  # Check Intercept
  expect_equal(
    fit$summary.fixed["Intercept", "mean"],
    0.5398535,
    tolerance = midtol
  )

  # Check SPDE
  expect_equal(
    fit$summary.random$field$mean[mesh$idx$loc[1:3]],
    c(-2.6003077, -0.2699909, 3.5188725),
    tolerance = midtol
  )

  pred_df <- pixels(mesh)
  expect_s4_class(pred_df, "SpatialPixelsDataFrame")
  pred <- predict(fit, pred_df, ~ exp(Intercept + field), n.samples = 5)
  expect_s4_class(pred, "SpatialPixelsDataFrame")
})




latent_spde2D_group_testdata <- function() {
  set.seed(123)

  # Load and reduce data set
  data(mrsea, package = "inlabru")
  mrsea <- local_mrsea_convert(mrsea, use_km = TRUE)
  coordnames(mrsea$points) <- c("Easting", "Northing")
  coordnames(mrsea$samplers) <- c("Easting", "Northing")

  mrsea$points <- mrsea$points[mrsea$points$season %in% c(1, 2), ]
  mrsea$samplers <- mrsea$samplers[mrsea$samplers$season %in% c(1, 2), ]

  # Integration points
  ips <- ipoints(mrsea$samplers, domain = mrsea$mesh, group = "season")

  # Run the model
  matern <- INLA::inla.spde2.pcmatern(mrsea$mesh,
    prior.sigma = c(0.1, 0.01),
    prior.range = c(10, 0.01)
  )

  cmp <-
    coordinates + season ~
    mySmooth(
      main = coordinates, model = matern,
      group = season, ngroup = 2
    ) + Intercept(1)
  fit <- lgcp(cmp, mrsea$points,
    ips = ips,
    options = list(
      control.inla = list(
        int.strategy = "eb"
      )
    )
  )

  list(
    mrsea = mrsea,
    matern = matern,
    cmp = cmp,
    fit = fit
  )
}

test_that("Latent models: SPDE with group parameter (spatiotemporal)", {
  skip_on_cran()
  local_bru_safe_inla()
  expect_warning(
    {
      data_ <- latent_spde2D_group_testdata()
    },
    "export to PROJ failed: generic error of unknown origin"
  )

  # Check Intercept
  expect_equal(
    data_$fit$summary.fixed["Intercept", "mean"],
    -2.206082,
    tolerance = midtol
  )

  # Check SPDE
  expect_equal(
    data_$fit$summary.random$mySmooth$mean[c(1, 250, 550)],
    c(-0.1618776, 0.7721959, 2.0314753),
    tolerance = midtol
  )
  expect_equal(
    data_$fit$summary.random$mySmooth$sd[c(1, 250, 550)],
    c(1.9784044, 0.7738195, 0.5835616),
    tolerance = midtol
  )
  expect_error(spde.posterior(data_$fit, "mySmooth", what = "range"), NA)
})
