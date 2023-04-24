local_bru_testthat_setup()

test_that("1D integration points can be generated", {
  local_bru_safe_inla()
  sf_obj1 <- sf::st_as_sf(data.frame(x=1:3, y=3:5), coords=c("x","y"))
  sf_obj2 <- sf::st_as_sf(data.frame(x=3:6, y=5:8), coords=c("x","y"))

  samplers <- list(data.frame(x=1:4), sf_obj1, sf_obj2)
  ips <- ipoints(c(0, 10), 3, name = "myDim")

  expect_s3_class(ips, "data.frame")
  expect_equal(nrow(ips), 3)
  expect_equal(ncol(ips), 3)
  expect_equal(names(ips), c("myDim", "weight", "group"))
  expect_equal(as.numeric(ips[1, ]), c(5 / 3, 10 / 3, 1))
  expect_equal(as.numeric(ips[2, ]), c(15 / 3, 10 / 3, 1))
  expect_equal(as.numeric(ips[3, ]), c(25 / 3, 10 / 3, 1))

  # mesh1D
  data(Poisson2_1D, package = "inlabru", envir = environment())
  x <- seq(0, 55, length = 50)
  mesh1D <- INLA::inla.mesh.1d(x, boundary = "free")
  matern <- INLA::inla.spde2.pcmatern(
    mesh1D,
    prior.range = c(1, 0.01),
    prior.sigma = c(0.1, 0.75)
  )

  mdl <- ~ spde1D(main = x, model = matern) + Intercept(1)

  ips <- apmaker(domain = list(x = mesh1D))

  fit1 <- bru(
    mdl,
    like(
      formula = x ~ .,
      family = "cp",
      data = pts2,
      domain = list(x = mesh1D),
      options = list(bru_compress_cp = FALSE)
    )
  )

  # 2d_lgcp_st
  library(inlabru)
  library(INLA)
  library(ggplot2)
  data(mrsea, package = "inlabru")
  matern <- INLA::inla.spde2.pcmatern(mrsea$mesh,
                                prior.sigma = c(0.1, 0.01),
                                prior.range = c(10, 0.01)
  )

  cmp <- coordinates + season ~ Intercept(1) +
    mySmooth(
      coordinates,
      model = matern,
      group = season,
      ngroup = 4
    )

  fit <- lgcp(cmp,
              data = mrsea$points,
              samplers = mrsea$samplers,
              domain = list(
                coordinates = mrsea$mesh,
                season = seq_len(4)
              )
  )

# gorillas
    data(gorillas, package = "inlabru")

    # Define SPDE prior
    matern <- INLA::inla.spde2.pcmatern(gorillas$mesh,
      prior.sigma = c(0.1, 0.01),
      prior.range = c(0.01, 0.01)
    )

    # Define domain of the LGCP as well as the model components (spatial SPDE
    # effect and Intercept)
    cmp <- coordinates ~ mySmooth(coordinates, model = matern) + Intercept(1)

    # Fit the model (with int.strategy="eb" to make the example take less time)
    fit <- lgcp(cmp, gorillas$nests,
      samplers = gorillas$boundary,
      domain = list(coordinates = gorillas$mesh),
      options = list(control.inla = list(int.strategy = "eb"))
    )

})

