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

  # dolphin
  envir <- new.env()
  data("mexdolphins", package = "dsm", envir = envir)
  mexdolphins <- as.list(envir)
  data.p4s <- "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
  data_crs <- sp::CRS(data.p4s)

  dset <- import.dsmdata(mexdolphins, covar.col = 8)
  mexdolphin <- as.spatial.dsdata(dset, cnames = c("x", "y"), crs = data_crs)

  # Target CRS
  target.p4s <- "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs +towgs84=0,0,0"
  target_crs <- sp::CRS(target.p4s)

  # Units to km
  mexdolphin$points <- fm_transform(mexdolphin$points, crs = target_crs)
  mexdolphin$samplers <- fm_transform(mexdolphin$samplers, crs = target_crs)
  mexdolphin$points$distance <- mexdolphin$points$distance / 1000
  mexdolphin$points$Effort <- mexdolphin$points$Effort / 1000
  mexdolphin$points$mid.x <- mexdolphin$points$mid.x / 1000
  mexdolphin$points$mid.y <- mexdolphin$points$mid.y / 1000
  mexdolphin$points$start.x <- mexdolphin$points$start.x / 1000
  mexdolphin$points$start.y <- mexdolphin$points$start.y / 1000
  mexdolphin$points$end.x <- mexdolphin$points$end.x / 1000
  mexdolphin$points$end.y <- mexdolphin$points$end.y / 1000
  mexdolphin$samplers$mid.x <- mexdolphin$samplers$mid.x / 1000
  mexdolphin$samplers$mid.y <- mexdolphin$samplers$mid.y / 1000
  mexdolphin$samplers$Effort <- mexdolphin$samplers$Effort / 1000
  mexdolphin$mesh$loc <- mexdolphin$mesh$loc / 1000

  # Set mesh crs
  mexdolphin$mesh$crs <- target_crs

  coordnames(mexdolphin$samplers) <- coordnames(mexdolphin$points)

  # Remove all-NA columns such as "distance" from samplers, since they may
  # interfere with normal usage.
  all_na <- vapply(
    names(mexdolphin$samplers),
    function(nm) {
      all(is.na(mexdolphin$samplers[[nm]]))
    },
    TRUE
  )
  mexdolphin$samplers <- mexdolphin$samplers[!all_na]

  ##### Prediction polygon
  polyloc <- as.data.frame(
    mexdolphin$mesh$loc[mexdolphin$mesh$segm$int$idx[, 1], c(1, 2)]
  )
  colnames(polyloc) <- c("x", "y")
  po <- sp::Polygon(polyloc, hole = FALSE)
  pos <- sp::Polygons(list(po), ID = "c")
  predpoly <- sp::SpatialPolygons(list(pos), proj4string = target_crs)
  df <- data.frame(weight = 1)
  rownames(df) <- "c"
  predpolyd <- sp::SpatialPolygonsDataFrame(predpoly, data = df)
  # plot(predpolyd)
  mexdolphin$ppoly <- predpolyd

##### Simulate a whole population #####
  distance <- seq(0, 8, length.out = 20)
  matern <- INLA::inla.spde2.pcmatern(
    mexdolphin$mesh,
    prior.range = c(50, 0.01),
    prior.sigma = c(2, 0.01)
  )
  cmps <-
    coordinates + distance ~
    Intercept(1) +
    df.lsigma(1) +
    spat(coordinates, model = matern)
  pred <- ~ log(1 - exp(-(distance / exp(df.lsigma))^-1)) + spat + Intercept + log(2)
  r <- lgcp(
    data = mexdolphin$points,
    samplers = cbind(mexdolphin$samplers, data.frame(weight = 1)),
    components = cmps,
    formula = pred,
    domain = list(
      coordinates = mexdolphin$mesh,
      distance = distance
    ),
    options = list(
      bru_verbose = 3,
      control.inla = list(int.strategy = "eb")
    )
  )
})
