test_that("fm_transform: geocentric globe transformation", {
  skip_on_cran()
  skip_if_not_installed("sf")
  skip_if_not_installed("sp")

  longlat_globe <- fm_crs("longlat_globe")
  globe <- fm_crs("globe")
  points0 <- sf::st_as_sf(
    as.data.frame(matrix(1:10, 5, 2)),
    coords = 1:2,
    crs = longlat_globe
  )

  points1 <- fm_transform(points0, crs = globe)
  points2 <- fm_transform(points1, crs = longlat_globe)

  expect_equal(
    range(rowSums(sf::st_coordinates(points1)^2)^0.5),
    6370.997 * c(1, 1)
  )
  expect_equal(
    sf::st_coordinates(points0),
    sf::st_coordinates(points2)
  )
})

test_that("fm_transform: geocentric globe transformation", {
  skip_on_cran()
  skip_if_not_installed("sf")
  skip_if_not_installed("sp")

  longlat_globe <- fm_crs("longlat_globe")
  longlat_sphere <- fm_crs("longlat_norm")
  globe <- fm_crs("globe")
  sphere <- fm_crs("sphere")
  points0 <- sf::st_as_sf(
    as.data.frame(matrix(1:10, 5, 2)),
    coords = 1:2,
    crs = longlat_globe
  )

  points1 <- fm_transform(points0, crs = sphere)
  points2 <- fm_transform(points1, crs = longlat_globe)

  expect_equal(
    range(rowSums(sf::st_coordinates(points1)^2)^0.5),
    c(1, 1)
  )
  expect_equal(
    sf::st_coordinates(points0),
    sf::st_coordinates(points2)
  )
})







test_that("fm_transform: geocentric globe transformation, sp involved", {
  skip_on_cran()
  skip_if_not_installed("sf")
  skip_if_not_installed("sp")

  longlat_globe <- fm_crs("longlat_globe")
  globe <- fm_crs("globe")
  #  warning(paste0(comment(fm_CRS(longlat_globe)), collapse = "\n"))
  points0_sp <-
    sp::SpatialPoints(
      matrix(1:10, 5, 2),
      proj4string = fm_CRS(longlat_globe)
    )
  points0 <- sf::st_as_sf(points0_sp)

  points1 <- fm_transform(points0, crs = globe)
  points2 <- fm_transform(points1, crs = longlat_globe)

  expect_equal(
    range(rowSums(sf::st_coordinates(points1)^2)^0.5),
    6370.997 * c(1, 1)
  )
  expect_equal(
    sf::st_coordinates(points0),
    sf::st_coordinates(points2)
  )
})

test_that("fm_transform: geocentric sphere transformation, sp involved", {
  skip_on_cran()
  skip_if_not_installed("sf")
  skip_if_not_installed("sp")

  longlat_globe <- fm_crs("longlat_globe")
  longlat_sphere <- fm_crs("longlat_norm")
  globe <- fm_crs("globe")
  sphere <- fm_crs("sphere")
  points0_sp <-
    sp::SpatialPoints(
      matrix(1:10, 5, 2),
      proj4string = fm_CRS(longlat_globe)
    )
  points0 <- sf::st_as_sf(points0_sp)

  points1 <- fm_transform(points0, crs = sphere)
  points2 <- fm_transform(points1, crs = longlat_globe)

  expect_equal(
    range(rowSums(sf::st_coordinates(points1)^2)^0.5),
    c(1, 1)
  )
  expect_equal(
    sf::st_coordinates(points0),
    sf::st_coordinates(points2)
  )
})
