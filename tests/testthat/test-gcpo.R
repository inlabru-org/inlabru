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

test_that("bru_block_gcpo returns correct structure for single likelihood", {
  skip_on_cran()
  local_bru_safe_inla()
  skip_if_not_installed("sf")
  
  cvpart <- cv_hex(
    gorillas_sf$boundary,
    cellsize = 0.5,
    n_group = 3,
    resolution = c(95, 80)
  )
  cvpart$block_ID <- seq_len(nrow(cvpart))
  cvpart$group <- NULL
  nblock <- nrow(cvpart)
  
  nests <-gorillas_sf$nests
  a <- sf::st_intersects(nests, cvpart)
  if (!all(vapply(a, function(x) (length(x) == 1), logical(1)))) {
    stop("Point in none or multiple polygons")
  }
  nests$.block <- unlist(a)
  
  fit <- lgcp(
    components = geometry ~ Intercept(1) ,
    data = gorillas_nests_major,
    samplers = cvpart,
    domain = list(geometry = gorillas_sf$mesh),
    control.gcpo = list(enable = TRUE, 
                        type.cv = "joint")
  )
  
  result <- bru_block_gcpo(fit)
  
  # correct structure
  expect_type(result, "list")
  expect_named(result, c("blocks", "gcpo"))
  
  # one entry per block
  expect_length(result$blocks[[1]], nblock)
  expect_length(result$gcpo, nblock)
  
  # scores on probability scale
  expect_true(all(result$gcpo > 0))

})

test_that("bru_block_gcpo returns correct structure for multiple likelihoods", {
  skip_on_cran()
  local_bru_safe_inla()
  skip_if_not_installed("sf")
  
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
  
  elev <- gorillas_sf_gcov()$elevation
  elev <- elev - mean(terra::values(elev), na.rm = TRUE)
  f.elev <- function(where) {
    v <- eval_spatial(elev, where, layer = "elevation")
    v
  }
  
  cmp <- ~ elev(f.elev(.data.), model = "linear") +
    Intercept(1)
  
  fml.major <- geometry ~ Intercept 
  fml.minor <- geometry ~ Intercept + elev
  
  lik_major <- bru_obs("cp",
                       formula = fml.major,
                       samplers = cvpart,
                       data = gorillas_nests_major,
                       domain = list(geometry = gorillas_sf$mesh),
                       control.gcpo = list(
                         enable = TRUE,
                         type.cv = "joint"
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
                         type.cv = "joint"
                       ),
                       tag = "minor"
  )
  
  fit <- bru(cmp, lik_minor, lik_major)
  
  result <- bru_block_gcpo(fit)
  
  #correct structure
  expect_type(result, "list")
  expect_named(result, c("blocks", "gcpo"))
  expect_named(result$blocks, c("minor", "major"))
  expect_named(result$gcpo, c("minor", "major"))
  
  # correct length entry per block
  expect_length(result$blocks[[1]], nblock)
  expect_length(result$gcpo[[1]], nblock)
})

test_that("bru_gcpo_table returns correct structure for multiple fits", {
  skip_on_cran()
  local_bru_safe_inla()
  skip_if_not_installed("sf")
  
  cvpart <- cv_hex(
    gorillas_sf$boundary,
    cellsize = 0.5,
    n_group = 3,
    resolution = c(95, 80)
  )
  cvpart$block_ID <- seq_len(nrow(cvpart))
  cvpart$group <- NULL
  nblock <- nrow(cvpart)
  
  nests <- gorillas_sf$nests
  
  # blocks 
  a <- sf::st_intersects(nests, cvpart)
  if (!all(vapply(a, function(x) (length(x) == 1), logical(1)))) {
    stop("Point in none or multiple polygons")
  }
  nests$.block <- unlist(a)
  
  # covariate
  elev <- gorillas_sf_gcov()$elevation
  elev <- elev - mean(terra::values(elev), na.rm = TRUE)
  f.elev <- function(where) {
    v <- eval_spatial(elev, where, layer = "elevation")
    v
  }
  
  # models
  fit1 <- lgcp(
    components =geometry ~ Intercept(1), 
    data = nests,
    samplers = cvpart,
    domain = list(geometry = gorillas_sf$mesh),
    control.gcpo = list(enable = TRUE, type.cv = "joint")
  )
  
  fit2 <- lgcp(
    components = geometry ~ Intercept(1) + 
      elev(f.elev(.data.), model = "linear"),
    data = nests,
    samplers = cvpart,
    domain = list(geometry = gorillas_sf$mesh),
    control.gcpo = list(enable = TRUE, type.cv = "joint")
  )
  
  # via named list
  df1 <- bru_gcpo_table(list(intercept = fit1, full = fit2))
  expect_s3_class(df1, "data.frame")
  expect_named(df1, c("block", "intercept", "full"))
  expect_equal(nrow(df1), nblock)
  
  # via dots
  df2 <- bru_gcpo_table(intercept = fit1, full = fit2)
  expect_identical(df1, df2)
  
  # via mixed
  df3 <- bru_gcpo_table(list(intercept = fit1), full = fit2)
  expect_identical(df1, df3)
  
  # summary
  gcpo_summary <- data.frame(
    gcpo = colSums(log(df1[, c("intercept","full")])),
    row.names = c("intercept","full")
  )
  expect_named(gcpo_summary, "gcpo")
  expect_equal(nrow(gcpo_summary), 2L)
})
