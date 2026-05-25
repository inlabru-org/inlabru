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

test_that("bru_block_gcpo works for regular gaussian models", {
  skip_on_cran()
  local_bru_safe_inla()

  set.seed(123)
  n <- 50
  df <- data.frame(x = runif(n))
  df$y <- 2 + 3 * df$x + rnorm(n, sd = 0.5)

  gcpo_opts <- list(
    enable         = TRUE,
    type.cv        = "joint",
    num.level.sets = 2
  )

  # fit full model with auto-generated groups
  fit_full <- bru(
    y ~ Intercept(1) + x,
    family       = "gaussian",
    data         = df,
    control.gcpo = gcpo_opts
  )

  # basic structure
  expect_false(is.null(fit_full$gcpo$gcpo))
  expect_length(fit_full$gcpo$gcpo, n)

  # bru_block_gcpo_single and bru_block_gcpo give same result
  res_single <- bru_block_gcpo_single(fit_full)
  res <- bru_block_gcpo(fit_full)
  expect_identical(res_single, res)

  # correct structure
  expect_named(res, c("blocks", "gcpo"))
  expect_type(res$gcpo, "double")
  expect_length(res$gcpo, n)

  # scores on probability scale
  expect_true(all(res$gcpo > 0))

  # fit null model with fixed groups from full model
  fixed_groups <- fit_full$gcpo$groups
  fit_null <- bru(
    y ~ Intercept(1),
    family = "gaussian",
    data = df,
    control.gcpo = list(
      enable  = TRUE,
      type.cv = "joint",
      groups  = fixed_groups
    )
  )

  res_null <- bru_block_gcpo(fit_null)
  expect_named(res_null, c("blocks", "gcpo"))
  expect_length(res_null$gcpo, n)

  # bru_gcpo_table
  tbl <- bru_gcpo_table(full = fit_full, null = fit_null)

  expect_s3_class(tbl, "data.frame")
  expect_named(tbl, c("block", "full", "null"))
  expect_equal(nrow(tbl), n)

  # log scores and summary
  tbl[, -1] <- log(tbl[, -1])
  tbl$gcpo_summary <- rowSums(tbl[, -1])
  expect_equal(ncol(tbl), 4L)

  gcpo_summary <- data.frame(
    gcpo      = colSums(tbl[, c("full", "null")]),
    row.names = c("full", "null")
  )
  expect_named(gcpo_summary, "gcpo")
  expect_equal(nrow(gcpo_summary), 2L)

  # full model should have higher total log score
  expect_gt(gcpo_summary["full", "gcpo"], gcpo_summary["null", "gcpo"])
})

test_that("bru_block_gcpo returns correct structure for single likelihood with
          cp family", {
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
  a <- sf::st_intersects(nests, cvpart)
  if (!all(vapply(a, function(x) (length(x) == 1), logical(1)))) {
    stop("Point in none or multiple polygons")
  }
  nests$.block <- unlist(a)

  fit <- lgcp(
    components = geometry ~ Intercept(1),
    data = nests,
    samplers = cvpart,
    domain = list(geometry = gorillas_sf$mesh),
    control.gcpo = list(
      enable = TRUE,
      type.cv = "joint"
    )
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

test_that("bru_block_gcpo returns correct structure for multiple likelihoods
          with cp family", {
  skip_on_cran()
  local_bru_safe_inla()
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")

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

  # correct structure
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
  skip_if_not_installed("terra")

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
    components = geometry ~ Intercept(1),
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
    gcpo = colSums(log(df1[, c("intercept", "full")])),
    row.names = c("intercept", "full")
  )
  expect_named(gcpo_summary, "gcpo")
  expect_equal(nrow(gcpo_summary), 2L)
})

test_that("bru_block_gcpo works for mixed joint model (gaussian + cp)", {
  skip_on_cran()
  local_bru_safe_inla()
  skip_if_not_installed("sf")
  skip_if_not_installed("fmesher")

  set.seed(123)

  # unit square mesh
  loc.domain <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), ncol = 2, byrow = TRUE)
  mesh <- fmesher::fm_mesh_2d_inla(
    loc.domain = loc.domain,
    max.edge   = c(0.2, 0.3),
    offset     = c(0.1, 0.1)
  )

  # unit square boundary
  boundary <- sf::st_as_sf(sf::st_sfc(
    sf::st_polygon(list(
      matrix(c(0, 0, 1, 0, 1, 1, 0, 1, 0, 0), ncol = 2, byrow = TRUE)
    ))
  ))

  # block grid
  cvpart <- cv_hex(boundary, cellsize = 0.1, n_group = 1)
  cvpart$block_ID <- seq_len(nrow(cvpart))
  cvpart$group <- NULL
  nblock <- nrow(cvpart)

  # point pattern data
  n_cp <- 30
  pts <- sf::st_as_sf(
    data.frame(x = runif(n_cp), y = runif(n_cp)),
    coords = c("x", "y")
  )
  a <- sf::st_intersects(pts, cvpart)
  if (any(lengths(a) != 1)) stop("Point in none or multiple polygons")
  pts$.block <- unlist(a)

  # gaussian data
  n_g <- 20
  df_gauss <- data.frame(
    y = rnorm(n_g, mean = 2, sd = 0.5),
    x = runif(n_g)
  )

  # likelihoods
  lik_gauss <- bru_obs(
    y ~ Intercept,
    family       = "gaussian",
    data         = df_gauss,
    control.gcpo = list(enable = TRUE, type.cv = "joint"),
    tag          = "gauss"
  )
  lik_cp <- bru_obs(
    geometry ~ Intercept,
    family       = "cp",
    data         = pts,
    samplers     = cvpart,
    domain       = list(geometry = mesh),
    control.gcpo = list(enable = TRUE, type.cv = "joint"),
    tag          = "cp"
  )

  fit <- bru(~ Intercept(1), lik_gauss, lik_cp)

  # raw gcpo vector should exist
  expect_false(is.null(fit$gcpo$gcpo))

  # test bru_block_gcpo
  result <- bru_block_gcpo(fit)

  # correct top-level structure
  expect_named(result, c("blocks", "gcpo"))

  # two likelihoods
  expect_length(result$blocks, 2L)
  expect_named(result$blocks, c("gauss", "cp"))
  expect_named(result$gcpo, c("gauss", "cp"))

  # gauss: groups path - one score per observation
  expect_length(result$gcpo$gauss, n_g)

  # cp: block path - one score per spatial block
  expect_length(result$gcpo$cp, nblock)
})
