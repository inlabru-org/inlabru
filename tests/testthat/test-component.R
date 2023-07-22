test_that("Component construction: linear model", {
  df <- data.frame(x = 1:10, response = 1:10)

  llik <- like_list(list(like(formula = response ~ ., data = df)))

  # Using label as input:
  cmp0 <- component_list(
    list(component("x", x)),
    lhoods = llik
  )[["x"]]

  expect_equal(cmp0$label, "x")
  expect_equal(cmp0$main$model, "linear")
  expect_equal(as.character(cmp0$main$input$input), "x")

  # Using label as input:
  cmp0 <- component_list(
    ~x,
    lhoods = llik
  )[["x"]]

  expect_equal(cmp0$label, "x")
  expect_equal(cmp0$main$model, "linear")
  expect_equal(as.character(cmp0$main$input$input), "x")

  cmp <- component_list(
    ~ beta(main = x, model = "linear", values = 1),
    lhoods = llik
  )[["beta"]]

  expect_equal(cmp$label, "beta")
  expect_equal(cmp$main$model, "linear")
  expect_equal(as.character(cmp$main$input$input), "x")

  # Covariate mapping
  df <- data.frame(x = 1:10)
  inp <- input_eval(cmp, data = df)
  expect_equal(
    inp,
    list(
      mapper = list(
        main = 1:10,
        group = 1,
        replicate = 1
      ),
      scale = NULL
    )
  )

  idx <- index_eval(cmp, inla_f = FALSE)
  expect_type(idx, "list")
  expect_equal(names(idx)[1], "beta")
  expect_equal(names(idx)[2], "beta.group")
  expect_equal(names(idx)[3], "beta.repl")
  expect_equal(idx$beta, 1)

  # A-matrix
  comp_lin <- comp_lin_eval(cmp, input = inp, inla_f = FALSE)
  A <- ibm_jacobian(comp_lin)
  expect_s4_class(A, "dgCMatrix")
  expect_equal(nrow(A), 10)
  expect_equal(ncol(A), 1)
  expect_equal(as.vector(A), 1:10)

  # Value
  v <- evaluate_effect_single_state(comp_lin, state = 2)
  expect_equal(v, 2 * df$x, ignore_attr = TRUE)

  v <- evaluate_effect_single_state(comp_lin, state = 2, A = A)
  expect_equal(v, 2 * df$x, ignore_attr = TRUE)

  cmps <- component_list(list(cmp))
  inps <- input_eval(cmps, data = df)
  v <- evaluate_effect_multi_state(
    cmps,
    input = inps,
    state = list(list(beta = 2))
  )
  expect_equal(v[[1]][["beta"]], 2 * df$x, ignore_attr = TRUE)

  v <- evaluate_effect_multi_state(
    cmps,
    input = inps,
    state = list(list(beta = 2))
  )
  expect_equal(v[[1]][["beta"]], 2 * df$x, ignore_attr = TRUE)
})



test_that("Component construction: duplicate detection", {
  expect_error(
    component_list(
      ~ -1 +
        beta(main = x, model = "linear", values = 1) +
        beta(main = x, model = "linear", values = 2)
    ),
    regexp = "Duplicated component labels detected: 'beta'"
  )
})



test_that("Component construction: offset", {
  cmp <- component_list(~ -1 + something(a, model = "offset"))
  inp <- input_eval(cmp, data = data.frame(a = 11:15))
  val <- evaluate_effect_single_state(cmp,
    input = inp,
    state = NULL
  )
  expect_equal(
    val$something,
    11:15,
    ignore_attr = TRUE
  )
})



test_that("Component construction: terra", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")

  f <- system.file("ex/elev.tif", package = "terra")
  r <- terra::rast(f)

  data <- sf::st_sf(
    data.frame(
      geometry = sf::st_sfc(
        sf::st_point(cbind(6, 50)),
        crs = "epsg:4326"
      ),
      response = 1
    )
  )

  expect_equal(eval_spatial(r, data), 406)

  llik <- like_list(list(like(formula = response ~ ., data = data)))

  cmp <- component_list(
    ~ -1 + something(eval_spatial(r, geometry), model = "linear"),
    lhoods = llik
  )
  inp <- input_eval(cmp, data = data)
  comp_lin <- comp_lin_eval(cmp,
    input = inp,
    state = list(something = 2),
    inla_f = FALSE
  )
  val <- evaluate_effect_single_state(
    comp_lin,
    state = list(something = 2)
  )
  expect_equal(
    val$something,
    2 * 406
  )

  cmp <- component_list(~ -1 + something(r, model = "linear"),
    lhoods = llik
  )
  inp <- input_eval(cmp, data = data)
  comp_lin <- comp_lin_eval(cmp,
    input = inp,
    state = list(something = 2),
    inla_f = FALSE
  )
  val <- evaluate_effect_single_state(
    comp_lin,
    state = list(something = 2)
  )
  expect_equal(
    val$something,
    2 * 406
  )

  cmp <- component_list(~ -1 + something(r, model = "linear", main_layer = 1),
    lhoods = llik
  )
  inp <- input_eval(cmp, data = data)
  comp_lin <- comp_lin_eval(cmp,
    input = inp,
    state = list(something = 2),
    inla_f = FALSE
  )
  val <- evaluate_effect_single_state(
    comp_lin,
    state = list(something = 2)
  )
  expect_equal(
    val$something,
    2 * 406
  )

  cmp <- component_list(~ -1 + something(r, model = "linear", main_layer = "elevation"),
    lhoods = llik
  )
  inp <- input_eval(cmp, data = data)
  comp_lin <- comp_lin_eval(cmp,
    input = inp,
    state = list(something = 2),
    inla_f = FALSE
  )
  val <- evaluate_effect_single_state(
    comp_lin,
    state = list(something = 2)
  )
  expect_equal(
    val$something,
    2 * 406
  )

  expect_error(
    component_list(~ something(r, model = "linear", main_layer = 2),
      lhoods = llik
    ),
    NULL
  )

  expect_error(
    component_list(~ something(r, model = "linear", main_layer = "elev"),
      lhoods = llik
    ),
    NULL
  )
})




test_that("Component construction: default index/mesh/mapping construction", {
  skip_on_cran()
  local_bru_safe_inla()

  lik <- like("gaussian",
    formula = y ~ .,
    data = data.frame(x = c(1, 1.5, 2, NA, 4), y = 11:15),
    include = "effect"
  )

  cmp1 <- component_list(~ effect(c(1, 1.5, 2, NA, 4), model = "iid") - 1)
  cmp2 <- add_mappers(cmp1, lhoods = like_list(list(lik)))
  expect_equal(
    ibm_values(cmp2$effect$mapper, multi = 1)$main,
    sort(unique(lik$data$x), na.last = NA)
  )
  expect_equal(
    ibm_eval(cmp2$effect$mapper$mappers$mapper$mappers$main,
      input = c(1, NA, 4),
      state = c(11, 12, 13, 14)
    ),
    c(11, 0, 14)
  )

  cmp1 <- component_list(~ effect(x, model = "rw2") - 1)
  cmp2 <- add_mappers(cmp1, lhoods = like_list(list(lik)))
  expect_equal(
    ibm_values(cmp2$effect$mapper, multi = 1)$main,
    sort(unique(lik$data$x), na.last = NA)
  )

  mesh1 <- fm_mesh_1d(
    sort(unique(lik$data$x), na.last = NA)
  )
  expect_error(
    component_list(
      ~ effect(x, model = "rw2", mapper = mesh1) - 1
    ),
    regexp = "Unknown mapper"
  )

  cmp1 <- component_list(
    ~ effect(x,
      model = "rw2",
      mapper = bru_mapper(mesh1, indexed = FALSE)
    ) - 1
  )
  cmp2 <- add_mappers(cmp1, lhoods = like_list(list(lik)))
  expect_equal(
    ibm_values(cmp2$effect$mapper, multi = 1)$main,
    sort(unique(lik$data$x), na.last = NA)
  )

  cmp1 <- component_list(
    ~ effect(x,
      model = "rw2",
      mapper = bru_mapper(mesh1, indexed = TRUE)
    ) - 1
  )
  cmp2 <- add_mappers(cmp1, lhoods = like_list(list(lik)))
  expect_equal(
    ibm_values(cmp2$effect$mapper, multi = 1)$main,
    seq_along(sort(unique(lik$data$x), na.last = NA))
  )
})



test_that("Component construction: main iid factor construction", {
  skip_on_cran()

  lik <- like("gaussian",
    formula = y ~ .,
    data = data.frame(x = as.factor(c(1, 1.5, 2, 3, 4)), y = 11:15),
    include = "effect"
  )

  cmp1 <- component_list(~ effect(as.factor(c(1, 1.5, 2, 3, 4)),
    model = "iid"
  ) - 1)
  cmp2 <- add_mappers(cmp1, lhoods = like_list(list(lik)))
  expect_equal(
    ibm_values(cmp2$effect$mapper, multi = 1)$main,
    as.character(lik$data$x)
  )

  expect_equal(
    ibm_eval(cmp2$effect$mapper$mappers$mapper$mappers$main,
      input = as.factor(c(1, NA, 4)),
      state = c(11, 12, 13, 14, 15)
    ),
    c(11, 0, 15)
  )
})

test_that("Component construction: group iid factor construction", {
  skip_on_cran()

  lik <- like("gaussian",
    formula = y ~ .,
    data = data.frame(x = as.factor(c(1, 1.5, 2, 3, 4)), y = 11:15),
    include = "effect"
  )

  cmp1 <- component_list(
    ~ -1 +
      effect(rep(1, 5),
        group = x,
        model = "iid",
        control.group = list(model = "iid")
      )
  )
  cmp2 <- add_mappers(cmp1, lhoods = like_list(list(lik)))
  expect_equal(
    ibm_values(cmp2$effect$mapper, multi = 1)$group,
    as.numeric(lik$data$x)
  )

  expect_equal(
    ibm_eval(cmp2$effect$mapper$mappers$mapper$mappers$group,
      input = as.factor(c(1, NA, 4)),
      state = c(11, 12, 13, 14, 15)
    ),
    c(11, 0, 15)
  )

  local_bru_safe_inla()
  expect_no_error(bru(cmp2, lik))
})

test_that("Component construction: replicate iid factor construction", {
  skip_on_cran()

  lik <- like("gaussian",
    formula = y ~ .,
    data = data.frame(x = as.factor(c(1, 1.5, 2, 3, 4)), y = 11:15),
    include = "effect"
  )

  cmp1 <- component_list(
    ~ -1 +
      effect(rep(1, 5),
        replicate = x,
        model = "iid"
      )
  )
  cmp2 <- add_mappers(cmp1, lhoods = like_list(list(lik)))
  expect_equal(
    ibm_values(cmp2$effect$mapper, multi = 1)$replicate,
    as.numeric(lik$data$x)
  )

  expect_equal(
    ibm_eval(cmp2$effect$mapper$mappers$mapper$mappers$replicate,
      input = as.factor(c(1, NA, 4)),
      state = c(11, 12, 13, 14, 15)
    ),
    c(11, 0, 15)
  )

  local_bru_safe_inla()
  expect_no_error(bru(cmp2, lik))
})



test_that("Component construction: unsafe intercepts", {
  cmp <- component_list(~ something_unknown - 1)
  lik <- like(formula = response ~ ., data = data.frame(response = 1:5))
  expect_warning(
    object = {
      model <- bru_model(cmp, like_list(list(lik)))
    },
    "All covariate evaluations for 'something_unknown' are NULL"
  )
})

test_that("Component construction: deprecated arguments", {
  expect_warning(
    component_list(~ something(map = a)),
    "Use of 'map' is deprecated"
  )
  skip_on_cran()
  local_bru_safe_inla()
  expect_warning(
    bru(~ something(map = a),
      formula = response ~ .,
      data = data.frame(a = 1:5, response = 11:15),
      options = list(bru_run = FALSE)
    ),
    "Use of 'map' is deprecated"
  )
})
