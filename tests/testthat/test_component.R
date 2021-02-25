local_bru_testthat_setup()

test_that("Component construction: linear model", {
  df <- data.frame(x = 1:10)

  cmp <- component_list(
    ~ beta(main = x, model = "linear", values = 1),
    lhoods = list(list(data = df))
  )[[1]]

  expect_equal(cmp$label, "beta")
  expect_equal(cmp$main$model, "linear")
  expect_equal(as.character(cmp$main$input$input), "x")

  # Covariate mapping
  df <- data.frame(x = 1:10)
  expect_equal(
    input_eval(cmp, data = df),
    list(
      main = 1:10,
      group = rep(1, 10),
      replicate = rep(1, 10)
    )
  )

  idx <- index_eval(cmp)
  expect_type(idx, "list")
  expect_equal(names(idx)[1], "beta")
  expect_equal(names(idx)[2], "beta.group")
  expect_equal(names(idx)[3], "beta.repl")
  expect_equal(idx$beta, 1)

  # A-matrix
  A <- amatrix_eval(cmp, data = df)
  expect_s4_class(A, "dgCMatrix")
  expect_equal(nrow(A), 10)
  expect_equal(ncol(A), 1)
  expect_equal(as.vector(A), 1:10)

  # Value
  v <- evaluate_effect_single(cmp, data = df, state = 2)
  expect_equal(v, 2 * df$x, ignore_attr = TRUE)

  v <- evaluate_effect_single(cmp, data = df, state = 2, A = A)
  expect_equal(v, 2 * df$x, ignore_attr = TRUE)

  cmps <- component_list(list(cmp))
  v <- evaluate_effect_multi(
    cmps,
    data = df,
    state = list(list(beta = 2))
  )
  expect_equal(v[[1]][["beta"]], 2 * df$x, ignore_attr = TRUE)

  v <- evaluate_effect_multi(
    cmps,
    data = NULL,
    state = list(list(beta = 2)),
    A = list(beta = A)
  )
  expect_equal(v[[1]][["beta"]], 2 * df$x, ignore_attr = TRUE)
})



test_that("Component construction: offset", {
  cmp <- component_list(~ something(a, model = "offset"))
  val <- evaluate_effect_single(cmp[["something"]],
    data = data.frame(a = 11:15),
    state = NULL
  )
  expect_equal(
    val,
    11:15,
    ignore_attr = TRUE
  )
})




test_that("Component construction: default mesh/mapping construction", {
  skip_on_cran()
  local_bru_safe_inla()
  
  lik <- like("gaussian",
    formula = y ~ .,
    data = data.frame(x = c(1, 1.5, 2, 3, 4), y = 11:15),
    include = "effect"
  )

  cmp1 <- component_list(~ effect(c(1, 1.5, 2, 3, 4), model = "iid") - 1)
  cmp2 <- add_mappers(cmp1, lhoods = list(lik))
  expect_equal(ibm_values(cmp2$effect$mapper, multi = 1)$main, lik$data$x)

  cmp1 <- component_list(~ effect(x, model = "rw2") - 1)
  cmp2 <- add_mappers(cmp1, lhoods = list(lik))
  expect_equal(ibm_values(cmp2$effect$mapper, multi = 1)$main, lik$data$x)

  mesh1 <- INLA::inla.mesh.1d(lik$data$x)
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
  cmp2 <- add_mappers(cmp1, lhoods = list(lik))
  expect_equal(ibm_values(cmp2$effect$mapper, multi = 1)$main, lik$data$x)

  cmp1 <- component_list(
    ~ effect(x,
      model = "rw2",
      mapper = bru_mapper(mesh1, indexed = TRUE)
    ) - 1
  )
  cmp2 <- add_mappers(cmp1, lhoods = list(lik))
  expect_equal(ibm_values(cmp2$effect$mapper, multi = 1)$main, seq_along(lik$data$x))
})



test_that("Component construction: unsafe intercepts", {
  cmp <- component_list(~ something_unknown - 1)
  lik <- like(formula = response ~ ., data = data.frame(response = 1:5))
  expect_warning(
    {
      model <- bru_model(cmp, list(lik))
    },
    "All covariate evaluations for 'something_unknown' are NULL"
  )
})

test_that("Component construction: deprecated arguments", {
  expect_warning(
    component_list(~ something(map = a)),
    "Use of 'map' is deprecated"
  )
  expect_warning(
    bru(~ something(map = a),
      formula = response ~ .,
      data = data.frame(a = 1:5),
      options = list(bru_run = FALSE)
    ),
    "Use of 'map' is deprecated"
  )
})
