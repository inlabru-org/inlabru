local_bru_testthat_setup()

# Test for data input with mismatching input/output sizes (e.g. for regional
# integration models, etc)

test_that("Component construction: default mesh/mapping construction, data is list", {
  skip_on_cran()
  local_bru_safe_inla()

  lik <- like("gaussian",
              formula = y ~ .,
              data = list(x = c(1, 1.5, 2, 3, 4), y = 11:15),
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



test_that("Component construction: unsafe intercepts, data is list", {
  cmp <- component_list(~ something_unknown - 1)
  lik <- like(formula = response ~ ., data = list(response = 1:5))
  expect_warning(
    object = {
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
  skip_on_cran()
  local_bru_safe_inla()
  expect_warning(
    bru(~ something(map = a),
        formula = response ~ .,
        data = list(a = 1:5),
        options = list(bru_run = FALSE)
    ),
    "Use of 'map' is deprecated"
  )
})
