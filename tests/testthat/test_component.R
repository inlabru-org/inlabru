context("Component construction (test_component.R)")
test_that("Component construction: linear model", {
  df <- data.frame(x = 1:10)

  cmp <- component(~ beta(main = x, model = "linear", values = 1),
    lhoods = list(data = df)
  )[[1]]

  expect_equal(cmp$label, "beta")
  expect_equal(cmp$main$model, "linear")
  expect_equal(as.character(cmp$main$input$input), "x")

  # Covariate mapping
  df <- data.frame(x = 1:10)
  expect_equal(input_eval(cmp, part = "main", data = df), 1:10)

  # TODO: The following must be preceeded by an auto-model initialisation based on

  # Index generator
  # likelihood data.
  idx <- index_eval(cmp)
  expect_is(idx, "list")
  expect_equal(names(idx)[1], "beta")
  expect_equal(names(idx)[2], "beta.group")
  expect_equal(names(idx)[3], "beta.repl")
  expect_equal(idx$beta, 1)

  # A-matrix
  A <- amatrix_eval(cmp, data = df)
  expect_is(A, "dgCMatrix")
  expect_equal(nrow(A), 10)
  expect_equal(ncol(A), 1)
  expect_equal(as.vector(A), 1:10)

  # Value
  v <- evaluate_effect_single(cmp, data = df, state = 2)
  expect_equivalent(v, 2 * df$x)

  v <- evaluate_effect_single(cmp, data = df, state = 2, A)
  expect_equivalent(v, 2 * df$x)

  cmps <- component(list(cmp))
  v <- evaluate_effect_multi(
    cmps,
    data = df,
    state = list(list(beta = 2))
  )
  expect_equivalent(v[[1]][["beta"]], 2 * df$x)

  v <- evaluate_effect_multi(
    cmps,
    data = NULL,
    state = list(list(beta = 2)),
    A = list(beta = A)
  )
  expect_equivalent(v[[1]][["beta"]], 2 * df$x)
})




# test_that("Component construction: default mesh/mapping construction", {
#
#   cmp <- list(label = "testlabel")
#   class(cmp) <- c("component", "list")
#
#   # Check for failure/success on valid/invalid inputs
#   expect_error(make.default.mesh(cmp,
#                                  model = NULL, model.type = "iid",
#                                  fvals = NULL
#   ))
#   expect_error(make.default.mesh(cmp,
#                                  model = NULL, model.type = "iid",
#                                  fvals = list()
#   ))
#   expect_error(
#     make.default.mesh(cmp,
#                       model = NULL, model.type = "iid",
#                       fvals = list(n = 2)
#     ),
#     NA
#   )
#
#   expect_error(make.default.mesh(cmp,
#                                  model = NULL, model.type = "seasonal",
#                                  fvals = NULL
#   ))
#   expect_error(make.default.mesh(cmp,
#                                  model = NULL, model.type = "seasonal",
#                                  fvals = list()
#   ))
#   expect_error(
#     make.default.mesh(cmp,
#                       model = NULL, model.type = "seasonal",
#                       fvals = list(season.length = 2)
#     ),
#     NA
#   )
#
#   for (model.type in c("rw1", "rw2", "ar", "ar1", "ou")) {
#     expect_error(make.default.mesh(cmp,
#                                    model = NULL, model.type = model.type,
#                                    fvals = NULL
#     ))
#     expect_error(make.default.mesh(cmp,
#                                    model = NULL, model.type = model.type,
#                                    fvals = list()
#     ))
#     expect_error(
#       make.default.mesh(cmp,
#                         model = NULL, model.type = model.type,
#                         fvals = list(values = 1:2)
#       ),
#       NA
#     )
#   }
# })
