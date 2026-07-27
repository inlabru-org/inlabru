test_that("Expr mapper", {
  # testthat::skip("Expr mapper test disabled")
  mapper <-
    bm_expr(
      expr = rlang::quo(cos(x) + y[c(1, 1, 2)] + z),
      labels = list(root = "latent", derived = "effect", suffix = "_latent")
    )

  expect_equal(ibm_n(mapper), NA_integer_)
  state <- list(x = 1:3, y = 4:5)
  expect_equal(ibm_n(mapper, state = state), 5)
  expect_equal(
    ibm_n(mapper, state = state, multi = TRUE),
    c(x = 3, y = 2)
  )
  expect_null(ibm_values(mapper))

  dat <- data.frame(z = 11:13)

  val <- ibm_eval(
    mapper,
    input = list(),
    state = state,
    data = list(data = dat)
  )
  val_reference <- cos(state$x) + state$y[c(1, 1, 2)] + dat$z
  expect_equal(val, val_reference)

  A <- ibm_jacobian(
    mapper,
    input = list(),
    state = state,
    data = list(data = dat)
  )
  A_reference <- Matrix::sparseMatrix(
    i = c(1:3, 1:3),
    j = c(1:3, 4, 4, 5),
    x = c(-sin(state$x), rep(1, 3)),
    dims = c(3, 5)
  )
  expect_equal(A, A_reference, tolerance = lowtol)
})
