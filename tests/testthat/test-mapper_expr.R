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

  val <- ibm_eval(mapper, input = list(data = dat), state = state)
  expect_equal(
    val,
    c(15.5403023058681403322, 15.583853163452858, 17.0100075033995571)
  )

  A <- ibm_jacobian(mapper, input = list(data = dat), state = state)

  A <- ibm_jacobian(mapper, input = list(data = as.list(dat)), state = state)

  A <- Matrix::sparseMatrix(
    i = 1:3,
    j = c(2, 7, 12),
    x = 1 * 2,
    dims = c(3, 12)
  )
  expect_equal(ibm_jacobian(mapper, olist_data), A)
})
