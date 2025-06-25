test_that("Additivity: Additive predictor with nonlinear components", {
  skip_on_cran()
  local_bru_safe_inla()

  withr::local_seed(12345L)

  df <- tibble::tibble(
    u = 2,
    x = qexp(pnorm(u)),
    y = rpois(10, exp(x))
  )

  # Linear
  cmp <- ~ x(1, mean.linear = 0, prec.linear = 1) + 0
  form <- y ~ x
  fit <- bru(
    components = cmp,
    bru_obs(
      formula = form,
      family = "poisson",
      data = df
    )
  )

  # The latent variable should match x
  expect_equal(
    fit$summary.fixed["x", "mean"],
    3.821398,
    tolerance = midtol
  )

  # Nonlinear
  cmp <- ~ x(
    1,
    mean.linear = 0,
    prec.linear = 1,
    marginal = bm_marginal(qexp, pexp, dexp, rate = 1)
  ) + 0
  form <- y ~ x
  fit <- bru(
    components = cmp,
    bru_obs(
      formula = form,
      family = "poisson",
      data = df
    )
  )

  # The latent variable should match u
  expect_equal(
    fit$summary.fixed["x", "mean"],
    2.007725,
    tolerance = midtol
  )
})
