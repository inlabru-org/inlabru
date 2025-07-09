test_that("Response and predictor mismatch handling", {
  skip_on_cran()
  local_bru_safe_inla()
  n <- 10
  df <- data.frame(Y = rpois(n, lambda = 3))

  cmpA <- ~ -1 + beta(1, model = "linear")
  cmpB <- ~ -1 + beta(rep(1, NROW(.data.)), model = "linear")

  lik1 <- bru_obs("poisson",
    formula = Y ~ .,
    data = df
  )

  expect_no_error({
    fit1A <- bru(components = cmpA, lik1)
  })

  expect_no_error({
    fit1B <- bru(components = cmpB, lik1)
  })

  lik2 <- bru_obs("poisson",
    formula = Y ~ beta,
    data = df
  )

  expect_no_error({
    fit2A <- bru(components = cmpA, lik2)
  })
  expect_no_error({
    fit2B <- bru(components = cmpB, lik2)
  })

  lik3 <- bru_obs("poisson",
    formula = Y ~ c(beta, beta),
    data = df
  )

  expect_error(
    {
      fit3A <- bru(components = cmpA, lik3)
    },
    paste0(
      "Number of rows \\(2\\) in the predictor for component 'beta' ",
      "does not match the length implied by the response data \\(10\\)"
    )
  )
  expect_error(
    {
      fit3B <- bru(components = cmpB, lik3)
    },
    paste0(
      "Number of rows \\(20\\) in the predictor for component 'beta' ",
      "does not match the length implied by the response data \\(10\\)"
    )
  )


  # INLA prior to mdata/surv stack support will give a different error message
  # than the one below, so we skip the test for old INLA versions.
  skip_if(utils::packageVersion("INLA") <= "24.06.26")

  expect_error(
    {
      fit4 <- bru(y ~ 0 + comp(1:3),
        data = data.frame(y = rnorm(2)),
        family = "gaussian"
      )
    },
    paste0(
      "The total number of response values \\(N=2\\) and predictor ",
      "values \\(N=3\\) do not match.\n",
      "  This is likely due to a mistake in the component or predictor ",
      "constructions."
    )
  )
})


test_that("Complex list data handling", {
  skip_on_cran()
  local_bru_safe_inla()
  withr::local_seed(12345L)
  n <- 6
  m <- 3
  data <- list(
    x1 = rnorm(m * n),
    weight = runif(m * n),
    .block = rep(1:n, m),
    x2 = rnorm(n)
  )
  agg <- bm_logsumexp()
  resp_data <- with(data, data.frame(Y = rpois(n, lambda = exp(
    1 +
      ibm_eval(
        agg,
        input = list(weight = weight, block = .block),
        state = x1
      ) + x2
  ))))
  cmpA <- ~ 0 + Intercept(1) + x1 + x2

  lik1 <- bru_obs(
    "poisson",
    formula = Y ~ Intercept + ibm_eval(
      agg,
      input = list(weight = weight, block = .block),
      state = x1
    ) + x2,
    data = data,
    allow_combine = TRUE,
    response_data = resp_data
  )

  expect_no_error({
    fit <- bru(components = cmpA, lik1)
  })
})
