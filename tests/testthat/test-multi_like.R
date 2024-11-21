test_that("Multiple likelihoods: basic model", {
  skip_on_cran()
  local_bru_safe_inla()

  set.seed(123L)

  lik1 <- bru_obs("gaussian",
    formula = y ~ .,
    data = data.frame(
      x = rep(c(1, 1.5, 2, 3, 4), 2),
      y = rep(c(11, 12, 13, 14, 12), 2) + rnorm(10, sd = 0.05)
    ),
    include = c("int1", "effect"),
    # Checks that control.family is handled
    control.family = list(hyper = list(prec = list(fixed = TRUE)))
  )
  lik2 <- bru_obs("poisson",
    formula = y ~ .,
    data = data.frame(
      x = c(2, 2.5, 3, 4, 5),
      y = ceiling(exp(c(13, 13.5, 14, 12, 12) - 10))
    ),
    include = c("int2", "effect")
  )

  cmp1 <- bru_component_list(~ effect(x, model = "rw2", scale.model = TRUE) - 1)
  cmp2 <- add_mappers(cmp1, lhoods = list(lik1, lik2))
  expect_equal(
    ibm_values(cmp2$effect$mapper, multi = 1)$main,
    sort(union(lik1$data$x, lik2$data$x))
  )

  cmp <- bru_component_list(~ -1 +
    effect(x,
      model = "rw2",
      values = seq(1, 5, by = 0.25),
      scale.model = TRUE
    ) +
    int1(1) + int2(1))

  fit <- bru(cmp, lik1, lik2)

  expect_equal(nrow(fit$summary.hyperpar), 1)
  expect_equal(fit$summary.hyperpar["Precision for effect", "mean"],
    2.0459,
    tolerance = midtol
  )
})


test_that("Predictor indexing", {
  skip_on_cran()
  local_bru_safe_inla()

  fit <- bru(
    ~ 0 + x,
    bru_obs(
      y ~ .,
      data = data.frame(x = 1:3, y = 1:3 + rnorm(3)),
      tag = "A"
    ),
    bru_obs(
      y ~ .,
      data = data.frame(x = 1:4, y = c(NA, NA, 3:4) + rnorm(4)),
      tag = "B"
    )
  )

  expect_equal(
    bru_index(fit),
    1L:7L
  )
  expect_equal(
    bru_index(fit, "A"),
    1L:3L
  )
  expect_equal(
    bru_index(fit, "B"),
    4L:7L
  )
  expect_equal(
    bru_index(fit, c("B", "A")),
    c(4L:7L, 1L:3L)
  )
  expect_equal(
    bru_index(fit, what = "observed"),
    c(1L:3L, 6L:7L)
  )
  expect_equal(
    bru_index(fit, what = "missing"),
    4L:5L
  )
})
