local_bru_testthat_setup()

test_that("Multiple likelihoods: basic model", {
  skip_on_cran()
  local_bru_safe_inla()

  set.seed(123L)

  lik1 <- like("gaussian",
    formula = y ~ .,
    data = data.frame(
      x = rep(c(1, 1.5, 2, 3, 4), 2),
      y = rep(c(11, 12, 13, 14, 12), 2) + rnorm(10, sd = 0.05)
    ),
    include = c("int1", "effect"),
    # Checks that control.family is handled
    control.family = list(hyper = list(prec = list(fixed = TRUE)))
  )
  lik2 <- like("poisson",
    formula = y ~ .,
    data = data.frame(
      x = c(2, 2.5, 3, 4, 5),
      y = ceiling(exp(c(13, 13.5, 14, 12, 12) - 10))
    ),
    include = c("int2", "effect")
  )

  cmp1 <- component_list(~ effect(x, model = "rw2", scale.model = TRUE) - 1)
  cmp2 <- add_mappers(cmp1, lhoods = list(lik1, lik2))
  expect_equal(
    ibm_values(cmp2$effect$mapper, multi = 1)$main,
    sort(union(lik1$data$x, lik2$data$x))
  )

  cmp <- component_list(~ -1 +
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
