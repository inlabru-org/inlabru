local_bru_testthat_setup()

test_that("bru: inla copy feature", {
  skip_on_cran()
  local_bru_safe_inla()

  # Seed influences data as well as predict()!
  set.seed(123L)

  df1 <- data.frame(x = cos(1:100))
  df2 <- data.frame(x = sin(1:100))
  df1 <- within(df1, y <- 1 + exp(2 * x) + rnorm(length(x), mean = 0, sd = 0.1))
  df2 <- within(df2, y <- 1 + (6 * x) + rnorm(length(x), mean = 0, sd = 0.1))

  cmp <- ~ +1 + myLin1(x,
    model = "rw1",
    mapper = bru_mapper(INLA::inla.mesh.1d(seq(-1, 1, length.out = 100)),
      indexed = FALSE
    )
  ) +
    myLin2(x, copy = "myLin1", fixed = FALSE)
  cmps <- component_list(cmp)

  fit <- bru(
    cmp,
    like(y ~ Intercept + exp(myLin1), family = "gaussian", data = df1, exclude = "myLin2"),
    like(y ~ Intercept + (myLin2), family = "gaussian", data = df2, exclude = "myLin1"),
    options = list(control.inla = list(int.strategy = "eb"))
  )


  expect_equal(
    fit$summary.fixed["Intercept", "mean"],
    1,
    tolerance = midtol
  )
  expect_equal(
    fit$summary.hyperpar["Beta for myLin2", "mean"],
    3,
    tolerance = midtol
  )

  pr <- predict(
    fit,
    data.frame(x = c(0.5, 1)),
    ~myLin2,
    n.samples = 500,
    seed = 1L
  )

  expect_equal(pr[, "mean"], c(3, 6), tolerance = midtol)
})
