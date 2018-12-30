
test_that("bru: linear component", {
  
  # Seed influences data as well as predict()!
  set.seed(123)
  
  input.df <- data.frame(x=cos(1:100))
  input.df <- within(input.df, y <- 5 + 2*x + rnorm(100, mean=0, sd=0.1))

  fit = bru(y ~ myLin(map = x, model = "linear") + Intercept, "gaussian", input.df)
  
  expect_equal(fit$summary.fixed["myLin", "mean"], 2.02, midtol)
  expect_equal(fit$summary.fixed["myLin", "sd"], 0.0126, midtol)
  
  pr = predict(fit, data.frame(x=c(1,2)), ~ myLin + 1, seed = 1)
  
  expect_equal(pr[,"mean"], c(3.019780, 5.039559), midtol)
  
})


test_that("bru: factor component", {
  
  # Seed influences data as well as predict()!
  set.seed(123)
  
  # Factorial data
  input.df <- data.frame(x = c(rep("Alpha", 10), rep("Beta", 10)))
  input.df$y <- c(rep(1, 10), rep(-2, 10)) + rnorm(20, mean=0, sd=0.1)
  
  # Fit the model
  fit = bru(y ~ fac(map = x, model = "factor") - Intercept, "gaussian", input.df)
  
  # Check fixed effect results
  expect_equal(fit$summary.fixed["facAlpha", "mean"], 1.007462, midtol)
  expect_equal(fit$summary.fixed["facAlpha", "sd"], 0.03185, midtol)
  expect_equal(fit$summary.fixed["facBeta", "mean"], -1.979136, midtol)
  expect_equal(fit$summary.fixed["facBeta", "sd"], 0.03185, midtol)
  
  # Check if prediction works
  pr = predict(fit, data.frame(x=c("Alpha","Beta")), ~ fac + 1, seed = 1)
  expect_equal(pr[, "mean"], c(2.00522, -0.97751), midtol)
  
})
  