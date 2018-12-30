
test_that("Component construction: linear model", {
  
  # Only check the following if the new backend is present
  if (exists("component")){
    
    cmp = component(~ beta(map = x, model = "linear"))[[1]]
    
    expect_equal(cmp$label, "beta")
    expect_equal(cmp$model, "linear")
    expect_equal(as.character(cmp$map), "x")
    
    # Covariate mapping
    df = data.frame(x = 1:10)
    expect_equal(map(cmp, data = df), 1:10)
    
    # Index generator
    idx = index(cmp, data = df)
    expect_is(idx, "data.frame")
    expect_equal(colnames(idx), "beta")
    expect_equal(idx$beta, 1:10)
    
    # A-matrix
    A = amatrix(cmp, data = df)
    expect_is(A, "ddiMatrix")
    expect_equal(nrow(A), 10)
    expect_equal(ncol(A), 10)
    expect_equal(diag(A), rep(1,10))
    
    # Value
    v = value(cmp, data = df, state = 2)
    expect_equal(v, 2 * df$x)
    
    v = value(cmp, data = df, state = list(beta = 2))
    expect_equal(v, 2 * df$x)
    
    v = value(cmp, data = df, state = list(beta = 2), A = A)
    expect_equal(v, 2 * df$x)
    
    # Inference
    set.seed(1)
    input.df <- data.frame(x=cos(1:100))
    input.df <- within(input.df, y <- 5 + 2*x + rnorm(100, mean=0, sd=0.1))
    
    # Fit a Gaussian likelihood model
    fit = bru(y ~ myLin(map = x, model = "linear") + Intercept, "gaussian", input.df)
    
    expect_equal(fit$summary.fixed["myLin", "mean"], 2.02, midtol)
    expect_equal(fit$summary.fixed["myLin", "sd"], 0.0126, midtol)
    
    pr = predict(fit, data.frame(x=c(1,2)), ~ myLin + 1)
  }
  
})
