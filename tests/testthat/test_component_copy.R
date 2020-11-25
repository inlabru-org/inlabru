local_bru_testthat_setup()

test_that("Component copy feature", {
  skip_on_cran()
  local_bru_safe_inla()
  
  mydata <- data.frame(x1 = rep(1:4, times = 2),
                       x2 = rep(c(1,2), each = 4))
  mydata <- within(mydata, y <- rpois(8, exp(x1^0.5 + x2^0.5*2 - 1)))
  
  cmp <- y ~ -1 + xx(x1, model = "rw2") + xx(x2, copy = "xx")

  inlaform <- y ~ 1 +
    f(x1, model = "rw2", values=1:4, scale.model = TRUE) +
    f(x2, copy = "x1", fixed = FALSE)
  fit <- INLA::inla(formula = inlaform, data = mydata, family = "poisson")
  
})
