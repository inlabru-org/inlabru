basic_intercept_testdata <- function() {
  set.seed(123)
  data.frame(Intercept = 1,
             y = rnorm(100))
}

basic_fixed_effect_testdata <- function() {
  cbind(basic_intercept_testdata(),
        data.frame(x1 = rnorm(100)))
}
