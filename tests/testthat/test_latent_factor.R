context("Latent models - factor (test_latent_factor.R)")

test_that("bru: factor component", {

  # Seed influences data as well as predict()!
  set.seed(123)

  # Factorial data
  input.df <- data.frame(
    x1 = factor(rep(c("Alpha", "Beta", "Mu"), c(70, 60, 70))),
    x2 = factor(rep(c("Delta", "Gamma", "H1", "H2", "H3", "H4", "H5"),
                    c(75, 75, 10, 10, 10, 10, 10)))
    )
  input.df$y <-
    c(1, -2, 0)[as.numeric(input.df$x1)] +
    c(1, -2, 1, -2, 1, -2, 1)[as.numeric(input.df$x2)] +
    rnorm(200, mean = 0, sd = 0.1)
  
  # Fit the model
#  fit0 <- inla(formula = y ~ 1 + x1 + f(x2, model = "iid"),
#             family = "gaussian",
#             data = input.df,
#             control.inla = list(int.strategy = "eb",
#                                 h = 0.005),
#             num.threads = "1:1"
#  )

  # Fit the model
  fit <- bru(components = ~ Intercept(main = 1) +
               fac1(main = x1, model = "factor") +
               fac2(main = x2, model = "iid"),
             formula = y ~ Intercept + fac1 + fac2,
             family = "gaussian",
             data = input.df,
             options = list(verbose = FALSE,
                            control.inla = list(int.strategy = "eb",
                                                h = 0.005),
                            num.threads = "1:1",
                            control.fixed = list(expand.factor.strategy = "inla"))
  )

#  fit0$summary.hyperpar
#  fit$summary.hyperpar
#  rbind(
#    cbind(ID = rownames(fit0$summary.fixed), fit0$summary.fixed),
#    fit0$summary.random$x2
#  )
#  rbind(
#    cbind(ID = rownames(fit$summary.fixed), fit$summary.fixed),
#    fit$summary.random$fac1,
#    fit$summary.random$fac2
#  )

  # Check factoreffect results
  expect_equal(fit$summary.random$fac1$mean,
               c(-3.072353, -1.092120),
               midtol)
  expect_equal(fit$summary.random$fac1$sd,
               c(0.04371715, 0.05016696),
               midtol)
  expect_equal(fit$summary.random$fac2$mean,
               c(1.213684, -1.723465,
                 1.310543, -1.671742,  1.302144, -1.700765,  1.271106),
               midtol)
  expect_equal(fit$summary.random$fac2$sd,
               c(0.5222028, 0.5208762,
                 0.5213029, 0.5213029, 0.5213029, 0.5213029, 0.5213029),
               midtol)
  
  # Check if prediction works
  pr <- predict(fit, data.frame(x1 = factor(c("Alpha", "Beta", "Mu")),
                                x2 = factor(c("Delta", "Delta", "Gamma"))),
                ~ fac1 + fac2, seed = 1)
  expect_equal(pr[, "mean"], c(1.168671, -1.904676, -2.857816), midtol)
  expect_equal(pr[, "sd"], c(0.5027916, 0.5041648, 0.5058944), midtol)
})
