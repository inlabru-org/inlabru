context("Latent models - AR1 (test_latent_ar1.R)")

latent_ar1_testdata <- function() {
  data1 <- data.frame(time = c(1,2,3,5,4),
                      obs = c(1,2,1,4,2))
  data1 <- data.frame(time = c(1,2,3.5,5,4),
                      obs = c(1,2,1,4,2))
  data2 <- data1[5:1, , drop = FALSE]
  data3 <- data1[c(2,5,4,3,1), , drop = FALSE]
  
  cmp <- obs ~ time(time, model = "ar1") - Intercept
  fit1 <- bru(cmp,
              data = data1, family = "gaussian"
  )
  fit2 <- bru(cmp,
              data = data2, family = "gaussian"
  )
  fit3 <- bru(cmp,
              data = data3, family = "gaussian"
  )
  
  list(
    data = list(data1 = data1, data2 = data2, data3 = data3),
    cmp = cmp,
    fit = list(fit1 = fit1, fit2 = fit2, fit3 = fit3)
  )
}

fit <- latent_ar1_testdata()

cbind(
  fit$fit[[1]]$summary.random$time$mean,
  fit$fit[[2]]$summary.random$time$mean,
  fit$fit[[3]]$summary.random$time$mean
)
cbind(
  fit$fit[[1]]$summary.random$time$sd,
  fit$fit[[2]]$summary.random$time$sd,
  fit$fit[[3]]$summary.random$time$sd
)
cbind(
  fit$fit[[1]]$summary.hyperpar[,"mean"],
  fit$fit[[2]]$summary.hyperpar[,"mean"],
  fit$fit[[3]]$summary.hyperpar[,"mean"]
)

fit$fit_inla <- list(
  INLA::inla(obs ~ f(time, model="ar1") - 1,
             data = fit$data[[1]], family = "gaussian"),
  INLA::inla(obs ~ f(time, model="ar1") - 1,
             data = fit$data[[2]], family = "gaussian"),
  INLA::inla(obs ~ f(time, model="ar1") - 1,
             data = fit$data[[3]], family = "gaussian")
)

cbind(
  fit$fit_inla[[1]]$summary.hyperpar[,"mean"],
  fit$fit_inla[[2]]$summary.hyperpar[,"mean"],
  fit$fit_inla[[3]]$summary.hyperpar[,"mean"]
)
cbind(
  fit$fit_inla[[1]]$summary.random$time$mean,
  fit$fit_inla[[2]]$summary.random$time$mean,
  fit$fit_inla[[3]]$summary.random$time$mean
)


test_that("Latent models: AR1 bru ordering", {
  skip_on_cran()
  library(INLA)
  data <- latent_ar1_testdata()
  
  # Check Intercept
  expect_equal(data$fit[[1]]$summary.fixed["Intercept", "mean"], 5.982867, midtol)
  
  # Check SPDE
  expect_equal(
    data$fit$summary.random$field$mean[c(1, 25, 50)],
    c(-5.195488, -4.480026, -7.150024), midtol
  )
})
