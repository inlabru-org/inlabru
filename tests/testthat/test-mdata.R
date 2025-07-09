test_that("mdata", {
  local_bru_safe_inla()
  skip_if(utils::packageVersion("INLA") <= "24.06.02")

  sim.poisson <- function(prob, m) {
    stopifnot(length(prob) == length(m) && length(prob) > 0)
    n <- length(m)
    y <- numeric(n)
    event <- (runif(n) < prob)
    idx.zero <- which(event)
    idx.non.zero <- which(!event)
    y[idx.zero] <- 0
    y[idx.non.zero] <- rpois(length(idx.non.zero), lambda = m[idx.non.zero])
    y
  }
  ## chose link-function to use for the zero-inflation probability
  link.simple <- "logit"
  inv.link <- INLA::inla.link.invlogit
  ## link.simple <- "probit"
  ## inv.link <- inla.link.invprobit
  ## link.simple <- "cloglog"
  ## inv.link <- inla.link.invcloglog
  n <- 1000
  z <- rnorm(n, sd = 0.3)
  x <- rnorm(n, sd = 0.2)
  xx <- rnorm(n, sd = 0.3)
  zz <- rnorm(n, sd = 0.2)
  E <- runif(n, min = 0.8, max = 1 / 0.8)
  beta <- c(1, 1.1, 2.1, 0, -2, 1.2, 2.2, 0)
  eta2 <- beta[1] + beta[2] * xx + beta[3] * zz + beta[4] * xx * zz
  eta1 <- beta[5] + beta[6] * x + beta[7] * z + beta[8] * x * z
  prob <- inv.link(eta1)
  m <- E * exp(eta2)
  ok <- FALSE
  while (!ok) {
    y <- sim.poisson(prob, m)
    ok <- !all(y == 0)
  }
  ## head(data.frame(y, E, x, z, xx, zz))
  suppressWarnings(
    r <- INLA::inla(
      INLA::inla.mdata(
        cbind(y, E),
        cbind(1, x, z, x * z)
      ) ~ 1 + xx + zz + xx * zz,
      family = "0poisson",
      data = data.frame(y, E, x, z, xx, zz),
      control.fixed = list(prec = 1, prec.intercept = 1),
      control.compute = list(cpo = TRUE),
      control.family = list(
        link.simple = link.simple,
        hyper = list(
          beta1 = list(param = c(0, 1)),
          beta2 = list(param = c(0, 1)),
          beta3 = list(param = c(0, 1)),
          beta4 = list(param = c(0, 1)),
          beta5 = list(param = c(0, 1))
        )
      )
    )
  )
  #  y2 <- INLA::inla.mdata(cbind(y, -1, E/2),
  #                         cbind(1, x, z, x*z))
  #  r2 <- INLA::inla(
  #    y2 ~ 1 + xx + zz + xx*zz,
  #    family = "0poisson",
  #    data = data.frame(y, E, x, z, xx, zz),
  #    control.fixed = list(prec = 1, prec.intercept = 1),
  #    control.compute = list(cpo = TRUE),
  #    control.family = list(
  #      list(link.simple = link.simple,
  #                          hyper = list(beta1 = list(param = c(0, 1)),
  #                                       beta2 = list(param = c(0, 1)),
  #                                       beta3 = list(param = c(0, 1)),
  #                                       beta4 = list(param = c(0, 1)),
  #                                       beta5 = list(param = c(0, 1)))))
  #  )
  #

  r_bru <- bru(
    ~ 0 + fix(
      ~ 1 + xx + zz + xx * zz,
      model = "fixed",
      hyper = list(prec = list(initial = 0, fixed = TRUE))
    ),
    bru_obs(
      INLA::inla.mdata(cbind(y, E), cbind(1, x, z, x * z)) ~ .,
      family = "0poisson",
      data = data.frame(y, E, x, z, xx, zz),
      control.family = list(
        link.simple = link.simple,
        hyper = list(
          beta1 = list(param = c(0, 1)),
          beta2 = list(param = c(0, 1)),
          beta3 = list(param = c(0, 1)),
          beta4 = list(param = c(0, 1)),
          beta5 = list(param = c(0, 1))
        )
      )
    ),
    options = list(control.compute = list(cpo = TRUE))
  )

  expect_equal(
    r_bru$summary.fixed$mean,
    r$summary.random$fix$mean,
    tolerance = lowtol
  )
  expect_equal(
    r_bru$summary.fixed$sd,
    r$summary.random$fix$sd,
    tolerance = lowtol
  )


  if (FALSE) {
    rr <- inla(
      inla.mdata(cbind(y, E), cbind(1, xx, zz, xx * zz)) ~ 1 + x + z + x * z,
      family = "0poissonS",
      data = data.frame(y, E, x, z, xx, zz),
      control.fixed = list(prec = 1, prec.intercept = 1),
      control.compute = list(cpo = TRUE),
      ## in this case we need to define link.simple as the main link
      control.family = list(
        control.link = list(model = link.simple),
        hyper = list(
          beta1 = list(param = c(0, 1)),
          beta2 = list(param = c(0, 1)),
          beta3 = list(param = c(0, 1)),
          beta4 = list(param = c(0, 1)),
          beta5 = list(param = c(0, 1))
        )
      )
    )
    summary(r)
    summary(rr)
    res <- cbind(
      "beta" = beta,
      "0poisson" = c(r$summary.fixed$mean, r$summary.hyperpar$mean),
      "0poissonS" = c(rr$summary.hyperpar$mean, rr$summary.fixed$mean)
    )
    res <- cbind(res,
      diff = (res[, 2] - beta),
      diffS = (res[, 3] - beta),
      "diff/sd" = (res[, 2] - beta) /
        c(r$summary.fixed$sd, r$summary.hyperpar$sd),
      "diffS/sd" = (res[, 3] - beta) /
        c(rr$summary.hyperpar$sd, rr$summary.fixed$sd)
    )
    mm <- nrow(res) %/% 2
    rownames(res) <-
      c(paste0("beta", 1:mm, ".poisson"), paste0("beta", 1:mm, ".prob"))
    print(round(dig = 2, res))
  }
})



test_that("surv", {
  local_bru_safe_inla()
  skip_if(utils::packageVersion("INLA") <= "24.06.26")

  df <- data.frame(
    time = 1:4,
    event = c(1, 1, 0, 1),
    xx = rnorm(4),
    Intercept = 1
  )

  r_inla <- INLA::inla(
    INLA::inla.surv(time, event) ~ 0 +
      f(Intercept, model = "linear", prec.linear = 0) +
      f(xx, model = "linear"),
    family = "weibullsurv",
    data = df
  )

  r_bru <- bru(
    INLA::inla.surv(time, event) ~ 0 + Intercept(1, prec.linear = 0) + xx(xx),
    family = "weibullsurv",
    data = df
  )

  expect_equal(
    r_bru$summary.fixed[, "mean"],
    r_inla$summary.fixed[, "mean"],
    tolerance = lowtol
  )

  expect_equal(
    r_bru$summary.fixed[, "sd"],
    r_inla$summary.fixed[, "sd"],
    tolerance = midtol
  )

  expect_equal(
    r_bru$summary.hyperpar[, "mean"],
    r_inla$summary.hyperpar[, "mean"],
    tolerance = lowtol
  )

  expect_equal(
    r_bru$summary.hyperpar[, "sd"],
    r_inla$summary.hyperpar[, "sd"],
    tolerance = midtol
  )

  # Multi-family
  stk1.est <- INLA::inla.stack(
    data = list(.dummy = rep(1, 4), link = 1),
    A = list(1),
    effects = list(list(Intercept = 1, xx = df$xx)),
    responses = list(with(data = df, expr = INLA::inla.surv(time, event)))
  )
  stk1.pred <- INLA::inla.stack(
    data = list(.dummy = rep(1, 4), link = 1),
    A = list(1),
    effects = list(list(Intercept = 1, xx = df$xx)),
    responses = list(with(data = df, expr = INLA::inla.surv(time * NA, event)))
  )
  stk2 <- INLA::inla.stack(
    data = list(.dummy = rep(1, 4), link = 2),
    A = list(1),
    effects = list(list(Intercept = 1, xx = df$xx)),
    responses = list(with(data = df, expr = INLA::inla.surv(time, event)))
  )
  stk <- INLA::inla.stack(
    INLA::inla.stack(stk1.est, stk1.pred),
    stk2,
    multi.family = TRUE
  )


  stk_data <- INLA::inla.stack.data(stk, .response.name = "response")
  r2_inla <- INLA::inla(
    response ~ 0 + f(Intercept, model = "linear", prec.linear = 0) +
      f(xx, model = "linear"),
    family = c("exponentialsurv", "weibullsurv"),
    data = stk_data,
    control.predictor = list(A = INLA::inla.stack.A(stk), link = stk_data$link)
  )

  r2_bru <- bru(
    ~ 0 + Intercept(Intercept, prec.linear = 0) + xx(xx),
    bru_obs(
      INLA::inla.surv(time, event) ~ .,
      family = "exponentialsurv",
      data = rbind(
        df,
        data.frame(time = NA, event = df$event, xx = df$xx, Intercept = 1)
      )
    ),
    bru_obs(
      INLA::inla.surv(time, event) ~ .,
      family = "weibullsurv",
      data = df,
    )
  )

  expect_equal(
    r2_bru$summary.fixed[, c("mean", "sd")],
    r2_inla$summary.fixed[, c("mean", "sd")],
    tolerance = lowtol
  )

  expect_equal(
    r2_bru$summary.hyperpar[, c("mean", "sd")],
    r2_inla$summary.hyperpar[, c("mean", "sd")],
    tolerance = midtol
  )
})
