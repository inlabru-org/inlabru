test_that("Mexdolphin: Hazard rate detection function", {
  skip_on_cran()
  local_bru_safe_inla()
  skip_if_not(bru_safe_sp())

  mexdolphin <- inlabru::mexdolphin_sp()

  sig <- function(x) bru_forward_transformation(qexp, x, rate = 1 / 8)
  hr <- function(distance, sigma) {
    1 - exp(-(distance / sigma)^-1)
  }
  log_hr <- function(distance, sigma) {
    log1p(-exp(-(distance / sigma)^-1))
  }
  cmp <- ~ sig_theta(1, prec.linear = 1) + Intercept(1)
  form <- distance ~ log_hr(distance, sig(sig_theta)) + Intercept
  form_list <- list(distance = distance) ~
    log(hr(distance, sig(sig_theta))) + Intercept

  pts <-
    mexdolphin$points

  fit <- bru(
    components = cmp,
    bru_obs(
      formula = form,
      family = "cp",
      data = pts,
      domain = list(distance = fm_mesh_1d(seq(0, 8, by = 0.1)))
    ),
    options = list(
      bru_verbose = 0,
      bru_compress_cp = TRUE,
      bru_max_iter = 10,
      verbose = FALSE,
      bru_initial = list(Intercept = 0, sig_theta = -1),
      inla.mode = "compact"
    )
  )

  expect_equal(
    fit$summary.fixed["sig_theta", "mean"],
    -0.4653159,
    tolerance = midtol
  )
  expect_equal(
    fit$summary.fixed["sig_theta", "sd"],
    0.3583088,
    tolerance = midtol
  )
  expect_equal(
    fit$summary.fixed["Intercept", "mean"],
    2.2593613,
    tolerance = midtol
  )
  expect_equal(
    fit$summary.fixed["Intercept", "sd"],
    0.2720904,
    tolerance = midtol
  )

  fit_list <- bru(
    components = cmp,
    bru_obs(
      formula = form_list,
      family = "cp",
      data = pts,
      domain = list(distance = fm_mesh_1d(seq(0, 8, by = 0.1)))
    ),
    options = list(
      bru_verbose = 0,
      bru_compress_cp = TRUE,
      bru_max_iter = 10,
      verbose = FALSE,
      bru_initial = list(Intercept = 0, sig_theta = -1),
      inla.mode = "compact"
    )
  )

  expect_equal(
    fit_list$summary.fixed["sig_theta", "mean"],
    -0.4653159,
    tolerance = midtol
  )
  expect_equal(
    fit_list$summary.fixed["sig_theta", "sd"],
    0.3583088,
    tolerance = midtol
  )
  expect_equal(
    fit_list$summary.fixed["Intercept", "mean"],
    2.2593613,
    tolerance = midtol
  )
  expect_equal(
    fit_list$summary.fixed["Intercept", "sd"],
    0.2720904,
    tolerance = midtol
  )

  ##########
  # pred <- predict(
  #   fit,
  #   newdata = data.frame(distance = seq(-8, 8, by = 0.01)),
  #   formula = ~ exp(Intercept + log_hr(abs(distance), sig(sig_theta)))
  # )
  #
  # library(ggplot2)
  # library(tidyverse)
  # ggplot(data.frame(distance = c(
  #   pts$distance,
  #   -pts$distance
  # ))) +
  #   geom_density(aes(distance, after_stat(count))) +
  #   geom_line(aes(distance, est, col = "Plugin"),
  #             data = data.frame(distance = seq(-8, 8, by = 0.001)) %>%
  #               mutate(est = hr(
  #                 abs(distance),
  #                 sig(fit$summary.fixed["sig_theta", "mean"])
  #               ) *
  #                 exp(fit$summary.fixed["Intercept", "mean"]))
  #   ) +
  #   geom_line(aes(distance, mean, col = "Pred"),
  #             data = pred
  #   ) +
  #   geom_line(aes(distance, est, col = "1"),
  #             data = data.frame(distance = seq(-8, 8, by = 0.001)) %>%
  #               mutate(est = hr(
  #                 abs(distance),
  #                 1
  #               ) *
  #                 exp(2.3))
  #   )
  #
  # pred <- predict(
  #   fit,
  #   data = ips,
  #   formula = ~ exp(Intercept + log_hr(abs(distance), sig(sig_theta)))
  # )
  # ggplot(data.frame(pts)) +
  #   stat_ecdf(aes(distance)) +
  #   geom_line(
  #     aes(distance,
  #         cumsum(hr(distance, sig(fit$summary.fixed["sig_theta", "mean"]))) /
  #           sum(hr(distance, sig(fit$summary.fixed["sig_theta", "mean"]))),
  #         col = "Plugin"
  #     ),
  #     data = data.frame(ips) %>% arrange(distance)
  #   ) +
  #   geom_line(
  #     aes(distance,
  #         cumsum(mean) / sum(mean),
  #         col = "Pred"
  #     ),
  #     data = pred %>% arrange(distance)
  #   ) +
  #   geom_line(
  #     aes(distance,
  #         cumsum(hr(distance, 1)) /
  #           sum(hr(distance, 1)),
  #         col = "1"
  #     ),
  #     data = data.frame(ips) %>% arrange(distance)
  #   )
  ##########
})


test_that("Marginal parameter transformation", {
  skip_on_cran()
  local_bru_safe_inla()
  data(mexdolphin_sf, package = "inlabru", envir = environment())

  hr <- function(distance, sigma) {
    1 - exp(-(distance / sigma)^-1)
  }
  log_hr <- function(distance, sigma) {
    log1p(-exp(-(distance / sigma)^-1))
  }
  cmp <- bru_component_list(~
    sigma(
      1,
      prec.linear = 1,
      marginal = bru_mapper_marginal(
        qfun = qexp,
        pfun = pexp,
        dfun = dexp,
        rate = 1 / 8
      )
    ) + Intercept(1))
  form <- distance ~ log_hr(distance, sigma = sigma) + Intercept

  pts <- mexdolphin_sf$points

  fit <- bru(
    components = cmp,
    bru_obs(
      formula = form,
      family = "cp",
      data = pts,
      domain = list(distance = fm_mesh_1d(seq(0, 8, by = 0.1)))
    ),
    options = list(
      bru_verbose = 0,
      bru_compress_cp = TRUE,
      bru_max_iter = 10,
      verbose = FALSE,
      bru_initial = list(
        Intercept = 0,
        sigma = ibm_eval(cmp$sigma$marginal, state = 1, inverse = TRUE)
      ),
      inla.mode = "compact"
    )
  )

  expect_equal(
    fit$summary.fixed["sigma", "mean"],
    -0.4653159,
    tolerance = midtol
  )
  expect_equal(
    fit$summary.fixed["sigma", "sd"],
    0.3583088,
    tolerance = midtol
  )
  expect_equal(
    fit$summary.fixed["Intercept", "mean"],
    2.2593613,
    tolerance = midtol
  )
  expect_equal(
    fit$summary.fixed["Intercept", "sd"],
    0.2720904,
    tolerance = midtol
  )
})
