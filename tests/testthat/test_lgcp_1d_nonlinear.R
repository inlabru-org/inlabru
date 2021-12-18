local_bru_testthat_setup()

test_that("Mexdolphin: Hazard rate detection function", {
  skip_on_cran()
  local_bru_safe_inla()
  data(mexdolphin, package = "inlabru", envir = environment())

  hr <- function(distance, lsig) {
    1 - exp(-(distance / (exp(lsig)))^-1)
  }

  cmp <- ~ lsig(1) + Intercept(1)
  form <- distance ~ log(hr(distance, lsig)) + Intercept
  ips <- ipoints(INLA::inla.mesh.1d(seq(0, 8, by = 0.1)), name = "distance")

  fit <- lgcp(
    components = cmp,
    mexdolphin$points,
    ips = ips,
    formula = form,
    options = list(
      bru_verbose = 2,
      control.inla = list(int.strategy = "auto"),
      bru_initial = list(lsig = -1) # A difficult starting point
    )
  )

  #  ggplot(data.frame(distance = c(
  #    mexdolphin$points$distance,
  #    -mexdolphin$points$distance))) +
  #    geom_density(aes(distance, after_stat(count))) +
  #    geom_line(aes(distance, est),
  #              data = data.frame(distance = seq(-8,8, by = 0.01)) %>%
  #                mutate(est = hr(abs(distance),
  #                                fit$summary.fixed["lsig","mean"]) *
  #                         exp(fit$summary.fixed["Intercept","mean"])))
  #
  #  plot(ips$distance, hr(ips$distance, fit$summary.fixed["lsig", "mean"]))

  expect_equal(fit$summary.fixed["lsig", "mean"], 1.06, tolerance = midtol)
  expect_equal(fit$summary.fixed["lsig", "sd"], 0.5183252, tolerance = midtol)
  expect_equal(fit$summary.fixed["Intercept", "mean"], 2.29, tolerance = midtol)
  expect_equal(fit$summary.fixed["Intercept", "sd"], 0.2900139, tolerance = midtol)
})


# timings <- function() {
#   data(mexdolphin, package = "inlabru", envir = environment())
#
#   hr <- function(distance, lsig) {
#     1 - exp(-(distance / (exp(lsig)))^-1)
#   }
#   cmp <- ~ lsig(1) + Intercept(1)
#   form <- distance ~ log(hr(distance, lsig)) + Intercept
#   ips <- ipoints(INLA::inla.mesh.1d(seq(0, 8, by = 0.1)), name = "distance")
#   local_bru_options_set(bru_verbose = FALSE)
#   bench::mark(
#     pandemic = {
#       local_bru_options_set(bru_method = list(taylor = "pandemic"))
#       fit <- lgcp(
#         components = cmp,
#         mexdolphin$points,
#         ips = ips,
#         formula = form,
#         options = list(control.inla = list(int.strategy = "auto"))
#       )
#     },
#     check = FALSE,
#     min_iterations = 2
#   )
# }
#
# timings()
# # A tibble: 2 x 13
# expression      min   median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result memory time         gc
# <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list> <list> <list>       <list>
#   1 legacy        1.36s    1.39s     0.722        NA     2.17     2     6      2.77s <NULL> <NULL> <bch:tm [2]> <tibble [2 × 3]>
#   2 pandemic      1.36s    1.37s     0.732        NA     2.19     2     6      2.73s <NULL> <NULL> <bch:tm [2]> <tibble [2 × 3]>
#
# Unit: seconds
# expr      min       lq     mean   median       uq      max neval cld
# legacy 1.447667 1.474479 1.507082 1.496711 1.520239 1.899379   100   b
# pandemic 1.252985 1.271926 1.291995 1.288727 1.307383 1.402191   100  a
