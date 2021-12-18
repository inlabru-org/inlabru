local_bru_testthat_setup()

test_that("Linearisation", {
  skip_on_cran()
  local_bru_safe_inla()

  data <- data.frame(x = 1:10)
  data <- within(data, {
    y <- exp(x / 5) - 4 + rnorm(length(x), sd = 0.1)
    z <- rpois(length(x), exp(x / 5) + 4)
  })

  cmp <- ~ -1 + x + Int_y(1) + Int_z(1)
  lhoods <-
    like_list(
      like(
        formula = y ~ exp(x) + Int_y_latent, data = data,
        include = c("x", "Int_y")
      ),
      like(
        formula = z ~ exp(x) + Int_z, data = data, family = "poisson",
        include = c("x", "Int_z")
      )
    )
  model <- bru_model(component_list(cmp, lhoods), lhoods)


  idx <- evaluate_index(model, lhoods)
  A <- evaluate_A(model, lhoods, inla_f = FALSE)
  lin0 <- bru_compute_linearisation.bru_model(model, lhoods,
    state = list(Int_y = 0, Int_z = 0, x = 0), A = A
  )
  lin <- bru_compute_linearisation.bru_model(model, lhoods,
    state = list(x = 1 / 5, Int_y = -4, Int_z = 4), A = A
  )

  stks0 <-
    lapply(
      seq_along(lhoods),
      function(lh_idx) {
        lh <- lhoods[[lh_idx]]
        lin_A <- lin0[[lh_idx]]$A
        INLA::inla.stack(
          list(
            BRU.response = lh$data[[lh$response]],
            BRU.E = lh[["E"]],
            BRU.Ntrials = lh[["Ntrials"]],
            BRU.offset = as.vector(lin0[[lh_idx]]$offset)
          ),
          A = lapply(names(lin_A), function(nm) {
            lin_A[[nm]][, idx[["inla_subset"]][[nm]], drop = FALSE]
          }),
          effects = idx[["idx_inla"]][names(lin0[[lh_idx]]$A)]
        )
      }
    )

  stk0 <-
    do.call(
      inlabru::inla.stack.mjoin,
      c(stks0, list(compress = TRUE, remove.unused = FALSE))
    )

  stk0_ <- bru_make_stack.bru_like_list(lhoods, lin0, idx)

  expect_s3_class(stk0, "inla.data.stack")

  expect_error(
    object = {
      fit <- bru(
        components = cmp,
        lhoods,
        options = list(
          control.inla = list(int.strategy = "eb"),
          bru_verbose = FALSE,
          bru_method = list(
            taylor = "pandemic",
            search = "all"
          )
        )
      )
    },
    NA
  )
})
