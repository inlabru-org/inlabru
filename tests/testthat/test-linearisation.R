test_that("Linearisation", {
  skip_on_cran()
  local_bru_safe_inla()

  data <- data.frame(x = 1:10)
  data <- within(data, {
    y <- exp(x / 5) - 4 + rnorm(length(x), sd = 0.1)
    z <- rpois(length(x), exp(x / 5) + 4)
  })

  cmp <- ~ -1 + x + Int_y(1) + Int_z(1)
  # NOTE: Need explicit include specification to allow mapper construction for Int_y
  lhoods <-
    like_list(
      like(
        formula = y ~ exp(x) + Int_y_latent, data = data
      ),
      like(
        formula = z ~ exp(x) + Int_z, data = data, family = "poisson"
      )
    )
  lhoods <- bru_used_update(lhoods, names(component_list(cmp)))

  used <- bru_used(lhoods[[1]])
  expect_equal(used[["effect"]], "x")
  expect_equal(used[["latent"]], "Int_y")

  used <- bru_used(lhoods[[2]])
  expect_equal(used[["effect"]], c("x", "Int_z"))
  expect_equal(used[["latent"]], character(0))

  used <- bru_used(lhoods)
  expect_equal(used[["effect"]], c("x", "Int_z"))
  expect_equal(used[["latent"]], "Int_y")

  model <- bru_model(component_list(cmp, lhoods), lhoods)

  idx <- evaluate_index(model, lhoods)
  inp <- evaluate_inputs(model, lhoods, inla_f = FALSE)
  comp_lin <- evaluate_comp_lin(model, input = inp, state = NULL)
  lin0 <- bru_compute_linearisation.bru_model(
    model,
    lhoods = lhoods,
    input = inp,
    state = list(Int_y = 0, Int_z = 0, x = 0),
    comp_simple = comp_lin
  )
  lin <- bru_compute_linearisation.bru_model(
    model,
    lhoods = lhoods,
    input = inp,
    state = list(x = 1 / 5, Int_y = -4, Int_z = 4),
    comp_simple = comp_lin
  )

  stks0 <-
    lapply(
      seq_along(lhoods),
      function(lh_idx) {
        lh <- lhoods[[lh_idx]]
        lin_off <- ibm_eval(lin0[[lh_idx]], multi = TRUE, inla_f = TRUE)
        lin_A <- ibm_jacobian(lin0[[lh_idx]], multi = TRUE, inla_f = TRUE)
        nms <- names(lin_A)
        INLA::inla.stack(
          list(
            BRU.response = lh$response_data[[lh$response]],
            BRU.E = lh[["E"]],
            BRU.Ntrials = lh[["Ntrials"]],
            BRU.weights = lh[["weights"]],
            BRU.offset = as.vector(lin_off)
          ),
          A = lapply(nms, function(nm) {
            lin_A[[nm]][, idx[["inla_subset"]][[nm]], drop = FALSE]
          }),
          effects = idx[["idx_inla"]][nms]
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
