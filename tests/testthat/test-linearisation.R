test_that("Linearisation", {
  skip_on_cran()
  local_bru_safe_inla()

  withr::local_seed(12345L)
  data <- data.frame(x = seq_len(10) / 1)
  data <- within(data, {
    y <- exp(x / 5) - 2 + rnorm(length(x), sd = 0.1)
    z <- rpois(length(x), (exp(x / 5) + exp(2)))
  })

  cmp <- ~ -1 + x + Int_y(1) + Int_z(1)
  lhoods <-
    c(
      bru_obs(
        formula = y ~ exp(x) + Int_y_latent,
        data = data
      ),
      bru_obs(
        formula = z ~ log(exp(x) + exp(Int_z_latent)),
        data = data,
        family = "poisson"
      )
    )
  lhoods <- bru_used_update(lhoods, names(bru_comp_list(cmp)))

  used <- bru_used(lhoods[[1]])
  expect_equal(used[["effect"]], "x")
  expect_equal(used[["latent"]], "Int_y")

  used <- bru_used(lhoods[[2]])
  expect_equal(used[["effect"]], "x")
  expect_equal(used[["latent"]], "Int_z")

  used <- bru_used(lhoods)
  expect_equal(used[["effect"]], "x")
  expect_equal(used[["latent"]], c("Int_y", "Int_z"))

  model <- bru_model(bru_comp_list(cmp), lhoods)
  lhoods <- model$lhoods

  idx <- bru_index(model, used = bru_used(lhoods))
  inp <- bru_input(model, lhoods, inla_f = TRUE)

  if (identical(bru_options_get("bru_method")$autodiff, "fullchain")) {
    lin0 <- ibm_as_taylor(model,
      input = inp,
      state = list(Int_y = 0, Int_z = 0, x = 0)
    )
    lin <- ibm_as_taylor(model,
      input = inp,
      state = list(x = 1 / 5, Int_y = -4, Int_z = log(4))
    )
  } else {
    comp_lin <- ibm_as_taylor(model, input = inp, state = NULL)
    lin0 <- bru_compute_linearisation(
      model,
      lhoods = lhoods,
      input = inp,
      state = list(Int_y = 0, Int_z = 0, x = 0),
      comp_simple = comp_lin
    )
    lin <- bru_compute_linearisation(
      model,
      lhoods = lhoods,
      input = inp,
      state = list(x = 1 / 5, Int_y = -4, Int_z = log(4)),
      comp_simple = comp_lin
    )
  }

  if (utils::packageVersion("INLA") > "24.06.02") {
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
              BRU.E = lh[["response_data"]][["BRU_E"]],
              BRU.Ntrials = lh[["response_data"]][["BRU_Ntrials"]],
              BRU.weights = lh[["response_data"]][["BRU_weights"]],
              BRU.scale = lh[["response_data"]][["BRU_scale"]],
              BRU.offset = as.vector(lin_off)
            ),
            A = lapply(nms, function(nm) {
              lin_A[[nm]][, idx[["inla_subset"]][[nm]], drop = FALSE]
            }),
            effects = idx[["idx_inla"]][nms],
            responses = list(lh$response_data[[lh$response]])
          )
        }
      )
  } else {
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
              BRU.E = lh[["response_data"]][["BRU_E"]],
              BRU.Ntrials = lh[["response_data"]][["BRU_Ntrials"]],
              BRU.weights = lh[["response_data"]][["BRU_weights"]],
              BRU.scale = lh[["response_data"]][["BRU_scale"]],
              BRU.offset = as.vector(lin_off)
            ),
            A = lapply(nms, function(nm) {
              lin_A[[nm]][, idx[["inla_subset"]][[nm]], drop = FALSE]
            }),
            effects = idx[["idx_inla"]][nms]
          )
        }
      )
  }

  stk0 <-
    do.call(
      bru_inla.stack.mjoin,
      c(stks0, list(compress = TRUE, remove.unused = FALSE))
    )

  stk0_ <- bru_make_stack(lhoods, lin0, idx)

  expect_s3_class(stk0, "inla.data.stack")

  expect_no_error(
    object = {
      fit <- bru(
        components = cmp,
        lhoods,
        options = list(
          #          bru_initial = list(
          #            x = 1 / 5, Int_y = -4, Int_z = log(4)
          #          ),
          control.inla = list(int.strategy = "eb"),
          bru_verbose = FALSE,
          bru_method = list(
            taylor = "pandemic",
            search = "all"
          )
        )
      )
    }
  )
})


test_that("Linearisation 2", {
  skip_on_cran()
  local_bru_safe_inla()

  withr::local_seed(12345L)
  n <- 3
  m <- 5
  N <- n * m
  data <- data.frame(
    x = seq_len(N),
    .block = rep(seq_len(N / n), each = n)
  )
  data <- within(data, {
    eta <- cos(x)
  })
  response_data <- data.frame(
    y = fm_block_logsumexp_eval(block = data$.block, values = data$eta)
  )
  X <- rnorm(N)

  #  profvis::profvis({
  #  bench::mark(
  A <- {
    local_bru_options_set(bru_method = list(
      autodiff = "pandemic",
      agg = "pandemic"
    ))
    model <- bru_model(
      components = bru_comp_list(~ -1 + x(x, model = "iid")),
      lhoods = bru_obs(
        formula = y ~ cos(x),
        is_rowwise = TRUE,
        data = data,
        response_data = response_data,
        family = "gaussian",
        aggregate = "logsumexp"
      )
    )

    used <- bru_used(model$lhoods)
    expect_equal(used[["effect"]], "x")

    idx <- bru_index(model, used = used)
    inp <- bru_input(model, model$lhoods, inla_f = TRUE)

    comp_lin <- ibm_as_taylor(model, input = inp, state = NULL)
    lin0 <- bru_compute_linearisation(
      model,
      lhoods = model$lhoods,
      input = inp,
      state = list(x = X),
      comp_simple = comp_lin
    )
    lin <- bru_compute_linearisation(
      model,
      lhoods = model$lhoods,
      input = inp,
      state = list(x = X / 5),
      comp_simple = comp_lin
    )
    list(lin0, lin)
  }
  #   ,
  B <- {
    local_bru_options_set(bru_method = list(
      autodiff = "fullchain",
      agg = "fullchain"
    ))
    model <-
      bru_model(
        components = ~ -1 + x(x, model = "iid"),
        bru_obs(
          formula = y ~ cos(x),
          is_rowwise = TRUE,
          data = data,
          response_data = response_data,
          family = "gaussian",
          aggregate = "logsumexp"
        )
      )

    used <- bru_used(model$lhoods)
    expect_equal(used[["effect"]], "x")

    idx <- bru_index(model, used = used)
    inp <- bru_input(model, model$lhoods, inla_f = TRUE)

    lin0 <- ibm_as_taylor(model,
      input = inp,
      state = list(x = X)
    )
    lin <- ibm_as_taylor(model,
      input = inp,
      state = list(x = X / 5)
    )
    list(lin0, lin)
  }
  #    ,
  C <- {
    local_bru_options_set(bru_method = list(
      autodiff = "fullchain",
      agg = "pandemic"
    ))
    model <- bru_model(
      components = ~ -1 + x(x, model = "iid"),
      bru_obs(
        formula = y ~ cos(x),
        #              is_rowwise = TRUE,
        data = data,
        response_data = response_data,
        family = "gaussian",
        aggregate = "logsumexp"
      )
    )

    used <- bru_used(model$lhoods)
    expect_equal(used[["effect"]], "x")

    idx <- bru_index(model, used = used)
    inp <- bru_input(model, model$lhoods, inla_f = TRUE)

    lin0 <- ibm_as_taylor(model,
      input = inp,
      state = list(x = X)
    )
    lin <- ibm_as_taylor(model,
      input = inp,
      state = list(x = X / 5)
    )
    list(lin0, lin)
  }
  #    ,
  #    check = FALSE
  #    )
  # })

  True <- -sin(X[1:n]) * exp(cos(X[1:n])) / sum(exp(cos(X[1:n])))
  expect_lt(max(abs(A[[1]][[1]]$jacobian$x[1, 1:n] - True)), 1e-5)
  expect_lt(max(abs(B[[1]][[1]]$jacobian$x[1, 1:n] - True)), 5e-7)
  expect_lt(max(abs(C[[1]][[1]]$jacobian$x[1, 1:n] - True)), 1e-5)
})
