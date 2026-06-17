library(dplyr)
library(ggplot2)

devtools::load_all()

timing <- NULL
nn <- c(1, 10, 100)
mm <- c(1, 3, 10,
        30, 100, 300, 1000,
        3000, 10000, 30000, 100000,
        300000, 1000000)
nm_max <- c(A = 3000, B = 1000000, C = 3000)
timer_max <- 10

n_df <- expand.grid(
  n = nn,
  m = mm
) |>
  arrange(desc(n * m)) |>
  dplyr::filter()

for (k in seq_len(nrow(n_df))) {
  n <- n_df[k, "n"]
  m <- n_df[k, "m"]
  N <- n * m
  if (N > max(nm_max)) {
    next
  }
  cat("n =", n, "m =", m, "N =", N, "\n")

  withr::local_seed(12345L)
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

  if (N <= nm_max["A"]) {
  A <- {
    local_bru_options_set(bru_method = list(
      autodiff = "pandemic",
      agg = "pandemic"
    ))
    model <- bru_model(
      ~ -1 + x(x, model = "iid"),
      bru_obs(
        formula = y ~ cos(x),
        is_rowwise = TRUE,
        data = data,
        response_data = response_data,
        family = "gaussian",
        aggregate = if (n == 1) NULL else "logsumexp"
      )
    )

    used <- bru_used(model$lhoods)
    expect_equal(used[["effect"]], "x")

    idx <- bru_index(model, used = used)
    inp <- bru_input(model, model$lhoods, inla_f = TRUE)

    iter <- 0
    toc <- proc.time()[1]
    while (toc[length(toc)] - toc[1] < timer_max) {
      iter <- iter + 1
      comp_lin <- ibm_as_taylor(model, input = inp, state = NULL)
      lin <- bru_compute_linearisation(
        model,
        lhoods = model$lhoods,
        input = inp,
        state = list(x = X),
        comp_simple = comp_lin
      )
      toc <- c(toc, proc.time()[1])
    }
    diff(toc)
  }
  A_iter <- iter
  if (length(A) > 1) {
    A_sd <- sd(A) / sqrt(iter)
  } else {
    A_sd <- NA
  }
  A <- mean(A)
  } else {
    A_iter <- NULL
    A_sd <- NULL
    A <- NULL
  }

  if (N <= nm_max["B"]) {
  B <- {
    local_bru_options_set(bru_method = list(
      autodiff = "fullchain",
      agg = "fullchain"
    ))
    model <- bru_model(
      ~ -1 + x(x, model = "iid"),
      bru_obs(
        formula = y ~ cos(x),
        is_rowwise = TRUE,
        data = data,
        response_data = response_data,
        family = "gaussian",
        aggregate = if (n == 1) NULL else "logsumexp"
      )
    )
    used <- bru_used(model$lhoods)
    expect_equal(used[["effect"]], "x")

    idx <- bru_index(model, used = used)
    inp <- bru_input(model, model$lhoods, inla_f = TRUE)

    iter <- 0
    toc <- proc.time()[1]
    while (toc[length(toc)] - toc[1] < timer_max) {
      iter <- iter + 1
      lin <- ibm_as_taylor(model,
                           input = inp,
                           state = list(x = X)
      )
      toc <- c(toc, proc.time()[1])
    }

    diff(toc)
  }
  B_iter <- iter
  if (length(B) > 1) {
    B_sd <- sd(B) / sqrt(iter)
  } else {
    B_sd <- NA
  }
  B <- mean(B)
  } else {
    B_iter <- NULL
    B_sd <- NULL
    B <- NULL
  }

  if ((N <= nm_max["C"]) || ((N <= nm_max["B"]) && (n == 1))) {
  C <- {
    local_bru_options_set(bru_method = list(
      autodiff = "fullchain",
      agg = "pandemic"
    ))
    model <- bru_model(
      ~ -1 + x(x, model = "iid"),
      bru_obs(
        formula = y ~ cos(x),
        is_rowwise = TRUE,
        data = data,
        response_data = response_data,
        family = "gaussian",
        aggregate = if (n == 1) NULL else "logsumexp"
      )
    )

    used <- bru_used(model$lhoods)
    expect_equal(used[["effect"]], "x")

    idx <- bru_index(model, used = used)
    inp <- bru_input(model, model$lhoods, inla_f = TRUE)

    iter <- 0
    toc <- proc.time()[1]
    while (toc[length(toc)] - toc[1] < timer_max) {
      iter <- iter + 1
      lin <- ibm_as_taylor(model,
                           input = inp,
                           state = list(x = X)
      )
      toc <- c(toc, proc.time()[1])
    }
    diff(toc)
  }
  C_iter <- iter
  if (length(C) > 1) {
    C_sd <- sd(C) / sqrt(iter)
  } else {
    C_sd <- NA
  }
  C <- mean(C)
  } else {
    C_iter <- NULL
    C_sd <- NULL
    C <- NULL
  }

  methods <- !c(is.null(A), is.null(B), is.null(C))

  timing_ <- tibble::tibble(
    n = n,
    m = m,
    N = n * m,
    Method = c("Monolith/Old", "Sequential/Chain", "Monolith/Chain")[methods],
    Iter = as.vector(c(A_iter, B_iter, C_iter)),
    Time = as.vector(c(A, B, C)),
    Time_sd = as.vector(c(A_sd, B_sd, C_sd))
  )

  timing <-
    bind_rows(
      timing, timing_
    )
}

timing <-
  timing |>
  group_by(n, m) |>
  mutate(
    Method_Reference = "Sequential/Chain",
    Reference = Time[Method == Method_Reference[1]]
  ) |>
  ungroup()

timing_est <- lapply(
  setNames(nm = unique(timing$Method)),
  function(r) {
    dat <- timing |> dplyr::filter(Method == r)
    nls(
      Time ~ (a1 * (n == nn[1]) + a2 * (n == nn[2]) + a3 * (n == nn[3])) +
        (b1 * (n == nn[1]) + b2 * (n == nn[2]) + b3 * (n == nn[3])) * N^d,
      data = dat,
      start = list(a1 = 0.01, a2 = 0.01, a3 = 0.01,
                   b1 = 0.01, b2 = 0.01, b3 = 0.01,
                   d = 1),
      weights = 1 / (dat$Time_sd)
    )
  }
)

timing <- timing |>
  group_by(Method) |>
  mutate(
    Time_Est = predict(timing_est[[Method[1]]], newdata = list(N = N, n = n))
  ) |>
  ungroup()

ggplot(
  timing |> mutate(Time_sd = ifelse(is.na(Time_sd), 0, Time_sd)),
  aes(
    x = N,
    y = Time,
    color = Method,
    fill = Method,
    shape = as.factor(n),
    linetype = as.factor(n)
  )
) +
  geom_ribbon(aes(ymin = pmax(0.01, Time - 2*Time_sd), ymax = Time + 2*Time_sd),
              alpha = 0.1) +
  geom_line() +
  geom_point() +
#  geom_line(aes(y = Time_Est), color = "black") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "N",
    y = "Seconds per evaluation",
    shape = "Aggregation",
    linetype = "Aggregation"
  ) +
  theme_minimal() +
  coord_cartesian(ylim = c(0.01, 5))
