#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom rlang .env

#' @title Plot inlabru convergence diagnostics
#'
#' @description
#' Draws four panels of convergence diagnostics for an iterated INLA method
#' estimation
#'
#' @param x a [bru] object, typically a result from [bru()] for a nonlinear
#' predictor model
#' @param from,to integer values for the range of iterations to plot.
#' Default `from = 1` (start from the first iteration) and `to = NULL` (end at
#' the last iteration).
#' Set `from = 0` to include the initial linearisation point in the track plot.
#' @param type `r lifecycle::badge("experimental")` character; "bru" (default)
#'   for iterative nonlinear inlabru convergence diagnostics plots, or "inla"
#'   for INLA optimiser trace plots.
#' @return A ggplot object with four panels of convergence diagnostics:
#' - `Tracks`: Mode and linearisation values for each effect
#' - `Mode - Lin`: Difference between mode and linearisation values for each
#'   effect
#' - `|Change| / sd`: Absolute change in mode and linearisation values
#' divided by the standard deviation for each effect
#' - `Change & sd`: Absolute change in mode and linearisation values
#' and standard deviation for each effect
#'
#' For multidimensional components, only the overall average, maximum, and
#' minimum values are shown.
#' @seealso [bru()]
#'
#' @details Requires the "dplyr", "ggplot2", "magrittr", and "patchwork"
#' packages to be installed.
#' @export
#' @examples
#' \dontrun{
#' fit <- bru(...)
#' bru_convergence_plot(fit)
#' }
bru_convergence_plot <- function(x, from = 1, to = NULL, type = NULL) {
  stopifnot(inherits(x, "bru"))
  x <- bru_check_object_bru(x)

  type <- match.arg(type, c("bru", "inla"))
  if (identical(type, "bru")) {
    make_bru_track_plots(x, from = from, to = to)[["default"]]
  } else if (identical(type, "inla")) {
    make_inla_track_plots(x, from = from, to = to)[["default"]]
  }
}


make_bru_track_plots <- function(fit, from = 1, to = NULL) {
  needed <- c("dplyr", "ggplot2", "magrittr", "patchwork")
  are_installed <-
    vapply(
      needed,
      function(x) {
        requireNamespace(x, quietly = TRUE)
      },
      TRUE
    )
  if (any(!are_installed)) {
    stop(
      paste0(
        "Needed package(s) ",
        paste0("'", needed[!are_installed], "'", collapse = ", "),
        " not installed, but are needed by make_track_plots()"
      )
    )
  }

  track_data <- fit$bru_iinla$track
  if (!is.null(from)) {
    stopifnot(is.numeric(from))
    stopifnot(from >= 0)
  } else {
    from <- min(track_data$iteration)
  }
  if (!is.null(to)) {
    stopifnot(is.numeric(to))
    stopifnot(to >= 0)
  } else {
    to <- max(track_data$iteration)
  }

  track_data <- track_data %>%
    dplyr::group_by(.data$effect, .data$iteration) %>%
    dplyr::mutate(alpha_value = 1 / sqrt(max(.data$index))) %>%
    dplyr::ungroup()

  alpha_value_lin_scale <- 0.8

  lty_ <- factor(c("Mode", "Lin", "Mode-Lin", "SD"),
    levels = c("Mode", "Lin", "Mode-Lin", "SD")
  )
  names(lty_) <- levels(lty_)
  col_ <- factor(c("Max", "Mean", "Min"),
    levels = c("Max", "Mean", "Min")
  )
  names(col_) <- levels(col_)

  sc <- ggplot2::scale_linetype_discrete(
    name = "Quantity",
    breaks = names(lty_),
    labels = labels(lty_),
    drop = FALSE,
    na.translate = FALSE
  )
  # Colour blind friendliness, see https://davidmathlogic.com/colorblind/
  sc_minmax <- list(
    ggplot2::scale_color_discrete(
      name = "Aspect",
      type = c("#DC3220", "#000000", "#005AB5"),
      breaks = names(col_),
      labels = labels(col_),
      drop = FALSE,
      na.translate = FALSE
    ),
    ggplot2::scale_fill_discrete(
      name = "Aspect",
      type = c("#DC3220", "#000000", "#005AB5"),
      breaks = names(col_),
      labels = labels(col_),
      drop = FALSE,
      na.translate = FALSE
    )
  )

  pl_theme_abs <-
    list(
      ggplot2::facet_wrap(ggplot2::vars(.data$effect), scales = "free_y"),
      sc,
      sc_minmax #
      #      ggplot2::guides(color = "none")
    )
  pl_theme_rel <-
    list(
      ggplot2::facet_wrap(ggplot2::vars(.data$effect), scales = "free_y"),
      sc,
      sc_minmax # ,
      #      ggplot2::guides(color = "none")
    )
  pl_theme_norm <-
    list(
      ggplot2::facet_wrap(ggplot2::vars(.data$effect), scales = "free_y"),
      sc,
      sc_minmax # ,
      #      ggplot2::guides(color = "none")
    )

  pl_tracks <-
    ggplot2::ggplot(
      track_data %>%
        dplyr::filter(.data$iteration >= from, .data$iteration <= to) %>%
        dplyr::group_by(
          .data$effect,
          .data$iteration
        ) %>%
        dplyr::summarise(
          MaxMode = max(.data$mode),
          MeanMode = mean(.data$mode),
          MinMode = min(.data$mode),
          MaxLin = max(.data$new_linearisation),
          MeanLin = mean(.data$new_linearisation),
          MinLin = min(.data$new_linearisation),
          .groups = "drop"
        )
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        .data$MaxMode,
        lty = lty_["Mode"],
        color = col_["Max"]
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        .data$MinMode,
        lty = lty_["Mode"],
        color = col_["Min"]
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        .data$MeanMode,
        lty = lty_["Mode"],
        color = col_["Mean"]
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        .data$MaxLin,
        lty = lty_["Lin"],
        col = "Max"
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        .data$MinLin,
        lty = lty_["Lin"],
        col = "Min"
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        .data$MeanLin,
        lty = lty_["Lin"],
        col = "Mean"
      ),
      na.rm = TRUE
    ) +
    ggplot2::ylab("Value") +
    pl_theme_abs +
    ggplot2::ggtitle("Tracks")

  pl_mode_lin <-
    ggplot2::ggplot(
      track_data %>%
        dplyr::filter(
          .data$iteration >= max(1, from),
          .data$iteration <= to
        ) %>%
        dplyr::filter(
          is.finite(.data$sd),
          is.finite(.data$new_linearisation)
        ) %>%
        dplyr::group_by(
          .data$effect,
          .data$iteration
        ) %>%
        dplyr::summarise(
          MaxSD = max(.data$sd),
          MeanSD = mean(.data$sd),
          MinSD = min(.data$sd),
          MaxDiff = max(.data$mode - .data$new_linearisation),
          MeanDiff = mean(.data$mode - .data$new_linearisation),
          MinDiff = min(.data$mode - .data$new_linearisation),
          .groups = "drop"
        )
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        .data$iteration,
        ymin = -.data$MaxSD,
        ymax = .data$MaxSD,
        fill = col_["Max"],
        lty = lty_["SD"]
      ),
      alpha = 0.1,
      na.rm = TRUE
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        .data$iteration,
        ymin = -.data$MeanSD,
        ymax = .data$MeanSD,
        fill = col_["Mean"],
        lty = lty_["SD"]
      ),
      alpha = 0.1,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        .data$MaxDiff,
        col = col_["Max"],
        lty = lty_["Mode-Lin"]
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        .data$MinDiff,
        col = col_["Min"],
        lty = lty_["Mode-Lin"]
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        .data$MeanDiff,
        col = col_["Mean"],
        lty = lty_["Mode-Lin"]
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        .data$MaxSD,
        col = col_["Max"],
        lty = lty_["SD"]
      ),
      alpha = 0.5,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        .data$MinSD,
        col = col_["Min"],
        lty = lty_["SD"]
      ),
      alpha = 0.5
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        .data$MeanSD,
        col = col_["Mean"],
        lty = lty_["SD"]
      ),
      alpha = 0.5,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        -.data$MaxSD,
        col = col_["Max"],
        lty = lty_["SD"]
      ),
      alpha = 0.5,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        -.data$MinSD,
        col = col_["Min"],
        lty = lty_["SD"]
      ),
      alpha = 0.5,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        -.data$MeanSD,
        col = col_["Mean"],
        lty = lty_["SD"]
      ),
      alpha = 0.5,
      na.rm = TRUE
    ) +
    ggplot2::guides(alpha = "none", fill = "none") +
    pl_theme_rel +
    ggplot2::ggtitle("Mode - Lin")

  track_data_prev <- dplyr::filter(
    track_data,
    .data$iteration < max(.data$iteration)
  )
  track_data_prev$iteration <- track_data_prev$iteration + 1
  track_data <- dplyr::left_join(track_data, track_data_prev,
    by = c("effect", "iteration", "index"),
    suffix = c("", ".prev")
  )

  pl_relative_change <-
    ggplot2::ggplot(
      track_data %>%
        dplyr::filter(
          .data$iteration >= max(1, from),
          .data$iteration <= to
        ) %>%
        dplyr::group_by(
          .data$effect,
          .data$iteration
        ) %>%
        dplyr::mutate(sd = dplyr::if_else(is.finite(.data$sd),
          .data$sd, 1.0
        )) %>%
        dplyr::summarise(
          MaxMode = max(abs(.data$mode - .data$mode.prev) / .data$sd),
          MeanMode = mean(abs(.data$mode - .data$mode.prev) / .data$sd),
          RMSMode = mean((.data$mode - .data$mode.prev)^2 / .data$sd^2)^0.5,
          MaxLin = max(abs(.data$new_linearisation -
            .data$new_linearisation.prev) / .data$sd),
          MeanLin = mean(abs(.data$new_linearisation -
            .data$new_linearisation.prev) / .data$sd),
          RMSLin = mean((.data$new_linearisation -
            .data$new_linearisation.prev)^2 / .data$sd^2)^0.5,
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          MaxMode = dplyr::if_else(.data$MaxMode > 0, .data$MaxMode, NA),
          MeanMode = dplyr::if_else(.data$MeanMode > 0, .data$MeanMode, NA),
          RMSMode = dplyr::if_else(.data$RMSMode > 0, .data$RMSMode, NA),
          MaxLin = dplyr::if_else(.data$MaxLin > 0, .data$MaxLin, NA),
          MeanLin = dplyr::if_else(.data$MeanLin > 0, .data$MeanLin, NA),
          RMSLin = dplyr::if_else(.data$RMSLin > 0, .data$RMSLin, NA),
        )
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        .data$MaxMode,
        lty = lty_["Mode"],
        col = col_["Max"]
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        .data$MeanMode,
        lty = lty_["Mode"],
        col = col_["Mean"]
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        .data$MaxLin,
        lty = lty_["Lin"],
        col = col_["Max"]
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        .data$MeanLin,
        lty = lty_["Lin"],
        col = col_["Mean"]
      ),
      na.rm = TRUE
    ) +
    pl_theme_norm +
    ggplot2::scale_y_log10() +
    ggplot2::ggtitle("|Change| / sd") +
    ggplot2::geom_hline(yintercept = fit$bru_info$options$bru_method$rel_tol)

  pl_change <-
    ggplot2::ggplot(
      track_data %>%
        dplyr::filter(
          .data$iteration >= max(1, from),
          .data$iteration <= to
        ) %>%
        dplyr::filter(is.finite(.data$sd)) %>%
        dplyr::group_by(
          .data$effect,
          .data$iteration
        ) %>%
        dplyr::summarise(
          MaxSD = max(.data$sd),
          MeanSD = mean(.data$sd),
          MinSD = min(.data$sd),
          MaxMode = max(.data$mode - .data$mode.prev),
          MeanMode = mean(.data$mode - .data$mode.prev),
          MinMode = min(.data$mode - .data$mode.prev),
          MaxLin = max(.data$new_linearisation -
            .data$new_linearisation.prev),
          MeanLin = mean(.data$new_linearisation -
            .data$new_linearisation.prev),
          MinLin = min(.data$new_linearisation -
            .data$new_linearisation.prev),
          .groups = "drop"
        )
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        .data$iteration,
        ymin = -.data$MaxSD,
        ymax = .data$MaxSD,
        fill = col_["Max"],
        lty = lty_["SD"]
      ),
      alpha = 0.1,
      na.rm = TRUE
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        .data$iteration,
        ymin = -.data$MeanSD,
        ymax = .data$MeanSD,
        fill = col_["Mean"],
        lty = lty_["SD"]
      ),
      alpha = 0.1,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        .data$MaxMode,
        col = col_["Max"],
        lty = lty_["Mode"]
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        .data$MinMode,
        col = col_["Min"],
        lty = lty_["Mode"]
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        .data$MeanMode,
        col = col_["Mean"],
        lty = lty_["Mode"]
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        .data$MaxLin,
        col = col_["Max"],
        lty = lty_["Lin"]
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        .data$MinLin,
        col = col_["Min"],
        lty = lty_["Lin"]
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        .data$MeanLin,
        col = col_["Mean"],
        lty = lty_["Lin"]
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        .data$MaxSD,
        col = col_["Max"],
        lty = lty_["SD"]
      ),
      alpha = 0.5,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        .data$MinSD,
        col = col_["Min"],
        lty = lty_["SD"]
      ),
      alpha = 0.5,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        .data$MeanSD,
        col = col_["Mean"],
        lty = lty_["SD"]
      ),
      alpha = 0.5,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        -.data$MaxSD,
        col = col_["Max"],
        lty = lty_["SD"]
      ),
      alpha = 0.5,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        -.data$MinSD,
        col = col_["Min"],
        lty = lty_["SD"]
      ),
      alpha = 0.5,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        -.data$MeanSD,
        col = col_["Mean"],
        lty = lty_["SD"]
      ),
      alpha = 0.5,
      na.rm = TRUE
    ) +
    ggplot2::guides(alpha = "none", fill = "none") +
    pl_theme_rel +
    ggplot2::ggtitle("Change & sd")

  pl_combined <-
    ((
      pl_tracks +
        ggplot2::geom_line(
          ggplot2::aes(.data$iteration, NA_real_, lty = lty_["Mode-Lin"]),
          na.rm = TRUE
        ) +
        ggplot2::geom_line(
          ggplot2::aes(.data$iteration, NA_real_, lty = lty_["SD"]),
          na.rm = TRUE
        ) |
        pl_mode_lin + ggplot2::guides(linetype = "none", color = "none")
    ) /
      (
        pl_relative_change +
          ggplot2::guides(linetype = "none", color = "none") |
          pl_change + ggplot2::guides(linetype = "none", color = "none")
      )
    ) +
      patchwork::plot_layout(guides = "collect") &
      ggplot2::theme(legend.position = "bottom")


  list(
    tracks = pl_tracks,
    mode_lin = pl_mode_lin,
    relative_change = pl_relative_change,
    change = pl_change,
    default = pl_combined
  )
}

make_inla_track_plots <- function(fit, from = 1, to = NULL) {
  needed <- c("dplyr", "ggplot2", "magrittr", "patchwork")
  are_installed <-
    vapply(
      needed,
      function(x) {
        requireNamespace(x, quietly = TRUE)
      },
      TRUE
    )
  if (any(!are_installed)) {
    stop(
      paste0(
        "Needed package(s) ",
        paste0("'", needed[!are_installed], "'", collapse = ", "),
        " not installed, but are needed by make_track_plots()"
      )
    )
  }

  track_data <- fit$bru_iinla$inla_track
  if (!is.null(from)) {
    stopifnot(is.numeric(from))
    stopifnot(from >= 0)
  } else {
    from <- min(track_data$iteration)
  }
  if (!is.null(to)) {
    stopifnot(is.numeric(to))
    stopifnot(to >= 0)
  } else {
    to <- max(track_data$iteration)
  }

  colnames(track_data$theta) <- paste0("theta", seq_len(ncol(track_data$theta)))
  rownames(track_data$theta) <- NULL

  alpha_value_lin_scale <- 0.8

  pl_trace1a <-
    ggplot2::ggplot(
      track_data %>%
        dplyr::filter(.data$iteration >= from, .data$iteration <= to),
      ggplot2::aes(
        .data$nfunc_total,
        .data$f,
        col = as.factor(.data$iteration)
      )
    ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::ggtitle("target function value")

  tracks <- cbind(
    track_data %>%
      dplyr::select(.data$iteration, .data$f, .data$nfunc, .data$nfunc_total),
    tibble::as_tibble(track_data$theta)
  )
  tracks <- tracks %>%
    tidyr::pivot_longer(
      cols = -c(.data$iteration, .data$f, .data$nfunc, .data$nfunc_total),
      names_to = "theta",
      values_to = "value"
    )
  pl_trace1b <-
    ggplot2::ggplot(
      tracks %>% dplyr::filter(.data$iteration >= from, .data$iteration <= to),
      ggplot2::aes(
        .data$nfunc_total,
        .data$value,
        col = .data$theta,
        lty = .data$theta
      )
    ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::ggtitle("theta values")

  pl_trace1b <-
    ggplot2::ggplot(
      tracks,
      ggplot2::aes(
        .data$nfunc_total,
        .data$value,
        col = .data$theta,
        lty = .data$theta
      )
    ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::ggtitle("theta values")

  pl_trace_combined <-
    (
      pl_trace1a | pl_trace1b
      #        ggplot2::guides(linetype = "none", color = "none")
    ) +
      patchwork::plot_layout(guides = "collect") &
      ggplot2::theme(legend.position = "bottom")


  tr <- track_data %>%
    dplyr::select(.data$iteration, .data$f, .data$nfunc, .data$nfunc_total)
  th <- tibble::as_tibble(track_data$theta)
  tr1 <- cbind(tr[-nrow(tr), ], theta = th[-nrow(tr), ])
  tr2 <- cbind(tr[-nrow(tr), ], theta = th[-1, ])

  track1 <- tr1 %>% tidyr::pivot_longer(
    cols = -c(.data$iteration, .data$f, .data$nfunc, .data$nfunc_total),
    names_to = "theta",
    values_to = "value"
  )
  track2 <- tr2 %>% tidyr::pivot_longer(
    cols = -c(.data$iteration, .data$f, .data$nfunc, .data$nfunc_total),
    names_to = "theta",
    values_to = "value"
  )
  track1$value_next <- track2$value

  pl_trace2 <-
    ggplot2::ggplot(
      track1 %>% dplyr::filter(.data$iteration >= from, .data$iteration <= to)
    ) +
    ggplot2::geom_path(ggplot2::aes(.data$value, .data$value_next)) +
    ggplot2::facet_wrap(~ .data$theta) +
    ggplot2::ggtitle("theta")

  list(
    trace1a = pl_trace1a,
    trace1b = pl_trace1b,
    trace2 = pl_trace2,
    default = pl_trace_combined
  )
}



#' @title Plot inlabru iteration timings
#'
#' @description
#' Draws the time per iteration for preprocessing (including linearisation),
#' `inla()` calls, and
#' line search. Iteration `0` is the time used for defining the model structure.
#'
#' @param x a [bru] object, typically a result from [bru()] for a nonlinear
#' predictor model
#'
#' @details Requires the "ggplot2" package to be installed.
#' @export
#' @examples
#' \dontrun{
#' fit <- bru(...)
#' bru_timings_plot(fit)
#' }
bru_timings_plot <- function(x) {
  needed <- c("ggplot2")
  are_installed <-
    vapply(
      needed,
      function(x) {
        requireNamespace(x, quietly = TRUE)
      },
      TRUE
    )
  if (any(!are_installed)) {
    stop(
      paste0(
        "Needed package(s) ",
        paste0("'", needed[!are_installed], "'", collapse = ", "),
        " not installed, but are needed by bru_timings_plot()"
      )
    )
  }

  timings <- bru_timings(x)

  ggplot2::ggplot(
    tidyr::pivot_longer(
      dplyr::rename(timings, CPU = .data$Time),
      cols = c("CPU", "System", "Elapsed"),
      names_to = "Time",
      values_to = "Value"
    )
  ) +
    ggplot2::geom_point(ggplot2::aes(
      .data$Iteration,
      .data$Value,
      col = .data$Task,
      shape = .data$Task
    )) +
    ggplot2::geom_line(ggplot2::aes(
      .data$Iteration,
      .data$Value,
      col = .data$Task,
      linetype = .data$Time
    )) +
    ggplot2::scale_y_continuous() +
    ggplot2::ylab("Time (s)") +
    ggplot2::guides(
      color = ggplot2::guide_legend(title = "Task"),
      linetype = ggplot2::guide_legend(title = "Time")
    )
}
