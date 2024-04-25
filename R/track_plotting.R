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
#'
#' @details Requires the "dplyr", "ggplot2", "magrittr", and "patchwork"
#' packages to be installed.
#' @export
#' @examples
#' \dontrun{
#' fit <- bru(...)
#' bru_convergence_plot(fit)
#' }
bru_convergence_plot <- function(x) {
  stopifnot(inherits(x, "bru"))
  x <- bru_check_object_bru(x)
  make_track_plots(x)[["default"]]
}


make_track_plots <- function(fit) {
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
  track_data <-
    fit$bru_iinla$track %>%
    dplyr::left_join(
      fit$bru_iinla$track %>%
        dplyr::filter(.data$iteration == max(.data$iteration)) %>%
        dplyr::rename(
          mode_end = .data$mode,
          new_linearisation_end = .data$new_linearisation,
          sd_end = .data$sd
        ),
      by = c("effect", "index"),
    ) %>%
    dplyr::mutate(iteration = .data$iteration.x) %>%
    dplyr::left_join(
      fit$bru_iinla$track %>%
        dplyr::filter(.data$iteration == 1) %>%
        dplyr::rename(
          mode_start = .data$mode,
          new_linearisation_start = .data$new_linearisation,
          sd_start = .data$sd
        ),
      by = c("effect", "index"),
    ) %>%
    dplyr::mutate(iteration = .data$iteration.x)

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
  sc_minmax <- list(ggplot2::scale_color_discrete(
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
  ))

  pl_theme_abs <-
    list(
      ggplot2::facet_wrap(ggplot2::vars(.data$effect), scales = "free"),
      sc,
      sc_minmax #
      #      ggplot2::guides(color = "none")
    )
  pl_theme_rel <-
    list(
      ggplot2::facet_wrap(ggplot2::vars(.data$effect), scales = "free"),
      sc,
      sc_minmax # ,
      #      ggplot2::guides(color = "none")
    )
  pl_theme_norm <-
    list(
      ggplot2::facet_wrap(ggplot2::vars(.data$effect)),
      sc,
      sc_minmax # ,
      #      ggplot2::guides(color = "none")
    )

  pl_tracks <-
    ggplot2::ggplot(
      track_data %>%
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
        ymin = -MaxSD,
        ymax = MaxSD,
        fill = col_["Max"],
        lty = lty_["SD"]
      ),
      alpha = 0.1,
      na.rm = TRUE
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        .data$iteration,
        ymin = -MeanSD,
        ymax = MeanSD,
        fill = col_["Mean"],
        lty = lty_["SD"]
      ),
      alpha = 0.1,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        MaxDiff,
        col = col_["Max"],
        lty = lty_["Mode-Lin"]
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        MinDiff,
        col = col_["Min"],
        lty = lty_["Mode-Lin"]
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        MeanDiff,
        col = col_["Mean"],
        lty = lty_["Mode-Lin"]
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        MaxSD,
        col = col_["Max"],
        lty = lty_["SD"]
      ),
      alpha = 0.5,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        MinSD,
        col = col_["Min"],
        lty = lty_["SD"]
      ),
      alpha = 0.5
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        MeanSD,
        col = col_["Mean"],
        lty = lty_["SD"]
      ),
      alpha = 0.5,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        -MaxSD,
        col = col_["Max"],
        lty = lty_["SD"]
      ),
      alpha = 0.5,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        -MinSD,
        col = col_["Min"],
        lty = lty_["SD"]
      ),
      alpha = 0.5,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        -MeanSD,
        col = col_["Mean"],
        lty = lty_["SD"]
      ),
      alpha = 0.5,
      na.rm = TRUE
    ) +
    ggplot2::guides(alpha = "none", fill = "none") +
    pl_theme_rel +
    ggplot2::ggtitle("Mode - Lin")

  track_data_prev <- dplyr::filter(track_data, .data$iteration < max(.data$iteration))
  track_data_prev$iteration <- track_data_prev$iteration + 1
  track_data <- dplyr::left_join(track_data, track_data_prev,
    by = c("effect", "iteration", "index"),
    suffix = c("", ".prev")
  )

  pl_relative_change <-
    ggplot2::ggplot(
      track_data %>%
        dplyr::group_by(
          .data$effect,
          .data$iteration
        ) %>%
        dplyr::mutate(sd = dplyr::if_else(is.finite(.data$sd), .data$sd, 1.0)) %>%
        dplyr::summarise(
          MaxMode = max(abs(.data$mode - .data$mode.prev) / .data$sd),
          MeanMode = mean(abs(.data$mode - .data$mode.prev) / .data$sd),
          RMSMode = mean((.data$mode - .data$mode.prev)^2 / .data$sd^2)^0.5,
          MaxLin = max(abs(.data$new_linearisation - .data$new_linearisation.prev) / .data$sd),
          MeanLin = mean(abs(.data$new_linearisation - .data$new_linearisation.prev) / .data$sd),
          RMSLin = mean((.data$new_linearisation - .data$new_linearisation.prev)^2 / .data$sd^2)^0.5,
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
    ggplot2::ggtitle("|Change| / sd (Max and Mean)") +
    ggplot2::geom_hline(yintercept = fit$bru_info$options$bru_method$rel_tol)

  pl_change <-
    ggplot2::ggplot(
      track_data %>%
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
          MaxLin = max(.data$new_linearisation - .data$new_linearisation.prev),
          MeanLin = mean(.data$new_linearisation - .data$new_linearisation.prev),
          MinLin = min(.data$new_linearisation - .data$new_linearisation.prev),
          .groups = "drop"
        )
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        .data$iteration,
        ymin = -MaxSD,
        ymax = MaxSD,
        fill = col_["Max"],
        lty = lty_["SD"]
      ),
      alpha = 0.1,
      na.rm = TRUE
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        .data$iteration,
        ymin = -MeanSD,
        ymax = MeanSD,
        fill = col_["Mean"],
        lty = lty_["SD"]
      ),
      alpha = 0.1,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        MaxMode,
        col = col_["Max"],
        lty = lty_["Mode"]
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        MinMode,
        col = col_["Min"],
        lty = lty_["Mode"]
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        MeanMode,
        col = col_["Mean"],
        lty = lty_["Mode"]
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        MaxLin,
        col = col_["Max"],
        lty = lty_["Lin"]
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        MinLin,
        col = col_["Min"],
        lty = lty_["Lin"]
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        MeanLin,
        col = col_["Mean"],
        lty = lty_["Lin"]
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        MaxSD,
        col = col_["Max"],
        lty = lty_["SD"]
      ),
      alpha = 0.5,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        MinSD,
        col = col_["Min"],
        lty = lty_["SD"]
      ),
      alpha = 0.5,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        MeanSD,
        col = col_["Mean"],
        lty = lty_["SD"]
      ),
      alpha = 0.5,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        -MaxSD,
        col = col_["Max"],
        lty = lty_["SD"]
      ),
      alpha = 0.5,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        -MinSD,
        col = col_["Min"],
        lty = lty_["SD"]
      ),
      alpha = 0.5,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        .data$iteration,
        -MeanSD,
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
        pl_relative_change + ggplot2::guides(linetype = "none", color = "none") |
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

  ggplot2::ggplot(x[["bru_timings"]]) +
    ggplot2::geom_point(ggplot2::aes(
      .data$Iteration,
      .data$Time,
      col = .data$Task,
    )) +
    ggplot2::geom_line(ggplot2::aes(
      .data$Iteration,
      .data$Time,
      col = .data$Task,
    )) +
    ggplot2::scale_y_continuous()
}
