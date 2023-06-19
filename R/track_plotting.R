#' @importFrom magrittr "%>%"
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

  sc <- ggplot2::scale_linetype_discrete(
    name = "Lines",
    breaks = names(lty_),
    labels = labels(lty_),
    drop = FALSE
  )

  pl_theme_abs <-
    list(
      ggplot2::facet_wrap(ggplot2::vars(.data$effect), scales = "free"),
      sc,
      ggplot2::guides(color = "none")
    )
  pl_theme_rel <-
    list(
      ggplot2::facet_wrap(ggplot2::vars(.data$effect), scales = "free"),
      sc,
      ggplot2::guides(color = "none")
    )
  pl_theme_norm <-
    list(
      ggplot2::facet_wrap(ggplot2::vars(.data$effect)),
      sc,
      ggplot2::guides(color = "none")
    )

  pl1 <-
    ggplot2::ggplot(track_data) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      .data$mode,
      col = .data$index,
      group = factor(.data$index),
      lty = lty_["Mode"],
      alpha = .data$alpha_value
    )) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      .data$new_linearisation,
      col = .data$index,
      group = factor(.data$index),
      lty = lty_["Lin"],
      alpha = .data$alpha_value * alpha_value_lin_scale
    )) +
    ggplot2::guides(alpha = "none") +
    pl_theme_abs +
    ggplot2::ggtitle("Tracks")

  pl2 <-
    ggplot2::ggplot(track_data) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      .data$mode - .data$mode_end,
      col = .data$index,
      group = factor(.data$index),
      lty = lty_["Mode"]
    )) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      .data$new_linearisation - .data$mode_end,
      col = .data$index,
      group = factor(.data$index),
      lty = lty_["Lin"]
    )) +
    pl_theme_rel +
    ggplot2::ggtitle("Rel end mode")
  pl3 <-
    ggplot2::ggplot(
      track_data %>%
        dplyr::filter(
          is.finite(.data$sd),
          is.finite(.data$new_linearisation)
        )
    ) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      .data$mode - .data$new_linearisation,
      col = .data$index,
      group = factor(.data$index),
      lty = lty_["Mode-Lin"],
      alpha = .data$alpha_value
    )) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      -.data$sd,
      col = .data$index,
      group = factor(.data$index),
      lty = lty_["SD"],
      alpha = .data$alpha_value
    )) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      .data$sd,
      col = .data$index,
      group = factor(.data$index),
      lty = lty_["SD"],
      alpha = .data$alpha_value
    )) +
    ggplot2::guides(alpha = "none") +
    pl_theme_rel +
    ggplot2::ggtitle("Mode - Lin")
  pl4 <-
    ggplot2::ggplot(track_data) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      (.data$mode - .data$mode_end) / .data$sd,
      col = .data$index,
      group = factor(.data$index),
      lty = lty_["Mode"]
    )) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      (.data$new_linearisation - .data$mode_end) / .data$sd,
      col = .data$index,
      group = factor(.data$index),
      lty = lty_["Lin"]
    )) +
    pl_theme_norm +
    ggplot2::ggtitle("Rel end mode / sd")

  track_data_prev <- dplyr::filter(track_data, .data$iteration < max(.data$iteration))
  track_data_prev$iteration <- track_data_prev$iteration + 1
  track_data <- dplyr::left_join(track_data, track_data_prev,
    by = c("effect", "iteration", "index"),
    suffix = c("", ".prev")
  )

  pl5 <-
    ggplot2::ggplot(track_data) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      abs(.data$mode - .data$mode.prev) / .data$sd,
      col = .data$index,
      group = factor(.data$index),
      lty = lty_["Mode"]
    )) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      abs(.data$new_linearisation - .data$new_linearisation.prev) / .data$sd,
      col = .data$index,
      group = factor(.data$index),
      lty = lty_["Lin"]
    )) +
    pl_theme_norm +
    ggplot2::scale_y_log10() +
    ggplot2::ggtitle("|Change| / sd")

  pl5b <-
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
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      .data$MaxMode,
      lty = lty_["Mode"]
    )) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      .data$MeanMode,
      lty = lty_["Mode"]
    )) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      .data$MaxLin,
      lty = lty_["Lin"]
    )) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      .data$MeanLin,
      lty = lty_["Lin"]
    )) +
    pl_theme_norm +
    ggplot2::scale_y_log10() +
    ggplot2::ggtitle("|Change| / sd (Max and Mean)") +
    ggplot2::geom_hline(yintercept = fit$bru_info$options$bru_method$rel_tol)

  pl6 <-
    ggplot2::ggplot(
      track_data %>%
        dplyr::filter(is.finite(.data$sd))
    ) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      (.data$mode - .data$mode.prev),
      col = .data$index,
      group = factor(.data$index),
      lty = lty_["Mode"],
      alpha = .data$alpha_value
    )) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      (.data$new_linearisation - .data$new_linearisation.prev),
      col = .data$index,
      group = factor(.data$index),
      lty = lty_["Lin"],
      alpha = .data$alpha_value * alpha_value_lin_scale
    )) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      .data$sd,
      col = .data$index,
      group = factor(.data$index),
      lty = lty_["SD"],
      alpha = .data$alpha_value
    )) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      -.data$sd,
      col = .data$index,
      group = factor(.data$index),
      lty = lty_["SD"],
      alpha = .data$alpha_value
    )) +
    ggplot2::guides(alpha = "none") +
    pl_theme_rel +
    ggplot2::ggtitle("Change & sd")

  pl_combined <-
    ((pl1 |
      pl3 + ggplot2::guides(linetype = "none")) /
      (pl5b + ggplot2::guides(linetype = "none") |
        pl6 + ggplot2::guides(linetype = "none"))) +
      patchwork::plot_layout(guides = "collect") &
      ggplot2::theme(legend.position = "bottom")


  list(pl1, pl2, pl3, pl4, pl5, pl5b, pl6,
    default = pl_combined
  )
}

# (
#   (
#     pl1 |
#      pl2
#  ) / (
#    pl3 | pl4
#  )
# )
