# library(tidyverse)
# library(patchwork)

#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#' @importFrom rlang .env
#' @import patchwork

make_track_plots <- function(fit) {
  if (!requireNamespace("dplyr", quietly = TRUE) ||
    !requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("magrittr", quietly = TRUE)) {
    stop("One or more of 'dplyr', 'ggplot2', 'magrittr' not installed, but are needed by make_track_plots()")
  }
  track_data <-
    fit$bru_iinla$track %>%
    dplyr::left_join(fit$bru_iinla$track %>%
      dplyr::filter(.data$iteration == max(.data$iteration)) %>%
      dplyr::rename(
        mode_end = .data$mode,
        new_linearisation_end = .data$new_linearisation,
        sd_end = .data$sd
      ),
    by = c("effect", "index"),
    ) %>%
    dplyr::mutate(iteration = .data$iteration.x) %>%
    dplyr::left_join(fit$bru_iinla$track %>%
      dplyr::filter(.data$iteration == 1) %>%
      dplyr::rename(
        mode_start = .data$mode,
        new_linearisation_start = .data$new_linearisation,
        sd_start = .data$sd
      ),
    by = c("effect", "index"),
    ) %>%
    dplyr::mutate(iteration = .data$iteration.x)

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
      lty = lty_["Mode"]
    )) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      .data$new_linearisation,
      col = .data$index,
      group = factor(.data$index),
      lty = lty_["Lin"]
    )) +
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
    ggplot2::ggplot(track_data) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      .data$mode - .data$new_linearisation,
      col = .data$index,
      group = factor(.data$index),
      lty = lty_["Mode-Lin"]
    )) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      -.data$sd,
      col = .data$index,
      group = factor(.data$index),
      lty = lty_["SD"]
    )) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      .data$sd,
      col = .data$index,
      group = factor(.data$index),
      lty = lty_["SD"]
    )) +
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
      (.data$mode - .data$mode.prev) / .data$sd,
      col = .data$index,
      group = factor(.data$index),
      lty = lty_["Mode"]
    )) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      (.data$new_linearisation - .data$new_linearisation.prev) / .data$sd,
      col = .data$index,
      group = factor(.data$index),
      lty = lty_["Lin"]
    )) +
    pl_theme_norm +
    ggplot2::ggtitle("Change / sd")

  pl6 <-
    ggplot2::ggplot(track_data) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      (.data$mode - .data$mode.prev),
      col = .data$index,
      group = factor(.data$index),
      lty = lty_["Mode"]
    )) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      (.data$new_linearisation - .data$new_linearisation.prev),
      col = .data$index,
      group = factor(.data$index),
      lty = lty_["Lin"]
    )) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      .data$sd,
      col = .data$index,
      group = factor(.data$index),
      lty = lty_["SD"]
    )) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      -.data$sd,
      col = .data$index,
      group = factor(.data$index),
      lty = lty_["SD"]
    )) +
    pl_theme_rel +
    ggplot2::ggtitle("Change & sd")

  pl_combined <-
    ((pl1 |
      pl3 + ggplot2::guides(linetype = "none")) /
      (pl5 + ggplot2::guides(linetype = "none") |
        pl6 + ggplot2::guides(linetype = "none"))) +
      patchwork::plot_layout(guides = "collect") &
      ggplot2::theme(legend.position = "bottom")


  list(pl1, pl2, pl3, pl4, pl5, pl6,
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
