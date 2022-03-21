# library(tidyverse)
# library(patchwork)

#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#' @importFrom rlang .env

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

  pl_theme <-
    list(
      ggplot2::facet_wrap(ggplot2::vars(.data$effect), scales = "free"),
      ggplot2::guides(color = "none")
    )

  pl1 <-
    ggplot2::ggplot(track_data) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      .data$mode,
      col = .data$index,
      group = factor(.data$index),
      lty = "Mode"
    )) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      .data$new_linearisation,
      col = .data$index,
      group = factor(.data$index),
      lty = "Lin"
    )) +
    pl_theme +
    ggplot2::ggtitle("Tracks")
  pl2 <-
    ggplot2::ggplot(track_data) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      .data$mode - .data$mode_end,
      col = .data$index,
      group = factor(.data$index),
      lty = "Mode"
    )) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      .data$new_linearisation - .data$mode_end,
      col = .data$index,
      group = factor(.data$index),
      lty = "Lin"
    )) +
    pl_theme +
    ggplot2::ggtitle("Rel end mode")
  pl3 <-
    ggplot2::ggplot(track_data) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      .data$mode - .data$new_linearisation,
      col = .data$index,
      group = factor(.data$index),
      lty = "Mode-Lin"
    )) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      -.data$sd,
      col = .data$index,
      group = factor(.data$index),
      lty = "SD"
    )) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      .data$sd,
      col = .data$index,
      group = factor(.data$index),
      lty = "SD"
    )) +
    pl_theme +
    ggplot2::ggtitle("Mode - Lin")
  pl4 <-
    ggplot2::ggplot(track_data) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      (.data$mode - .data$mode_end) / .data$sd,
      col = .data$index,
      group = factor(.data$index),
      lty = "Mode"
    )) +
    ggplot2::geom_line(ggplot2::aes(
      .data$iteration,
      (.data$new_linearisation - .data$mode_end) / .data$sd,
      col = .data$index,
      group = factor(.data$index),
      lty = "Lin"
    )) +
    pl_theme +
    ggplot2::ggtitle("Rel end mode / sd")

  list(pl1, pl2, pl3, pl4)
}

# (
#   (
#     pl1 |
#      pl2
#  ) / (
#    pl3 | pl4
#  )
# )
