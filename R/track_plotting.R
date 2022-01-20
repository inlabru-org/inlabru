# library(tidyverse)
# library(patchwork)


make_track_plots <- function(fit) {
  if (!requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("magrittr", quietly = TRUE)) {
    stop("One or more of 'dplyr', 'ggplot2', 'magrittr' not installed, but are needed by make_track_plots()")
  }
  track_data <-
    fit$bru_iinla$track %>%
    dplyr::left_join(fit$bru_iinla$track %>%
      dplyr::filter(iteration == max(iteration)) %>%
      dplyr::rename(
        mode_end = mode,
        new_linearisation_end = new_linearisation,
        sd_end = sd
      ),
    by = c("effect", "index"),
    ) %>%
    dplyr::mutate(iteration = iteration.x) %>%
    dplyr::left_join(fit$bru_iinla$track %>%
      dplyr::filter(iteration == 1) %>%
      dplyr::rename(
        mode_start = mode,
        new_linearisation_start = new_linearisation,
        sd_start = sd
      ),
    by = c("effect", "index"),
    ) %>%
    dplyr::mutate(iteration = iteration.x)

  pl_theme <-
    list(
      ggplot2::facet_wrap(ggplot2::vars(effect), scales = "free"),
      ggplot2::guides(color = "none")
    )

  pl1 <-
    ggplot2::ggplot(track_data) +
    ggplot2::geom_line(ggplot2::aes(
      iteration,
      mode,
      col = index,
      group = factor(index),
      lty = "Mode"
    )) +
    ggplot2::geom_line(ggplot2::aes(
      iteration,
      new_linearisation,
      col = index,
      group = factor(index),
      lty = "Lin"
    )) +
    pl_theme +
    ggplot2::ggtitle("Tracks")
  pl2 <-
    ggplot2::ggplot(track_data) +
    ggplot2::geom_line(ggplot2::aes(
      iteration,
      mode - mode_end,
      col = index,
      group = factor(index),
      lty = "Mode"
    )) +
    ggplot2::geom_line(ggplot2::aes(
      iteration,
      new_linearisation - mode_end,
      col = index,
      group = factor(index),
      lty = "Lin"
    )) +
    pl_theme +
    ggplot2::ggtitle("Rel end mode")
  pl3 <-
    ggplot2::ggplot(track_data) +
    ggplot2::geom_line(ggplot2::aes(
      iteration,
      mode - new_linearisation,
      col = index,
      group = factor(index),
      lty = "Mode-Lin"
    )) +
    ggplot2::geom_line(ggplot2::aes(
      iteration,
      -sd,
      col = index,
      group = factor(index),
      lty = "SD"
    )) +
    ggplot2::geom_line(ggplot2::aes(
      iteration,
      sd,
      col = index,
      group = factor(index),
      lty = "SD"
    )) +
    pl_theme +
    ggplot2::ggtitle("Mode - Lin")
  pl4 <-
    ggplot2::ggplot(track_data) +
    ggplot2::geom_line(ggplot2::aes(
      iteration,
      (mode - mode_end) / sd,
      col = index,
      group = factor(index),
      lty = "Mode"
    )) +
    ggplot2::geom_line(ggplot2::aes(
      iteration,
      (new_linearisation - mode_end) / sd,
      col = index,
      group = factor(index),
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
