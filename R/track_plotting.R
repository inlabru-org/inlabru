# library(tidyverse)
# library(patchwork)

make_track_plots <- function(fit) {
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
      ggplot2::facet_wrap(vars(effect), scales = "free"),
      ggplot2::guides(color = "none")
    )

  pl1 <-
    ggplot2::ggplot(track_data) +
    ggplot2::geom_line(aes(
      iteration,
      mode,
      col = index,
      group = factor(index),
      lty = "Mode"
    )) +
    ggplot2::geom_line(aes(
      iteration,
      new_linearisation,
      col = index,
      group = factor(index),
      lty = "Lin"
    )) +
    pl_theme +
    ggtitle("Tracks")
  pl2 <-
    ggplot2::ggplot(track_data) +
    ggplot2::geom_line(aes(
      iteration,
      mode - mode_end,
      col = index,
      group = factor(index),
      lty = "Mode"
    )) +
    ggplot2::geom_line(aes(
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
    ggplot2::geom_line(aes(
      iteration,
      mode - new_linearisation,
      col = index,
      group = factor(index),
      lty = "Mode-Lin"
    )) +
    ggplot2::geom_line(aes(
      iteration,
      -sd,
      col = index,
      group = factor(index),
      lty = "SD"
    )) +
    ggplot2::geom_line(aes(
      iteration,
      sd,
      col = index,
      group = factor(index),
      lty = "SD"
    )) +
    pl_theme +
    ggtitle("Mode - Lin")
  pl4 <-
    ggplot2::ggplot(track_data) +
    ggplot2::geom_line(aes(
      iteration,
      (mode - mode_end) / sd,
      col = index,
      group = factor(index),
      lty = "Mode"
    )) +
    ggplot2::geom_line(aes(
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
