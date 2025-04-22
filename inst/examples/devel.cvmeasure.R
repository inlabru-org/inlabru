\donttest{
if (bru_safe_inla() &&
  require("ggplot2", quietly = TRUE) &&
  require("patchwork", quietly = TRUE) &&
  require("sn", quietly = TRUE) &&
  require("terra", quietly = TRUE) &&
  require("sf", quietly = TRUE) &&
  require("RColorBrewer", quietly = TRUE) &&
  require("dplyr", quietly = TRUE) &&
  require("magrittr", quietly = TRUE)) {
  # Load Gorilla data

  gorillas <- gorillas_sf
  gorillas$gcov <- gorillas_sf_gcov()

  # Use RColorBrewer

  # Fit a model with two components:
  # 1) A spatial smooth SPDE
  # 2) A spatial covariate effect (vegetation)

  pcmatern <- INLA::inla.spde2.pcmatern(
    gorillas$mesh,
    prior.sigma = c(0.1, 0.01),
    prior.range = c(0.01, 0.01)
  )

  cmp <- geometry ~ 0 +
    vegetation(gorillas$gcov$vegetation, model = "factor_contrast") +
    spde(geometry, model = pcmatern) -
    Intercept(1)

  fit <- lgcp(
    cmp,
    gorillas$nests,
    samplers = gorillas$boundary,
    domain = list(geometry = gorillas$mesh),
    options = list(control.inla = list(int.strategy = "eb"))
  )

  # Predict SPDE and vegetation at a grid covering the domain of interest
  pred_loc <- fm_pixels(
    gorillas$mesh,
    mask = gorillas$boundary,
    dims = c(200, 200),
    format = "sf"
  )
  pred <- predict(
    fit,
    pred_loc,
    ~ list(
      joint = spde + vegetation,
      field = spde,
      veg = vegetation
    )
  )

  pred_collect <-
    rbind(
      pred$joint %>% dplyr::mutate(component = "joint"),
      pred$field %>% dplyr::mutate(component = "field"),
      pred$veg %>% dplyr::mutate(component = "veg")
    ) %>%
    dplyr::mutate(var = sd^2)

  # Plot component mean

  ggplot(pred_collect) +
    geom_tile(aes(geometry = geometry, fill = mean), stat = "sf_coordinates") +
    coord_equal() +
    facet_wrap(~component, nrow = 1) +
    theme(legend.position = "bottom")

  # Plot component variance

  ggplot(pred_collect) +
    geom_tile(aes(geometry = geometry, fill = var), stat = "sf_coordinates") +
    coord_equal() +
    facet_wrap(~component, nrow = 1) +
    theme(legend.position = "bottom")

  # Calculate variance and correlation measure

  vm <- devel.cvmeasure(pred$joint, pred$field, pred$veg)
  # Compute nominal relative variance contributions; note that these can be
  # greater than 100%!
  vm <- dplyr::mutate(
    vm,
    var1_rel = if_else(var1 <= 0, NA, var1 / var.joint),
    var2_rel = if_else(var2 <= 0, NA, var2 / var.joint)
  )
  lprange <- range(vm$var.joint, vm$var1, vm$var2)
  vm <- tidyr::pivot_longer(vm,
    cols = c(var.joint, var1, var2, cov, cor, var1_rel, var2_rel),
    names_to = "component",
    values_to = "value"
  )

  # Variance contribution of the components

  csc <- scale_fill_gradientn(
    colours = brewer.pal(9, "YlOrRd"),
    limits = lprange
  )

  vm_ <- dplyr::filter(vm, component %in% c("var.joint", "var1", "var2"))
  ggplot(vm_) +
    geom_tile(aes(geometry = geometry, fill = value),
      stat = "sf_coordinates"
    ) +
    csc +
    coord_equal() +
    facet_wrap(~component, nrow = 1) +
    theme(legend.position = "bottom")

  # Relative variance contribution of the components
  # When bo

  vm_ <- dplyr::filter(vm, component %in% c("var1_rel", "var2_rel"))
  ggplot(vm_) +
    geom_tile(aes(geometry = geometry, fill = value),
      stat = "sf_coordinates"
    ) +
    scale_fill_gradientn(
      colours = brewer.pal(9, "YlOrRd"),
      trans = "log"
    ) +
    coord_equal() +
    facet_wrap(~component, nrow = 1)

  # Where both relative contributions are larger than 1, the posterior
  # correlations are strongly negative.

  # Covariance and correlation of field and vegetation

  vm_cov <- dplyr::filter(vm, component %in% "cov")
  vm_cor <- dplyr::filter(vm, component %in% "cor")
  (ggplot(vm_cov) +
    geom_tile(aes(geometry = geometry, fill = value),
      stat = "sf_coordinates"
    ) +
    scale_fill_gradientn(
      colours = rev(brewer.pal(9, "RdBu")),
      limits = c(-1, 1) * max(abs(vm_cov$value), na.rm = TRUE)
    ) +
    coord_equal() +
    theme(legend.position = "bottom") +
    ggtitle("Covariances") |
    ggplot(vm_cor) +
      geom_tile(aes(geometry = geometry, fill = value),
        stat = "sf_coordinates"
      ) +
      scale_fill_gradientn(
        colours = rev(brewer.pal(9, "RdBu")),
        limits = c(-1, 1) * max(abs(vm_cor$value), na.rm = TRUE)
      ) +
      coord_equal() +
      theme(legend.position = "bottom") +
      ggtitle("Correlations")
  )

  # Variance and correlation integrated over space

  vrt <- fm_vertices(gorillas$mesh, format = "sf")
  pred_vrt <- predict(
    fit,
    vrt,
    ~ list(
      joint = spde + vegetation,
      field = spde,
      veg = vegetation
    )
  )

  vm.int <- devel.cvmeasure(
    pred_vrt$joint,
    pred_vrt$field,
    pred_vrt$veg,
    samplers = fm_int(gorillas$mesh, gorillas$boundary),
    mesh = gorillas$mesh
  )
  vm.int
}
}
