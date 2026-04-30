#' 1D LGCP bin count simulation and comparison with data
#'
#' A common procedure of analyzing the distribution of 1D points is to chose a
#' binning and plot the data's histogram with respect to this binning. This
#' function compares the counts that the histogram calculates to simulations
#' from a 1D log Gaussian Cox process conditioned on the number of data samples.
#' For each bin this results in a median number of counts as well as a
#' confidence interval. If the LGCP is a plausible model for the observed points
#' then most of the histogram counts (number of points within a bin) should be
#' within the confidence intervals. Note that a proper comparison  is a multiple
#' testing problem which the function does not solve for you.
#'
#'
#' @export
#' @importFrom graphics hist
#' @param result A result object from a [bru()] or [lgcp()] call
#' @param predictor A formula describing the prediction of a 1D LGCP via
#'   [predict()].
#' @param observations A vector of observed values
#' @param breaks A vector of bin boundaries
#' @param nint Number of integration intervals per bin. Increase this if the
#'   bins are wide and the LGCP is not smooth.
#' @param probs numeric vector of probabilities with values in `[0,1]`
#' @param \dots arguments passed on to [predict.bru()]
#' @return An `data.frame` with a ggplot attribute `ggp`
#' @importFrom rlang .data
#' @examples
#' \dontrun{
#' if (require(ggplot2) && require(fmesher) && bru_safe_inla()) {
#'   # Load a point pattern
#'   data(Poisson2_1D)
#'
#'   # Take a look at the point (and frequency) data
#'
#'   ggplot(pts2) +
#'     geom_histogram(
#'       aes(x = x),
#'       binwidth = 55 / 20,
#'       boundary = 0,
#'       fill = NA,
#'       color = "black"
#'     ) +
#'     geom_point(aes(x), y = 0, pch = "|", cex = 4) +
#'     coord_fixed(ratio = 1)
#'
#'   # Fit an LGCP model
#'   x <- seq(0, 55, length.out = 50)
#'   mesh1D <- fmesher::fm_mesh_1d(x, boundary = "free")
#'   matern <- INLA::inla.spde2.pcmatern(mesh1D,
#'     prior.range = c(1, 0.01),
#'     prior.sigma = c(1, 0.01),
#'     constr = TRUE
#'   )
#'   mdl <- x ~ spde1D(x, model = matern) + Intercept(1)
#'   fit.spde <- lgcp(mdl, pts2, domain = list(x = mesh1D))
#'
#'   # Calculate bin statistics
#'   bc <- bincount(
#'     result = fit.spde,
#'     observations = pts2,
#'     breaks = seq(0, max(pts2), length.out = 12),
#'     predictor = x ~ exp(spde1D + Intercept)
#'   )
#'
#'   # Plot them!
#'   attributes(bc)$ggp
#' }
#' }
#'
bincount <- function(
  result,
  predictor,
  observations,
  breaks,
  nint = 20,
  probs = c(0.025, 0.5, 0.975),
  ...
) {
  # Sort probabilities
  probs <- sort(probs)

  # Filter out observations outside bins
  observations <- observations[
    (observations >= min(breaks)) &
      (observations <= max(breaks))
  ]

  # Number of ...
  nobs <- length(observations) # observations
  nbins <- length(breaks) - 1 # bins

  # Mid points of bins
  mid <- breaks[1:(length(breaks) - 1)] + diff(breaks) / 2

  # Create integration points
  nsub <- nint
  u <- rep(
    (seq_len(nsub) - 1L) / (nsub + 1),
    nbins
  )
  domain_breaks <-
    c(
      breaks[rep(seq_len(nbins), each = nsub)] *
        (1 - u) +
        breaks[rep(seq_len(nbins) + 1, each = nsub)] * u,
      breaks[length(breaks)]
    )
  points <- fmesher::fm_int(
    samplers = cbind(
      breaks[1:(length(breaks) - 1)],
      breaks[2:(length(breaks))]
    ),
    domain = fmesher::fm_mesh_1d(domain_breaks),
    int.args = list(
      method = "stable",
      nsub1 = 1L
    ),
    name = as.character(predictor)[2]
  )

  # Sampler
  smp <- generate(
    result,
    newdata = points,
    formula = predictor,
    ...,
    format = "matrix"
  )

  # Integrate per bin
  qq <- matrix(0.0, nrow = nbins, ncol(smp))
  for (k in seq_len(ncol(smp))) {
    qq[, k] <- fmesher::fm_block_logsumexp_eval(
      points$.block,
      weights = points$weight,
      values = log(smp[, k]),
      log = TRUE
    )
    # Normalize bin probabilities
    qq[, k] <- qq[, k] - max(qq[, k])
    qq[, k] <- exp(qq[, k] - log(sum(exp(qq[, k]))))
  }

  # For each bin calculate predictive interval
  pint <- list()
  for (k in 1:nbins) {
    xx <- 0:(nobs + 1L)
    cdff <- function(p) pbinom(xx, size = nobs, prob = p)
    zz <- apply(qq[k, , drop = FALSE], MARGIN = 2, cdff)
    zz <- rowMeans(zz)
    pint[[k]] <- vapply(probs, function(pr) xx[sum(zz < pr) + 1L], 0.0)
  }
  pint <- data.frame(do.call(rbind, pint))
  colnames(pint) <- paste0("q", probs)
  mxname <- paste0("q", max(probs))
  miname <- paste0("q", min(probs))
  mdname <- paste0("q", probs[1 + floor(length(probs) / 2)])

  # Append more information
  pint <- cbind(
    data.frame(
      mid = mid,
      width = diff(breaks),
      counts = hist(observations, breaks = breaks, plot = FALSE)$counts
    ),
    pint
  )
  pint$inside <-
    (pint$counts >= pint[[miname]]) &
      (pint$counts <= pint[[mxname]])

  ggp <- ggplot2::ggplot() +
    ggplot2::geom_crossbar(
      data = pint,
      mapping = ggplot2::aes(
        x = .data[["mid"]],
        y = .data[[mdname]],
        ymin = .data[[miname]],
        ymax = .data[[mxname]],
        fill = .data[["inside"]],
        color = .data[["inside"]]
      ),
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = pint,
      mapping = ggplot2::aes(
        x = .data[["mid"]],
        y = .data[[mdname]]
      ),
      shape = 95,
      size = 3,
      color = "blue"
    ) +
    ggplot2::geom_point(
      data = pint,
      mapping = ggplot2::aes(
        x = .data[["mid"]],
        y = .data[["counts"]]
      ),
      shape = 20,
      size = 2
    ) +
    ggplot2::xlab(all.vars(update.formula(predictor, . ~ 0))) +
    ggplot2::ylab("count")

  attr(pint, "ggp") <- ggp

  # Return
  pint
}


#' Variance and correlations measures for prediction components
#'
#' Calculates local and integrated variance and correlation measures as
#' introduced by Yuan et al. (2017).
#'
#' @export
#' @param joint A joint `prediction` of two latent model components.
#' @param prediction1 A `prediction` of the first component.
#' @param prediction2 A `prediction` of the second component.
#' @param samplers An `sf` or `SpatialPolygon` object describing the area for
#'   which to compute the cumulative variance measure.
#' @param mesh The [fmesher::fm_mesh_2d] for which the prediction was performed
#'   (required for cumulative Vmeasure).
#' @return Variance and correlations measures.
#'
#' @references
#' Y. Yuan, F. E. Bachl, F. Lindgren, D. L. Brochers, J. B. Illian, S. T.
#' Buckland, H. Rue, T. Gerrodette. 2017. Point process models for
#' spatio-temporal distance sampling data from a large-scale survey of blue
#' whales. <https://arxiv.org/abs/1604.06013>
#'
#' @examples
#' \donttest{
#' if (bru_safe_inla() &&
#'   require("ggplot2", quietly = TRUE) &&
#'   require("patchwork", quietly = TRUE) &&
#'   requireNamespace("sn", quietly = TRUE) &&
#'   bru_safe_terra(quietly = TRUE) &&
#'   require("sf", quietly = TRUE) &&
#'   require("RColorBrewer", quietly = TRUE) &&
#'   require("dplyr", quietly = TRUE)) {
#'   # Load Gorilla data
#'
#'   gorillas <- gorillas_sf
#'   gorillas$gcov <- gorillas_sf_gcov()
#'
#'   # Use RColorBrewer
#'
#'   # Fit a model with two components:
#'   # 1) A spatial smooth SPDE
#'   # 2) A spatial covariate effect (vegetation)
#'
#'   pcmatern <- INLA::inla.spde2.pcmatern(
#'     gorillas$mesh,
#'     prior.sigma = c(0.1, 0.01),
#'     prior.range = c(0.01, 0.01)
#'   )
#'
#'   cmp <- geometry ~ 0 +
#'     vegetation(gorillas$gcov$vegetation, model = "factor_contrast") +
#'     spde(geometry, model = pcmatern) +
#'     Intercept(1)
#'
#'   fit <- lgcp(
#'     cmp,
#'     gorillas$nests,
#'     samplers = gorillas$boundary,
#'     domain = list(geometry = gorillas$mesh),
#'     options = list(control.inla = list(int.strategy = "eb"))
#'   )
#'
#'   # Predict SPDE and vegetation at a grid covering the domain of interest
#'   pred_loc <- fmesher::fm_pixels(
#'     gorillas$mesh,
#'     mask = gorillas$boundary,
#'     dims = c(200, 200),
#'     format = "sf"
#'   )
#'   pred <- predict(
#'     fit,
#'     pred_loc,
#'     ~ list(
#'       joint = spde + vegetation,
#'       field = spde,
#'       veg = vegetation
#'     )
#'   )
#'
#'   pred_collect <-
#'     rbind(
#'       pred$joint |> dplyr::mutate(component = "joint"),
#'       pred$field |> dplyr::mutate(component = "field"),
#'       pred$veg |> dplyr::mutate(component = "veg")
#'     ) |>
#'     dplyr::mutate(var = sd^2)
#'
#'   # Plot component mean
#'
#'   ggplot(pred_collect) +
#'     geom_tile(aes(geometry = geometry, fill = mean),
#'       stat = "sf_coordinates"
#'     ) +
#'     coord_equal() +
#'     facet_wrap(~component, nrow = 1) +
#'     theme(legend.position = "bottom")
#'
#'   # Plot component variance
#'
#'   ggplot(pred_collect) +
#'     geom_tile(aes(geometry = geometry, fill = var),
#'       stat = "sf_coordinates"
#'     ) +
#'     coord_equal() +
#'     facet_wrap(~component, nrow = 1) +
#'     theme(legend.position = "bottom")
#'
#'   # Calculate variance and correlation measure
#'
#'   vm <- devel.cvmeasure(pred$joint, pred$field, pred$veg)
#'   # Compute nominal relative variance contributions; note that these can be
#'   # greater than 100%!
#'   vm <- dplyr::mutate(
#'     vm,
#'     var1_rel = if_else(var1 <= 0, NA, var1 / var.joint),
#'     var2_rel = if_else(var2 <= 0, NA, var2 / var.joint)
#'   )
#'   lprange <- range(vm$var.joint, vm$var1, vm$var2)
#'   vm <- tidyr::pivot_longer(vm,
#'     cols = c(var.joint, var1, var2, cov, cor, var1_rel, var2_rel),
#'     names_to = "component",
#'     values_to = "value"
#'   )
#'
#'   # Variance contribution of the components
#'
#'   csc <- scale_fill_gradientn(
#'     colours = brewer.pal(9, "YlOrRd"),
#'     limits = lprange
#'   )
#'
#'   vm_ <- dplyr::filter(vm, component %in% c("var.joint", "var1", "var2"))
#'   ggplot(vm_) +
#'     geom_tile(aes(geometry = geometry, fill = value),
#'       stat = "sf_coordinates"
#'     ) +
#'     csc +
#'     coord_equal() +
#'     facet_wrap(~component, nrow = 1) +
#'     theme(legend.position = "bottom")
#'
#'   # Relative variance contribution of the components
#'   # When bo
#'
#'   vm_ <- dplyr::filter(vm, component %in% c("var1_rel", "var2_rel"))
#'   ggplot(vm_) +
#'     geom_tile(aes(geometry = geometry, fill = value),
#'       stat = "sf_coordinates"
#'     ) +
#'     scale_fill_gradientn(
#'       colours = brewer.pal(9, "YlOrRd"),
#'       trans = "log"
#'     ) +
#'     coord_equal() +
#'     facet_wrap(~component, nrow = 1)
#'
#'   # Where both relative contributions are larger than 1, the posterior
#'   # correlations are strongly negative.
#'
#'   # Covariance and correlation of field and vegetation
#'
#'   vm_cov <- dplyr::filter(vm, component %in% "cov")
#'   vm_cor <- dplyr::filter(vm, component %in% "cor")
#'   (ggplot(vm_cov) +
#'     geom_tile(aes(geometry = geometry, fill = value),
#'       stat = "sf_coordinates"
#'     ) +
#'     scale_fill_gradientn(
#'       colours = rev(brewer.pal(9, "RdBu")),
#'       limits = c(-1, 1) * max(abs(vm_cov$value), na.rm = TRUE)
#'     ) +
#'     coord_equal() +
#'     theme(legend.position = "bottom") +
#'     ggtitle("Covariances") |
#'     ggplot(vm_cor) +
#'       geom_tile(aes(geometry = geometry, fill = value),
#'         stat = "sf_coordinates"
#'       ) +
#'       scale_fill_gradientn(
#'         colours = rev(brewer.pal(9, "RdBu")),
#'         limits = c(-1, 1) * max(abs(vm_cor$value), na.rm = TRUE)
#'       ) +
#'       coord_equal() +
#'       theme(legend.position = "bottom") +
#'       ggtitle("Correlations")
#'   )
#'
#'   # Variance and correlation integrated over space
#'
#'   vrt <- fmesher::fm_vertices(gorillas$mesh, format = "sf")
#'   pred_vrt <- predict(
#'     fit,
#'     vrt,
#'     ~ list(
#'       joint = spde + vegetation,
#'       field = spde,
#'       veg = vegetation
#'     )
#'   )
#'
#'   vm.int <- devel.cvmeasure(
#'     pred_vrt$joint,
#'     pred_vrt$field,
#'     pred_vrt$veg,
#'     samplers = fmesher::fm_int(gorillas$mesh, gorillas$boundary),
#'     mesh = gorillas$mesh
#'   )
#'   vm.int
#' }
#' }
#'
devel.cvmeasure <- function(
  joint,
  prediction1,
  prediction2,
  samplers = NULL,
  mesh = NULL
) {
  # Covariance
  joint$cov <- (joint[["sd"]]^2 -
    prediction1[["sd"]]^2 -
    prediction2[["sd"]]^2) /
    2

  # Correlation
  corr <- function(joint, a, b) {
    ((joint - a - b) / (2 * sqrt(a * b)))
  }
  cor <- corr(joint[["sd"]]^2, prediction1[["sd"]]^2, prediction2[["sd"]]^2)
  if (any((cor > 1) | (cor < (-1)), na.rm = TRUE)) {
    cor[(cor > 1) | (cor < (-1))] <- NA
  }
  joint$cor <- cor

  if (!is.null(samplers)) {
    #
    # PRESENTED TO YOU BY HACKY McHACKERSON
    #

    wips <- fm_int(mesh, samplers)
    A <- fm_basis(mesh, loc = wips)

    weights <- wips$weight
    weights <- weights / sum(weights)
    vj <- sum(as.vector(A %*% joint[["sd"]]^2) * weights)
    v1 <- sum(as.vector(A %*% prediction1[["sd"]]^2) * weights)
    v2 <- sum(as.vector(A %*% prediction2[["sd"]]^2) * weights)
    cr <- corr(vj, v1, v2)
    ret <- data.frame(var.joint = vj, var1 = v1, var2 = v2, cor = cr)
  } else {
    tmp <- attributes(joint)
    ret <- joint[, c("cov", "cor")]
    ret$var.joint <- joint[["sd"]]^2
    ret$var1 <- prediction1[["sd"]]^2
    ret$var2 <- prediction2[["sd"]]^2
    attr(ret, "misc") <- tmp[["misc"]]
    attr(ret, "type") <- tmp[["type"]]
  }

  ret
}
