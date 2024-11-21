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
#'   mesh1D <- fm_mesh_1d(x, boundary = "free")
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
bincount <- function(result, predictor, observations, breaks, nint = 20,
                     probs = c(0.025, 0.5, 0.975), ...) {
  # Sort probabilities
  probs <- sort(probs)

  # Filter out observations outside bins
  observations <- observations[(observations >= min(breaks)) &
    (observations <= max(breaks))]

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
      breaks[rep(seq_len(nbins), each = nsub)] * (1 - u) +
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
  smp <- generate(result,
    newdata = points,
    formula = predictor,
    ...,
    format = "matrix"
  )

  # Integrate per bin
  qq <- matrix(0.0, nrow = nbins, ncol(smp))
  for (k in seq_len(ncol(smp))) {
    qq[, k] <- fmesher::fm_block_logsumexp_eval(points$.block,
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
    zz <- apply(zz, MARGIN = 1, mean)
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
      data = pint, mapping = ggplot2::aes(
        x = .data[["mid"]],
        y = .data[[mdname]]
      ),
      shape = 95, size = 3, color = "blue"
    ) +
    ggplot2::geom_point(
      data = pint, mapping = ggplot2::aes(
        x = .data[["mid"]],
        y = .data[["counts"]]
      ),
      shape = 20, size = 2
    ) +
    ggplot2::xlab(all.vars(update.formula(predictor, ~.0))) +
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
#' @param samplers A SpatialPolygon object describing the area for which to
#'   compute the cumulative variance measure.
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
#' @example inst/examples/devel.cvmeasure.R

devel.cvmeasure <- function(joint, prediction1, prediction2, samplers = NULL,
                            mesh = NULL) {
  # Covariance
  joint$cov <- (joint[["sd"]]^2 - prediction1[["sd"]]^2 -
    prediction2[["sd"]]^2) / 2

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

  return(ret)
}
