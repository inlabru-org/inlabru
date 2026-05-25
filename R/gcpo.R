# gcpo.R
#
# Functions for Group Cross-Validation Predictive Ordering (GCPO).
#
# These functions support INLA's control.gcpo machinery for computing
# leave-group-out cross-validation scores, and provide post-fit extraction
# of block-averaged GCPO scores.
#
# Internal helpers called from: bru_obs_family_cp(), bru_obs_family_cp_sp(),
#                                iinla()  (all in bru.inference.R)
# User-facing:                   bru_block_gcpo(), bru_gcpo_table()


# Benchmark results:
# > x <- sample(1:1000, size = 1000000, replace = TRUE)
# > bench::mark(
#     R = convert_group_cv_blocks_to_friends_list(x, method = "R"),
#     C = convert_group_cv_blocks_to_friends_list(x, method = "C")
#   )
# expression      min   median `itr/sec`
# R             4.23s    4.23s     0.236
# C            97.1ms  97.97ms    10.1
convert_group_cv_blocks_to_friends_list <- function(group_cv_block,
                                                    method = "C") {
  # Would like:
  # group_cv_friends <- lapply(seq_len(max(ips$.block)), function(i) {
  #   which(group_cv_block == i)
  # })
  # Current inla.group.cv interface (2025-08-28)
  # Inefficient implementation:
  # group_cv_friends <- lapply(seq_along(group_cv_block), function(i) {
  #   which(group_cv_block == group_cv_block[i])
  # })
  # Faster version:
  if (method == "R") {
    grp_index <- seq_len(max(group_cv_block))
    group_cv_friends_groups <- lapply(grp_index, function(i) {
      which(group_cv_block == i)
    })
    group_cv_friends <- lapply(group_cv_block, function(i) {
      group_cv_friends_groups[[i]]
    })
  } else if (method == "C") {
    group_cv_friends <-
      inlabru_group_cv_block_conversion(
        as.integer(group_cv_block),
        as.integer(max(group_cv_block)),
        per_node = as.logical(TRUE)
      )
  }
  group_cv_friends
}


#' @title Utility functions for bru observation model GCPO
#'
#' @description
#' Extract and combine `INLA::control.gcpo` options from `bru_obs` and
#' `bru_obs_list` objects, adjusting observation indices for multi-likelihood
#' models where the response vector is the concatenation of all likelihoods.
#'
#' @param x A `bru_obs` or `bru_obs_list` object.
#' @param \dots Further arguments passed to submethods.
#' @returns
#'   * `bru_obs_control_gcpo()` returns a list with `INLA::control.gcpo`
#'   options, with predictor/response variable indices unified for
#'   multi-observation models. If a `bru_obs` model has a `control.gcpo = NULL`
#'   argument, an empty list is returned for that model.
#' @seealso [bru_obs()], [bru_block_gcpo()]
#' @export
#' @keywords internal
#' @rdname bru_obs_gcpo
#' @name bru_obs_gcpo
bru_obs_control_gcpo <- function(x, ...) {
  UseMethod("bru_obs_control_gcpo")
}


#' @param index_offset integer; offset to add to indices in `control.gcpo`,
#'   equal to the total number of response rows from all preceding likelihoods.
#' @param index_length integer; number of response rows for this likelihood.
#'   Defaults to `bru_response_size(x)`.
#' @param force_weights logical; if `TRUE`, ensure that `control.gcpo$weights`
#'   is populated. Required when any likelihood in a `bru_obs_list` has
#'   non-NULL `control.gcpo$weights`.
#' @export
#' @rdname bru_obs_gcpo
bru_obs_control_gcpo.bru_obs <- function(x,
                                         index_offset,
                                         index_length = bru_response_size(x),
                                         force_weights,
                                         ...) {
  c.gcpo <- x[["control.gcpo"]]
  if (is.null(c.gcpo)) {
    return(list())
  }
  for (nm in intersect(
    c("groups", "selection", "group.selection", "friends"),
    names(c.gcpo)
  )) {
    if (nm == "groups") {
      c.gcpo[[nm]] <- lapply(c.gcpo[[nm]], function(v) {
        v$idx <- v$idx + index_offset
        v
      })
    } else {
      c.gcpo[[nm]] <- lapply(c.gcpo[[nm]], function(v) v + index_offset)
    }
  }
  if (!("friends" %in% names(c.gcpo)) && !("groups" %in% names(c.gcpo))) {
    c.gcpo[["friends"]] <- as.list(index_offset + seq_len(index_length))
  }
  if (force_weights && !("weights" %in% names(c.gcpo))) {
    c.gcpo[["weights"]] <- rep(1.0, index_length)
  } else if ("weights" %in% names(c.gcpo)) {
    if (length(c.gcpo[["weights"]]) != index_length) {
      stop(glue::glue(
        "Length of control.gcpo$weights ({length(c.gcpo[['weights']])}) ",
        "does not match response length ({index_length})"
      ))
    }
  }
  c.gcpo
}


#' @param control.gcpo list of `INLA::control.gcpo` default options to merge
#'   with per-likelihood settings.
#' @export
#' @rdname bru_obs_gcpo
bru_obs_control_gcpo.bru_obs_list <- function(x,
                                              control.gcpo = NULL,
                                              ...) {
  response_sizes <- bru_response_size(x)

  any_element <- vapply(
    c("groups", "selection", "group.selection", "friends", "weights"),
    function(nm) {
      any(vapply(x, function(lh) !is.null(lh[["control.gcpo"]][[nm]]), TRUE))
    },
    TRUE
  )
  all_element <- vapply(
    c("groups", "selection", "group.selection", "friends", "weights"),
    function(nm) {
      all(vapply(x, function(lh) !is.null(lh[["control.gcpo"]][[nm]]), TRUE))
    },
    TRUE
  )

  c.gcpo <- lapply(
    seq_along(x),
    function(k) {
      bru_obs_control_gcpo(
        x[[k]],
        index_offset = sum(response_sizes[seq_len(k - 1)]),
        index_length = response_sizes[k],
        force_weights = any_element["weights"]
      )
    }
  )

  # If given in one model, must be given in all models
  c.gcpo.combined <- list()
  for (nm in c("groups", "selection", "group.selection")) {
    if (any_element[nm]) {
      if (!all_element[nm]) {
        stop(glue::glue(
          "control.gcpo${nm} given in some, but not all, observation models"
        ))
      }
      c.gcpo.combined[[nm]] <-
        do.call("c", lapply(c.gcpo, function(lh) lh[[nm]]))
    }
  }

  # Combine friends only when present (not combined with groups)
  if (any(vapply(c.gcpo, function(lh) !is.null(lh[["friends"]]), TRUE))) {
    c.gcpo.combined[["friends"]] <-
      do.call("c", lapply(c.gcpo, function(lh) lh[["friends"]]))
  }
  if (any_element[["weights"]]) {
    c.gcpo.combined[["weights"]] <-
      unlist(lapply(c.gcpo, function(lh) lh[["weights"]]))
  }

  merge_list_check <- function(a, b, exclude = character(0)) {
    for (nm in setdiff(names(b), exclude)) {
      if (nm %in% names(a)) {
        if (identical(a[[nm]], b[[nm]])) {
          next
        }
        stop(glue::glue(
          "Cannot merge control.gcpo lists: ",
          "both contain element '{nm}' with conflicting content."
        ))
      }
      a[[nm]] <- b[[nm]]
    }
    a
  }

  exc <- names(c.gcpo.combined)
  if (!is.null(control.gcpo)) {
    c.gcpo.combined <- merge_list_check(
      c.gcpo.combined,
      control.gcpo,
      exclude = exc
    )
  }
  for (idx in seq_along(c.gcpo)) {
    c.gcpo.combined <- merge_list_check(
      c.gcpo.combined,
      c.gcpo[[idx]],
      exclude = exc
    )
  }

  c.gcpo.combined
}


# Path 1: BRU_block based extraction (lgcp / cp family)
bru_block_gcpo_from_block <- function(lh, gcpo_vec, offset) {
  blk_var <- lh[["response_data"]][["BRU_block"]]
  blocks <- unique(blk_var)
  gcpo <- vapply(
    blocks,
    function(b) mean(gcpo_vec[which(blk_var == b) + offset]),
    numeric(1)
  )
  list(blocks = blocks, gcpo = gcpo)
}


# Path 2: groups based extraction (all other families)
bru_block_gcpo_from_groups <- function(lh, gcpo_vec, offset) {
  n <- bru_response_size(lh)
  idx <- offset + seq_len(n)
  list(
    blocks = seq_len(n),
    gcpo   = gcpo_vec[idx]
  )
}


# Internal workhorse for a single bru fit.
# Dispatches per-likelihood to bru_block_gcpo_from_block() when BRU_block
# is present (cp/lgcp), and to bru_block_gcpo_from_groups() otherwise.
# This handles mixed joint models where one likelihood is cp and another
# is a standard family.
bru_block_gcpo_single <- function(fit) {
  fit <- bru_check_object_bru(fit)

  gcpo_vec <- fit[["gcpo"]][["gcpo"]]
  if (is.null(gcpo_vec)) {
    stop(
      "fit$gcpo$gcpo is NULL. Refit with ",
      "control.gcpo = list(enable = TRUE)."
    )
  }

  lhoods <- as_bru_obs_list(fit)
  nlhoods <- length(lhoods)

  # Cumulative offsets into the stacked GCPO vector
  n_per_lhood <- vapply(lhoods, function(lh) bru_response_size(lh), integer(1))
  index_offset <- c(0L, cumsum(n_per_lhood))

  # Sanitise likelihood names
  lhood_names <- names(lhoods)
  if (is.null(lhood_names) || all(is.na(lhood_names))) {
    lhood_names <- paste0("lhood", seq_len(nlhoods))
  }

  # Per-likelihood dispatch
  results <- lapply(seq_len(nlhoods), function(k) {
    lh <- lhoods[[k]]
    offset <- index_offset[k]
    if (!is.null(lh[["response_data"]][["BRU_block"]])) {
      bru_block_gcpo_from_block(lh, gcpo_vec, offset)
    } else {
      bru_block_gcpo_from_groups(lh, gcpo_vec, offset)
    }
  })
  names(results) <- lhood_names

  if (nlhoods == 1L) {
    return(list(
      blocks = setNames(list(results[[1L]]$blocks), lhood_names[1L]),
      gcpo   = results[[1L]]$gcpo
    ))
  }

  list(
    blocks = setNames(lapply(results, `[[`, "blocks"), lhood_names),
    gcpo   = setNames(lapply(results, `[[`, "gcpo"), lhood_names)
  )
}

#' Extract block-averaged GCPO scores from one or more fitted bru models
#'
#' @description
#' After fitting a model with `control.gcpo = list(enable = TRUE)`,
#' this function reads the raw per-observation GCPO vector from
#' `fit$gcpo$gcpo` and returns one score per block or observation group.
#'
#' Two cases are handled automatically:
#'
#' \describe{
#'   \item{Block-based (lgcp / `"cp"` family)}{When `BRU_block` is present
#'     in `response_data` (populated automatically by the `"cp"` family
#'     machinery from the samplers structure), INLA repeats the same GCPO
#'     value for every observation in a leave-out block. The mean within
#'     each block collapses the repeated values to one score per block.}
#'   \item{Groups-based (all other families)}{When `BRU_block` is absent,
#'     INLA builds groups internally (via `num.level.sets` or `friends`) or
#'     from a user-supplied `groups` list. One score per observation is
#'     returned directly from `fit$gcpo$gcpo`.}
#' }
#'
#' Scores are returned on the probability scale, consistent with
#' `fit$cpo$cpo`. For multi-likelihood models, cumulative row offsets
#' into the stacked INLA response vector are handled automatically.
#'
#' @param fit A fitted object of class `bru`, or a named list of such objects,
#'   or several such objects passed via `\dots`. Each must have been fitted
#'   with `control.gcpo = list(enable = TRUE)` so that `fit$gcpo$gcpo` is
#'   non-NULL.
#' @param \dots Additional named fitted `bru` objects.
#'
#' @return
#' For a single fit, a list with two elements:
#' \describe{
#'   \item{`blocks`}{Named list with one element per likelihood, each a
#'     vector of block labels (for `"cp"` models: unique `BRU_block` values;
#'     for other models: observation indices).}
#'   \item{`gcpo`}{GCPO scores on the probability scale, consistent with
#'     `fit$cpo$cpo`. A numeric vector for single-likelihood models; a named
#'     list of numeric vectors for multi-likelihood models.}
#' }
#' For multiple fits, a named list of such objects, one per fit.
#'
#' @examples
#' \donttest{
#' if (bru_safe_inla()) {
#'   # Block-based example (lgcp)
#'   cvpart <- cv_hex(
#'     gorillas_sf$boundary,
#'     cellsize = 0.5,
#'     n_group = 3,
#'     resolution = c(95, 80)
#'   )
#'   cvpart$block_ID <- seq_len(nrow(cvpart))
#'   cvpart$group <- NULL
#'
#'   nests <- gorillas_sf$nests
#'   a <- sf::st_intersects(nests, cvpart)
#'   nests$.block <- unlist(a)
#'
#'   fit1 <- lgcp(
#'     geometry ~ Intercept(1),
#'     data = nests,
#'     samplers = cvpart,
#'     domain = list(geometry = gorillas_sf$mesh),
#'     control.gcpo = list(enable = TRUE, type.cv = "joint")
#'   )
#'   result <- bru_block_gcpo(fit1)
#'   result$blocks
#'   result$gcpo
#'   sum(log(result$gcpo))
#'
#'   # Multiple models
#'   pcmatern <- INLA::inla.spde2.pcmatern(
#'     gorillas_sf$mesh,
#'     prior.sigma = c(1, 0.01),
#'     prior.range = c(0.1, 0.01)
#'   )
#'   fit2 <- lgcp(
#'     geometry ~ Intercept(1) + field(geometry, model = pcmatern),
#'     data = nests,
#'     samplers = cvpart,
#'     domain = list(geometry = gorillas_sf$mesh),
#'     control.gcpo = list(enable = TRUE, type.cv = "joint")
#'   )
#'   result <- bru_block_gcpo(list(Model1 = fit1, Model2 = fit2))
#'   str(result)
#' }
#' }
#' @seealso [bru()], [bru_obs()], [bru_obs_control_gcpo()], [bru_gcpo_table()]
#' @export
bru_block_gcpo <- function(fit, ...) {
  fits <- if (is.list(fit) && !inherits(fit, "bru")) {
    fit
  } else {
    c(list(fit), list(...))
  }
  fits <- Filter(Negate(is.null), fits)

  if (length(fits) == 1L) {
    return(bru_block_gcpo_single(fits[[1]]))
  }

  results <- lapply(fits, bru_block_gcpo_single)
  if (!is.null(names(fits))) {
    names(results) <- names(fits)
  }
  results
}


#' Compare block-averaged GCPO scores across multiple fitted bru models
#'
#' @description
#' Calls [bru_block_gcpo()] on each fit and combines the results into a
#' `data.frame` with one row per block and one column per model, making it
#' straightforward to compare GCPO scores across models block by block.
#'
#' @param fits A named list of fitted `bru` objects, or `NULL` if models are
#'   passed via `\dots`.
#' @param \dots Named fitted `bru` objects, used when `fits` is not a list.
#'
#' @return
#' For a single-likelihood model, a `data.frame` with columns:
#' \describe{
#'   \item{`block`}{Block labels from `response_data$BRU_block` (for `"cp"`
#'     models) or observation indices (for other models).}
#'   \item{one column per model}{GCPO scores on the probability scale, named
#'     after the elements of `fits` or `\dots`.}
#' }
#' For multi-likelihood models, a named list of such `data.frame`s, one per
#' likelihood.
#'
#' @examples
#' \donttest{
#' if (bru_safe_inla()) {
#'   cvpart <- cv_hex(
#'     gorillas_sf$boundary,
#'     cellsize = 0.5,
#'     n_group = 3,
#'     resolution = c(95, 80)
#'   )
#'   cvpart$block_ID <- seq_len(nrow(cvpart))
#'   cvpart$group <- NULL
#'
#'   nests <- gorillas_sf$nests
#'   a <- sf::st_intersects(nests, cvpart)
#'   nests$.block <- unlist(a)
#'
#'   fit1 <- lgcp(
#'     geometry ~ Intercept(1),
#'     data = nests,
#'     samplers = cvpart,
#'     domain = list(geometry = gorillas_sf$mesh),
#'     control.gcpo = list(enable = TRUE, type.cv = "joint")
#'   )
#'   pcmatern <- INLA::inla.spde2.pcmatern(
#'     gorillas_sf$mesh,
#'     prior.sigma = c(1, 0.01),
#'     prior.range = c(0.1, 0.01)
#'   )
#'   fit2 <- lgcp(
#'     geometry ~ Intercept(1) + field(geometry, model = pcmatern),
#'     data = nests,
#'     samplers = cvpart,
#'     domain = list(geometry = gorillas_sf$mesh),
#'     control.gcpo = list(enable = TRUE, type.cv = "joint")
#'   )
#'   gcpo_df <- bru_gcpo_table(Model1 = fit1, Model2 = fit2)
#'   names(gcpo_df)
#'   colSums(log(gcpo_df[, -1]))
#' }
#' }
#' @seealso [bru_block_gcpo()]
#' @export
bru_gcpo_table <- function(fits = NULL, ...) {
  dots <- list(...)

  if (is.null(fits)) {
    fits <- dots
  } else if (inherits(fits, "bru")) {
    fits <- c(list(fits), dots)
  } else {
    fits <- c(fits, dots)
  }

  if (length(fits) == 0L) {
    stop("No fitted bru objects supplied.")
  }
  if (is.null(names(fits)) || any(names(fits) == "")) {
    stop("All fitted bru objects must be named.")
  }
  if (length(fits) == 1L) {
    stop(
      "bru_gcpo_table() requires at least two named models for comparison. ",
      "Use bru_block_gcpo() for a single fit."
    )
  }

  results <- bru_block_gcpo(fits)

  # Detect single vs multi-likelihood from first result
  first <- results[[1]]
  nlhoods <- length(first$blocks)

  build_df <- function(lhood_idx) {
    blocks <- first$blocks[[lhood_idx]]
    df <- data.frame(block = blocks)
    for (nm in names(results)) {
      gcpo <- results[[nm]]$gcpo
      df[[nm]] <- if (nlhoods == 1L) gcpo else gcpo[[lhood_idx]]
    }
    df
  }

  if (nlhoods == 1L) {
    return(build_df(1L))
  }

  out <- lapply(seq_len(nlhoods), build_df)
  names(out) <- names(first$blocks)
  out
}
