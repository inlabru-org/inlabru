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
# User-facing:                   bru_block_gcpo()


# Benchmark results (see comment below for details):
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


#' Utility functions for bru observation model GCPO
#'
#' @description
#' Extract and combine `INLA::control.gcpo` options from `bru_obs` and
#' `bru_obs_list` objects, adjusting observation indices for multi-likelihood
#' models where the response vector is the concatenation of all likelihoods.
#'
#' @param x A `bru_obs` or `bru_obs_list` object.
#' @param \dots Further arguments passed to submethods.
#'
#' @return
#' * `bru_obs_control_gcpo()` returns a list with `INLA::control.gcpo` options,
#'   with predictor/response variable indices unified for multi-observation
#'   models. If a `bru_obs` model has a `NULL` `control.gcpo` argument, an
#'   empty list is returned for that model.
#'
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
  if (!("friends" %in% names(c.gcpo))) {
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
        index_length  = response_sizes[k],
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
  
  # Combine friends and weights across likelihoods
  c.gcpo.combined[["friends"]] <-
    do.call("c", lapply(c.gcpo, function(lh) lh[["friends"]]))
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


#' Extract block-averaged GCPO scores from a fitted bru model
#'
#' @description
#' Reads `fit$gcpo$gcpo` (the raw per-observation GCPO vector produced by
#' INLA) and averages within each block defined by the `BRU_block` column in
#' each likelihood's `response_data`, returning one score per block.
#'
#' Within each block, INLA repeats the same GCPO value for every observation
#' belonging to that leave-out group. Taking the mean collapses the repeated
#' values to one representative score per block, consistent with how
#' `fit$cpo$cpo` is reported on the probability scale.
#'
#' For multi-likelihood models the raw GCPO vector is the concatenation of all
#' likelihoods in order, so cumulative row offsets are computed automatically
#' from `nrow(lh$response_data)` for each likelihood.
#'
#' @param fit A fitted object of class `bru`, with `fit$gcpo$gcpo` non-NULL
#'   (i.e. fitted with
#'   `options = list(control.compute = list(control.gcpo = list(enable = TRUE)))`)
#'   and a `BRU_block` column present in `response_data` of each likelihood.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{`blocks`}{Named list, one element per likelihood, each a vector of
#'     the unique block labels found in that likelihood's
#'     `response_data$BRU_block`.}
#'   \item{`gcpo`}{For a single-likelihood model, a numeric vector of
#'     block-averaged GCPO scores on the probability scale (one per unique
#'     block, in the order they appear in `blocks[[1]]`). For a
#'     multi-likelihood model, a named list of such vectors, one per
#'     likelihood.}
#' }
#'
#' @seealso [bru_obs()], [bru()]
#' @export
bru_block_gcpo <- function(fit) {
  fit <- bru_check_object_bru(fit)
  
  gcpo_vec <- fit[["gcpo"]][["gcpo"]]
  if (is.null(gcpo_vec)) {
    stop(
      "fit$gcpo$gcpo is NULL. Refit with ",
      "options = list(control.compute = list(control.gcpo = list(enable = TRUE)))."
    )
  }
  
  lhoods  <- as_bru_obs_list(fit)
  nlhoods <- length(lhoods)
  
  # Unique block labels per likelihood
  block_per_lhood <- lapply(lhoods, function(lh) {
    blk <- lh[["response_data"]][["BRU_block"]]
    if (is.null(blk)) {
      stop(
        "response_data$BRU_block is NULL for likelihood '",
        lh[["tag"]], "'. ",
        "A BRU_block column is required for block GCPO extraction."
      )
    }
    unique(blk)
  })
  names(block_per_lhood) <- names(lhoods)
  
  # Cumulative row offsets into the stacked GCPO vector
  n_per_lhood  <- vapply(lhoods, function(lh) nrow(lh[["response_data"]]), 0L)
  index_offset <- c(0L, cumsum(n_per_lhood))
  
  # Average gcpo_vec over rows belonging to each block in one likelihood.
  # Within a block all values are identical (INLA repeats the group score);
  # mean() is used for robustness against floating point noise.
  average_blocks <- function(lhood_idx) {
    lh      <- lhoods[[lhood_idx]]
    blk_var <- lh[["response_data"]][["BRU_block"]]
    offset  <- index_offset[lhood_idx]
    vapply(
      block_per_lhood[[lhood_idx]],
      function(b) mean(gcpo_vec[which(blk_var == b) + offset]),
      numeric(1)
    )
  }
  
  if (nlhoods == 1L) {
    return(list(
      blocks = block_per_lhood,
      gcpo   = average_blocks(1L)
    ))
  }
  
  gcpo <- lapply(seq_len(nlhoods), average_blocks)
  names(gcpo) <- names(lhoods)
  
  list(blocks = block_per_lhood, gcpo = gcpo)
}